#!/usr/bin/env python3
import json,html,platform
from pathlib import Path
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parent
def main():
    d=pd.read_csv(ROOT/'results/sizes.csv');t=pd.read_csv(ROOT/'results/times.csv')
    d['charged_bytes_per_head']=d.packet_bytes+d.basis_bytes/16+8
    native=d[d.offset==0].copy()
    frames=d.groupby(['doc','offset','method']).agg(packet_bytes=('packet_bytes','sum'),basis_bytes=('basis_bytes','first'),heads=('head','count')).reset_index()
    assert (frames.heads==16).all()
    frames['total_bytes']=128+frames.basis_bytes+frames.packet_bytes
    frames.to_csv(ROOT/'results/layer_frames.csv',index=False)
    med=t.groupby(['doc','offset','head','method'])[['encode_us','decode_us','cycle_us']].median().reset_index()
    timing=med[med.cycle_us>=0].groupby('method')[['encode_us','decode_us','cycle_us']].median()
    def aggregate(data):
        a=data.groupby('method').agg(bytes_per_head=('charged_bytes_per_head','mean'),packet_bytes=('packet_bytes','mean'),predictor_metadata_bytes=('model_bytes','mean'),residual_field_bytes=('residual_bytes','mean'),basis_bytes_per_layer=('basis_bytes','first'),relative_l2=('relative_l2','mean'))
        a['bytes_per_layer']=a.bytes_per_head*16
        return a.join(timing)
    primary=aggregate(native);allstates=aggregate(d)
    primary.to_csv(ROOT/'results/native_summary.csv');allstates.to_csv(ROOT/'results/replay_summary.csv')
    selection=json.loads((ROOT/'bases/selection.json').read_text());sel={(r['method'],r['rank']):r for r in selection['rows']}
    comparison=[]
    for rank in [1,2,4,8,16,32,64,96]:
        H=primary.loc[f'hankel_{rank}'];E=primary.loc[f'energy_{rank}']
        comparison.append({'rank':rank,'theta':sel['hankel',rank]['theta'],
            'hankel_norm_ratio':sel['hankel',rank]['hankel_error']/sel['energy',rank]['hankel_error'],
            'hankel_bytes_per_head':float(H.bytes_per_head),'energy_bytes_per_head':float(E.bytes_per_head),
            'byte_saving_vs_energy_pct':float(100*(1-H.bytes_per_head/E.bytes_per_head))})
    comp=pd.DataFrame(comparison);comp.to_csv(ROOT/'results/matched_rank_comparison.csv',index=False)
    bydoc=native.groupby(['doc','method']).charged_bytes_per_head.mean().unstack();bydoc.to_csv(ROOT/'results/native_by_document.csv')
    h=d[d.method.str.startswith('hankel_')]
    optimistic=native[native.method.str.startswith('hankel_')].groupby(['doc','head']).packet_bytes.min().mean()+8
    checks=json.loads((ROOT/'results/checks.json').read_text());ind=json.loads((ROOT/'results/independent_checks.json').read_text());guards=json.loads((ROOT/'results/guard_checks.json').read_text())
    gt=pd.read_csv(ROOT/'results/guard_times.csv');gt=gt[gt.cycle_us>=0].groupby(['doc','offset','head','rank'])[['encode_us','decode_us','cycle_us']].median().groupby('rank').median()
    replay=pd.read_csv(ROOT/'results/replay.csv')
    report={'conclusion':'No meaningful byte advantage established for the frozen finite-Hankel/Stein predictor family. Rank96 expands; small-rank matched gains stay below 0.5%. Every Hankel packet loses to Codec A on every evaluated head, even before basis charging.',
      'scope':{'primary_native_heads':48,'documents':3,'layer':1,'additional_replayed_head_states':192,'new_full_model_inference':False,'independent_new_corpus':False},
      'primary':primary.reset_index().to_dict('records'),'matched_ranks':comparison,
      'hankel_beats_A_packet_count':int((h.packet_bytes<h.codec_A_bytes).sum()),'hankel_packet_count':len(h),
      'optimistic_rank_oracle_free_basis_bytes_per_head':float(optimistic),
      'checks':checks,'independent_checks':ind,'guard_checks':guards,'guard_timing':gt.reset_index().to_dict('records'),
      'max_replayed_vs_archived_readout_relative_error':float(replay.relative_readout_error.max()),
      'calibration_refit_seconds':selection['seconds'],
      'predeclared_meaningful_criterion_met':False}
    (ROOT/'results/summary.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    def table(df,fmt=3):return df.to_html(index=False,border=0,float_format=lambda x:f'{x:,.{fmt}f}',classes='data')
    chosen=['raw','codec_A','codec_B','local','row_ref','energy_2','hankel_2','cross_16','energy_96','hankel_96']
    main_table=primary.loc[chosen].reset_index()[['method','bytes_per_head','bytes_per_layer','cycle_us']]
    main_table.columns=['Method','Bytes/head equivalent','Bytes/layer snapshot','Median cycle, µs/head']
    rank_table=comp[['rank','theta','hankel_norm_ratio','energy_bytes_per_head','hankel_bytes_per_head','byte_saving_vs_energy_pct']].copy()
    rank_table.columns=['Rank','Selected θ','Hankel norm H/E','Energy bytes/head','Hankel bytes/head','Byte saving vs energy, %']
    bydoc_table=bydoc[['codec_A','codec_B','energy_2','hankel_2','energy_96','hankel_96']].reset_index()
    H96=primary.loc['hankel_96'].bytes_per_head;E96=primary.loc['energy_96'].bytes_per_head
    H2=primary.loc['hankel_2'].bytes_per_head;E2=primary.loc['energy_2'].bytes_per_head
    text=f'''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Does Hankel reduce exact GDN storage?</title><style>
body{{font:16px/1.6 system-ui,sans-serif;background:#edf1f5;color:#182a38;margin:0}}main{{max-width:1100px;margin:auto;background:white;padding:38px}}h1{{font-size:34px;line-height:1.2}}h2{{margin-top:36px;font-size:23px}}.meta{{color:#52697b;font-size:13px;letter-spacing:.06em}}.verdict{{padding:20px 24px;background:#fff3e4;border-left:5px solid #ac6500}}.note{{padding:16px;background:#edf5fa}}table{{border-collapse:collapse;width:100%;font-size:14px;margin:18px 0}}th,td{{padding:10px;text-align:right;border-bottom:1px solid #d5dde4}}th:first-child,td:first-child{{text-align:left}}th{{background:#edf3f7}}.scroll{{overflow-x:auto}}code,pre{{font-family:ui-monospace,monospace}}pre{{background:#f2f5f7;padding:16px;white-space:pre-wrap}}.small{{font-size:13px;color:#52697b}}@media print{{body{{background:white}}main{{padding:15px}}h2,table{{break-inside:avoid}}}}</style><main>
<p class="meta">DATAR CONSULTING · EXECUTED RESEARCH PILOT · 6 NOVEMBER 2026</p>
<h1>Does a Hankel predictor make exact GDN storage cheaper?</h1>
<p>Matched residual coding on archived Qwen3.5-0.8B Q4_0 states. The Hankel-selected predictor, energy predictor and cross-approximation predictor all use the same reversible residual format.</p>
<div class="verdict"><strong>The tested Hankel construction does not demonstrate a meaningful compression advantage.</strong>
<p>At the archived rank 96, the calibration surrogate’s Hankel error falls by {100*(1-comp.iloc[-1].hankel_norm_ratio):.2f}%, but total encoded storage rises by {100*(H96/E96-1):.3f}% versus the energy basis. The best matched-rank saving on the primary native snapshots is {comp.byte_saving_vs_energy_pct.max():.3f}%, at rank 2.</p>
<p>No Hankel candidate beats Codec A on any evaluated head, even before charging the shared basis. This rejects the declared 5% pilot criterion for this family and dataset.</p></div>
<h2>1. Primary measurements</h2><p>The primary data are <strong>three native layer snapshots, 16 heads each</strong>: atelier, coast and ledger, layer 1 after token 428. Each head is a 128 × 128 FP32 matrix. They were excluded from the two-text basis calibration, but they are historical diagnostic texts, not a fresh independent corpus.</p>
<div class="scroll">{table(main_table,2)}</div>
<p class="small">Bytes include the complete per-head packet, one FP32 basis shared by the 16 simultaneous heads where needed, and a 128-byte layer header. “Bytes/head equivalent” divides this self-contained layer snapshot by 16; no basis is amortized across tokens. Raw framed storage is 65,560 bytes/head versus 65,536 bytes of state data. Timing is a separate median over local replay cycles, not full-model inference.</p>
<h2>2. Does the Hankel objective predict encoded length?</h2><div class="scroll">{table(rank_table,4)}</div>
<p>At rank 96, both methods pay identical basis and coefficient costs. The Hankel residual itself is larger. At rank 2, the small mean gain is consistent across these three documents, but does not approach the declared 5% criterion. Rank 2 is an exploratory sweep result, not a newly validated deployment choice.</p>
<div class="scroll">{table(bydoc_table,2)}</div>
<p>An optimistic oracle that selects the shortest Hankel rank independently for every native head and gives every basis to the decoder for free still costs <strong>{optimistic:,.2f} bytes/head</strong>, versus {primary.loc['codec_A'].bytes_per_head:,.2f} for Codec A. This diagnostic is not a deployable format; it shows that rank selection and basis amortization alone cannot fix the measured deficit.</p>
<h2>3. Exactness and replay checks</h2>
<ul><li>{checks['roundtrips']:,} exact round trips across 240 head states and 27 configurations.</li>
<li>{checks['exact_transitions']:,} cycles preserve both updated FP32 state words and GDN readout words relative to the same portable reference kernel.</li>
<li>Six complete layer frames written, reloaded and reconstructed; four new-predictor frames independently parsed, totaling {ind['total_words']:,} identical FP32 words.</li>
<li>116 standalone checks, including signed zero, random IEEE words for the row-reference path, and malformed packet rejection. AddressSanitizer and UndefinedBehaviorSanitizer pass these checks; LeakSanitizer is disabled.</li></ul>
<p>Only the offset-zero matrices are native captured states. The other 192 head states come from 32-step replays using archived q/k/v/g/β controls. The portable kernel’s raw readouts differ slightly from the archived llama.cpp readouts: maximum relative L2 error {replay.relative_readout_error.max():.3g}. Codec equality is against the portable kernel itself. No fresh logits, language-quality score, GPU result or full-model latency is claimed.</p>
<h2>4. What was held fixed?</h2>
<p>The archived rank-96 Hankel basis is reused byte for byte, transposed only for storage layout. For each other rank, 41 energy/Stein interpolation bases are scored using the same complete finite-support scalar word-Hankel surrogate and calibration inputs. The chosen θ is frozen before encoding evaluation data.</p>
<p>This is the existing finite-surrogate Hankel selector: four action letters, two-action horizon, discount 19/20 and calibration queries. Its ordinary Hankel rank is at most 22. Codec value rank (1–96) is a different quantity. It is not the full Qwen sequence Hankel, and it is not a newly learned dynamic core-conditioned predictor.</p>
<p>For a shared basis B, coefficients C = BᵀS are evaluated in a fixed FP32 FMA order. Each length-128 coefficient vector is stored as 128 signed 16-bit integers and a signed power-of-two exponent. The decoder forms its prediction from exactly those stored coefficients and the FP32 basis. Cross approximation uses the same factor quantizer and FMA evaluation; its two factor vectors cost 516 bytes per component. Basis coefficients cost 258 bytes per component per head. All factors are charged.</p>
<p>Every new predictor feeds exactly the same magnitude-word residual packer: d = magnitude_bits(S) − magnitude_bits(P), signed zigzag coding, one globally selected low-bit width, patched high-bit exceptions and a separate original-sign bitmap. The encoder never substitutes a rounded floating-point subtraction for this exact residual. Codec A, adaptive Codec B and the local codec are retained as external reference formats.</p>
<pre>Residual bytes = 2,048 + 2,048 b
               + 2,048 · 1[exceptions &gt; 0]
               + ceil(exceptions · (32 − b) / 8)

Head packet = 16-byte header + predictor metadata + residual
Shared-basis layer = 128-byte header + 512r basis bytes
                   + sum of the 16 head packets</pre>
<p>The reference Codec B searches ranks 1, 2, 4, 8, 16 and 32 per head and can fall back to Codec A. Its adapted framing adds the same 16-byte outer header used here. The new fixed-rank cross predictors use the identical residual mechanism but do not adapt rank, which explains why the adaptive reference can be smaller.</p>
<h2>5. Fallback and runtime cost</h2>
<p>A guarded encoder was also run at ranks 2 and 96. It encodes both the Hankel candidate and Codec A, then stores the smaller packet. It selected Hankel <strong>{guards['selected_hankel']} times out of {guards['roundtrips']}</strong>. These checkpoints therefore have Codec A’s size and need no stored Hankel basis; the encoder still needs the basis and pays the unsuccessful predictor search.</p>
<div class="scroll">{table(gt.reset_index(),3)}</div>
<p class="small">Guard times are µs/head, measured in a separate sequential follow-up; they are descriptive, not a controlled statistical speedup claim. Main benchmark: one CPU pinned, one warmup and three repetitions, configurations rotated across records and trials. C++ with -O3 -march=native, no fast-math, implicit FP contraction disabled, explicit FMA used by predictors. Encoding includes metadata allocation; decoding and the cycle include prediction reconstruction. Calibration is an offline cost, measured at {selection['seconds']:.3f} seconds here.</p>
<h2>6. Mathematical interpretation and next decision</h2>
<p>The exact residual provides injectivity. The predictor only changes the number and distribution of residual bits. Reducing a future-readout norm does not, without further assumptions, minimize the encoded integer residual. The current experiment measures that mismatch directly.</p>
<p>For example, future probes can be insensitive to a component of state while an exact full-state codec must still retain it. Even small-amplitude FP32 values can contain many independent mantissa bits. A useful positive theorem would need explicit control of the word residual distribution, predictor description cost and small-magnitude behavior, rather than a Hankel singular-value tail alone.</p>
<p><strong>Decision:</strong> do not claim that this Hankel basis makes exact GDN storage materially cheaper. Retain this as a reproducible negative ablation. A different dynamic core or a selector trained directly on total coded length is a new candidate requiring a new frozen test; the present result neither validates nor rules out every such construction.</p>
<h2>7. Reproduction and evidence</h2>
<pre>python3 reproduce.py
# Optional: recompute all 41 calibration candidates per rank
python3 reproduce.py --refit
# Linux x86-64, g++, gcc, Python 3, NumPy, pandas; no model download.</pre>
<p>The bundle contains source code, the frozen contract, input provenance and hashes, deployed bases, all row-level measurements, independent decoder, exact sample frames and the report. Primary input artifacts: Hankel_GDN_Test_Reproduction.zip; Hankel_GDN_Test_Evidence_01_of_02.zip; Hankel_GDN_Test_Evidence_02_of_02.zip. Existing codec sources come from BTSP_Method2_GDN_Codec_Benchmark_20261105.zip. Their exact Library identifiers and file hashes are in provenance.json.</p>
<p class="small">Prepared for Vashista Nobaub, Datar Consulting, with AI assistance. These are executed local pilot measurements, not peer review, formal proof-assistant certification or a population-level performance guarantee.</p></main></html>'''
    (ROOT/'Hankel_Exact_Residual_Report.html').write_text(text)
    print(json.dumps({'native_H96':H96,'native_E96':E96,'native_H2':H2,'native_E2':E2,'max_matched_saving_percent':comp.byte_saving_vs_energy_pct.max(),'decision':'criterion not met'}))
if __name__=='__main__':main()
