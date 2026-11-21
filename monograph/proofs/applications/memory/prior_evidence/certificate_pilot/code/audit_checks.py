"""Independent integrity audit, separate from model training and certification.

Run before training, then with --require-complete after all 50 runs finish.
Checks architecture, causal execution, labels, split boundaries, state accounting,
parameter matching, pairing and validation-only checkpoint selection. It is not
a proof of statistical generalization or of floating-point error bounds.
"""
import argparse
import copy
import hashlib
import itertools
import json
from pathlib import Path

import numpy as np
import torch

import benchmark as bm

SEEDS = [101, 102, 103, 104, 105]


def manual_s3(word):
    state = [0, 1, 2]
    for token in word:
        if token == 1:
            state[0], state[1] = state[1], state[0]
        elif token == 2:
            state[1], state[2] = state[2], state[1]
    return tuple(state)


def paired_inputs(length, seed, count=512):
    first = np.random.default_rng(seed).integers(0, 3, (count, length), dtype=np.int64)
    first[:, 8:10] = [1, 2]
    second = first.copy()
    second[:, 8:10] = [2, 1]
    return first, second


def independent_state_replay(model, x, injection, routes, rank=8):
    """Separate NumPy/float64 loop, selecting hard modes directly."""
    with torch.no_grad():
        main_input = torch.tanh(model.main_input(model.embedding(x))).double().numpy()
        decay = (.999 * torch.sigmoid(model.decay_raw)).double().numpy()
        if model.has_residual:
            matrices = model.matrices().double().numpy()[:, :, :rank, :rank]
            write = model.write[:, :rank].double().numpy()
            u = injection.double().numpy()
            modes = routes.argmax(-1).numpy()
    batch, length = x.shape
    main = np.zeros((batch, model.main_dim))
    residual = np.zeros((batch, model.heads, rank))
    outputs = []
    for t in range(length):
        main = main * decay + main_input[:, t]
        if model.has_residual:
            for sample in range(batch):
                for head in range(model.heads):
                    residual[sample, head] = (
                        matrices[modes[sample, t], head] @ residual[sample, head]
                        + write[head] * u[sample, t, head]
                    )
            outputs.append(np.concatenate([main, residual.reshape(batch, -1)], axis=1))
        else:
            outputs.append(main.copy())
    return np.stack(outputs, axis=1)


def run(require_complete=False):
    torch.set_num_threads(1)
    root = Path(__file__).resolve().parents[1]
    report = {"checks": {}, "limitations": [
        "Integrity checks are distinct from the separately generated exact stored-coefficient certificate.",
        "256 bytes counts recurrent FP32 coordinates; weights, activations, buffers and framework overhead are excluded.",
        "Full-sequence forward retains intermediate states and is not an optimized streaming kernel.",
        "The stacked Euclidean contraction is a sufficient, stronger constraint than general mean-square stability.",
        "A negative outcome within this restricted parameterization does not refute every certificate-compatible architecture.",
        "No exact optimal-rank AAK solver or formal floating-point execution proof is performed by this audit.",
    ]}
    checks = report["checks"]

    words = np.array(list(itertools.product(range(3), repeat=7)), dtype=np.int64)
    expected = np.array([bm.LOOKUP[manual_s3(row)] for row in words])
    assert np.array_equal(bm.labels_s3(words)[:, -1], expected)
    assert len(set(expected)) == 6
    checks["s3_independent_exhaustive_length7"] = {"pass": True, "words": len(words)}

    paired_checks = {}
    for length, seed in [(32, 19001), (96, 19002)]:
        first, second = paired_inputs(length, seed)
        assert all(np.array_equal((first == a).sum(1), (second == a).sum(1)) for a in range(3))
        assert all(manual_s3(x) != manual_s3(y) for x, y in zip(first, second))
        assert np.array_equal(first[:, :8], second[:, :8])
        assert np.array_equal(first[:, 10:], second[:, 10:])
        paired_checks[str(length)] = {"pass": True, "pairs": len(first),
                                      "same_symbol_counts": True, "all_labels_differ": True,
                                      "common_suffix_length": length - 10}
    checks["order_counterfactual_pairs"] = paired_checks

    validation = np.random.default_rng(17001).integers(0, 3, (512, 32), dtype=np.int64)
    test = np.random.default_rng(18032).integers(0, 3, (1024, 32), dtype=np.int64)
    first, second = paired_inputs(32, 19001)
    held_out = set(row.astype(np.uint8).tobytes() for row in np.concatenate([validation, test, first, second]))
    overlap = {}
    for seed in SEEDS:
        rng = np.random.default_rng(10000 + seed)
        matches = 0
        digest = hashlib.sha256()
        for _ in range(1500):
            batch = rng.integers(0, 3, (64, 32), dtype=np.int64)
            digest.update(batch.astype(np.uint8).tobytes())
            matches += sum(row.astype(np.uint8).tobytes() in held_out for row in batch)
        assert matches == 0
        overlap[str(seed)] = {"pass": True, "train_sequences": 96000,
                              "exact_held_out_matches": matches,
                              "common_training_stream_sha256": digest.hexdigest()}
    checks["order_sequence_separation_and_common_streams"] = overlap

    encoded, vocab = bm.load_text()
    train_text = (root / "data/train.txt").read_text()
    assert vocab == ["<UNK>"] + sorted(set(train_text))
    text_checks = {}
    for split, data in encoded.items():
        raw = (root / "data" / f"{split}.txt").read_bytes()
        starts = np.arange(0, len(data) - 64, 64)
        targets = starts[:, None] + np.arange(1, 65)[None, :]
        assert targets.max() < len(data)
        assert len(np.unique(targets)) == targets.size
        x, y = bm.text_windows(data, starts[:2], 64)
        assert torch.equal(x[:, 1:], y[:, :-1])
        text_checks[split] = {"pass": True, "characters": len(data),
                             "sha256": hashlib.sha256(raw).hexdigest(),
                             "evaluation_targets": int(targets.size),
                             "unknown_characters": int((data == 0).sum())}
    checks["text_train_only_vocabulary_and_split_bounds"] = text_checks
    train_blocks = {train_text[i:i + 64] for i in range(len(train_text) - 63)}
    duplicate_counts = {}
    for split in ["validation", "test"]:
        held = (root / "data" / f"{split}.txt").read_text()
        duplicate_counts[split] = sum(held[i:i + 64] in train_blocks for i in range(len(held) - 63))
    checks["text_exact_64char_overlap_with_train"] = {
        "measured": True, "counts": duplicate_counts,
        "interpretation": "Reported without changing splits; repeated literary strings are not automatically leakage.",
    }

    architecture, parameter_matching, stability_isolation = {}, {}, {}
    for task, size, outputs in [("order", 3, 6), ("text", len(vocab), len(vocab))]:
        target = bm.make_model("multi_cert", size, outputs, task, 101)
        path = bm.make_model("multi_path", size, outputs, task, 101)
        assert target.state_dict().keys() == path.state_dict().keys()
        assert all(torch.equal(target.state_dict()[key], path.state_dict()[key]) for key in target.state_dict())
        stability_isolation[task] = {"pass": True, "identical_initial_raw_parameters": True,
                                     "difference": "Only transition normalization differs."}
        target_count = bm.count_params(target)
        for arm in bm.ARMS:
            model = bm.make_model(arm, size, outputs, task, 101).eval()
            torch.manual_seed(912)
            x = torch.randint(size, (3, 17))
            with torch.no_grad():
                full, main, residual, states, injection, routes = model(x, return_parts=True)
                prefix = model(x[:, :9])
            causal_error = float((full[:, :9] - prefix).abs().max())
            assert causal_error < 2e-6
            state_count = model.main_dim + (model.heads * model.rank if model.has_residual else 0)
            assert state_count == 64 and states.shape[-1] == 64
            manual_states = independent_state_replay(model, x, injection, routes)
            replay_error = float(np.abs(manual_states - states.numpy()).max())
            assert replay_error < 2e-5
            assert torch.allclose(full, main + residual, rtol=0, atol=1e-6)
            row = {"pass": True, "prefix_max_abs_error": causal_error,
                   "independent_float64_state_replay_max_error": replay_error,
                   "state_scalars": state_count, "state_bytes_fp32": state_count * 4,
                   "parameters": bm.count_params(model)}
            if model.has_residual:
                assert injection.abs().max() <= 1
                assert torch.equal(routes.sum(-1), torch.ones_like(routes[..., 0]))
                assert set(routes.unique().tolist()).issubset({0., 1.})
                with torch.no_grad():
                    residual_state = states[:, :, model.main_dim:].reshape(3, 17, 4, 8)
                    scalar = (residual_state * model.head_read).sum(-1)
                    independent_residual = scalar @ model.output_vectors.weight.T
                    scalar_error = float((independent_residual - residual).abs().max())
                    assert scalar_error < 2e-5
                    projected = model(x, project_rank=4, return_parts=True)
                    assert projected[3].shape[-1] == 48
                    assert torch.equal(projected[1], main)
                    assert torch.equal(projected[4], injection)
                    assert torch.equal(projected[5], routes)
                    mutation = copy.deepcopy(model)
                    mutation.write.add_(1.0)
                    mutated = mutation(x, return_parts=True)
                    assert torch.equal(mutated[1], main)
                    assert torch.equal(mutated[4], injection)
                    assert torch.equal(mutated[5], routes)
                    mat = model.matrices().double().numpy()
                row.update(scalar_additive_readout_max_error=scalar_error,
                           projection_preserves_main_routes_and_injections=True,
                           residual_mutation_preserves_main_routes_and_injections=True,
                           projected_state_scalars=48, projected_state_bytes_fp32=192)
                if arm == "multi_cert":
                    head_max = []
                    projected_max = []
                    for h in range(4):
                        gram = sum(a.T @ a for a in mat[:, h])
                        projected_gram = sum(a[:4, :4].T @ a[:4, :4] for a in mat[:, h])
                        head_max.append(float(np.linalg.eigvalsh(gram).max()))
                        projected_max.append(float(np.linalg.eigvalsh(projected_gram).max()))
                    assert max(head_max) < .99 and max(projected_max) < .99
                    row.update(initial_stack_gram_largest_eigenvalues=head_max,
                               initial_projected_gram_largest_eigenvalues=projected_max,
                               exact_stability_check="Performed separately by certify.py on selected effective matrices.")
            architecture[f"{task}_{arm}"] = row
            if arm == "single_param":
                difference = bm.count_params(model) - target_count
                assert abs(difference) <= (33 + outputs) / 2
                parameter_matching[task] = {"pass": True, "target": target_count,
                                            "control": bm.count_params(model), "difference": difference}
    checks["architecture_and_causal_replay"] = architecture
    checks["parameter_control_nearest_width"] = parameter_matching
    checks["path_vs_counting_stability_isolation"] = stability_isolation

    present, missing = [], []
    for task in ["text", "order"]:
        for arm in bm.ARMS:
            for seed in SEEDS:
                stem = f"{task}_{arm}_seed{seed}"
                path = root / "results" / f"{stem}.json"
                checkpoint = root / "checkpoints" / f"{stem}.pt"
                if not path.exists():
                    missing.append(stem)
                    continue
                result = json.loads(path.read_text())
                key = "bpc" if task == "text" else "loss"
                chosen = min(result["trace"], key=lambda item: item["validation"][key])
                assert result["selected_step"] == chosen["step"]
                assert result["validation_best"] == chosen["validation"][key]
                assert result["steps"] == (1800 if task == "text" else 1500)
                assert result["state_bytes_fp32"] == 256
                assert checkpoint.exists()
                present.append(stem)
    if require_complete:
        assert not missing, missing
    checks["validation_only_checkpoint_selection"] = {"pass": True, "present_runs": len(present),
                                                      "missing_runs": missing,
                                                      "completeness_required": require_complete}
    manifest_path = root / "results/run_manifest.json"
    if manifest_path.exists():
        manifest = json.loads(manifest_path.read_text())
        for filename, digest in manifest["hashes"].items():
            assert hashlib.sha256((root / filename).read_bytes()).hexdigest() == digest, filename
        checks["protocol_code_and_data_manifest_unchanged_since_launch"] = {
            "pass": True, "verified_sha256": manifest["hashes"]}
    report["source_sha256"] = {name: hashlib.sha256((root / name).read_bytes()).hexdigest()
                                for name in ["PROTOCOL.md", "code/benchmark.py", "code/audit_checks.py"]
                                if (root / name).exists()}
    report["all_checks_passed"] = True
    (root / "results").mkdir(exist_ok=True)
    dest = root / "results/audit.json"
    dest.write_text(json.dumps(report, indent=2))
    print(json.dumps({"all_checks_passed": True, "present_runs": len(present), "report": str(dest)}))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--require-complete", action="store_true")
    run(parser.parse_args().require_complete)
