"""Independent verifier of serialized rational certificates.

Does not import certify.py, its Dyadic class, its Stein solver or its LDL code.
Recomputes all decisive statements from Fraction-valued matrices and binds them
to the checkpoint bytes. No floating-point estimate is used to accept a bound.
"""
import argparse
import copy
import hashlib
import json
from fractions import Fraction
from pathlib import Path

import numpy as np
import torch

ROOT = Path(__file__).resolve().parents[1]
F = Fraction


def rational(obj):
    return F(int(obj['numerator']), int(obj['denominator']))


def unpack(obj):
    shape = tuple(obj['shape'])
    raw = np.asarray(obj['numerators'], dtype=object)
    assert raw.shape == shape
    denominator = 2 ** int(obj['denominator_power_of_two'])
    return np.asarray([F(int(x), denominator) for x in raw.flat], dtype=object).reshape(shape)


def exact_float_array(array):
    array = np.asarray(array)
    return np.asarray([F.from_float(float(x)) for x in array.flat], dtype=object).reshape(array.shape)


def zero(n, m=None):
    return np.full((n, n if m is None else m), F(0), dtype=object)


def quadratic(x, matrix):
    return (x.T @ matrix @ x).item()


def schur_positive(matrix):
    """Exact symmetric Schur elimination, separately written from LDL producer."""
    assert np.array_equal(matrix, matrix.T)
    current = matrix.copy()
    pivots = []
    while len(current):
        pivot = current[0, 0]
        assert pivot > 0, 'Nonpositive exact Schur pivot'
        pivots.append(pivot)
        column = current[1:, :1]
        current = current[1:, 1:] - column @ column.T / pivot
    return pivots


def stability(matrices, q, serialized):
    n = matrices.shape[-1]
    gram = sum((a.T @ a for a in matrices), zero(n))
    gap = np.eye(n, dtype=object) * q - gram
    pivots = schur_positive(gap)
    assert serialized['passed'] is True
    assert rational(serialized['q']) == q
    assert pivots == [rational(p) for p in serialized['ldl_pivots']]


def apply(matrices, x):
    return sum((a @ x @ a.T for a in matrices), zero(len(x)))


def verify(cert, checkpoint_dir):
    path = checkpoint_dir / cert['checkpoint']
    assert path.is_file()
    assert hashlib.sha256(path.read_bytes()).hexdigest() == cert['checkpoint_sha256']
    ck = torch.load(path, map_location='cpu', weights_only=True)
    assert ck['args']['arm'] == cert['arm'] == 'multi_cert'
    assert ck['args']['task'] == cert['task'] and ck['args']['seed'] == cert['seed']
    exported = cert['effective_coefficients_exact_dyadic']
    r = unpack(exported['R_modes_heads'])
    b = unpack(exported['b_heads'])
    c = unpack(exported['c_heads'])
    v = unpack(exported['v_output_by_head'])
    source = [ck['effective_matrices'], ck['state_dict']['write'],
              ck['state_dict']['head_read'], ck['state_dict']['output_vectors.weight']]
    for name, exact, tensor in zip(['R', 'b', 'c', 'v'], [r, b, c, v], source):
        assert tensor.dtype == torch.float32
        assert np.array_equal(exact, exact_float_array(tensor.numpy()))
        assert np.array_equal(exact, exact_float_array(cert['effective_coefficients_decimal'][name]))
    q = rational(cert['q'])
    assert q == F(99, 100)
    assert cert['projection_rank_per_head'] == 4
    assert cert['source_state_bytes'] == 256 and cert['projected_state_bytes'] == 192
    assert b.shape == c.shape == (4, 8)
    assert len(cert['heads']) == 4
    aggregate = F(0)
    for h, item in enumerate(cert['heads']):
        assert item['head'] == h
        n, k = item['original_rank'], item['reduced_rank']
        assert (n, k) == (8, 4)
        rh = r[:, h]
        stability(rh, q, item['full_stability'])
        stability(rh[:, :k, :k], q, item['projected_stability'])
        e = np.full((len(rh), n + k, n + k), F(0), dtype=object)
        e[:, :n, :n], e[:, n:, n:] = rh, rh[:, :k, :k]
        B = np.concatenate([b[h], b[h, :k]])[:, None]
        C = np.concatenate([c[h], -c[h, :k]])[:, None]
        g = unpack(item['stein_candidates']['G0'])
        s = unpack(item['stein_candidates']['S0'])
        assert g.shape == s.shape == (n + k, n + k)
        assert np.array_equal(g, g.T) and np.array_equal(s, s.T)
        d1 = B @ B.T + apply(e, g) - g
        d2 = g + apply(e, s) - s
        a1 = sum(abs(x) for x in d1.flat)
        a2 = sum(abs(x) for x in d2.flat)
        cnorm = (C.T @ C).item()
        bnorm = (B.T @ B).item()
        center = quadratic(C, s)
        correction = cnorm * (a1 / (1 - q) ** 2 + a2 / (1 - q))
        upper = center + correction
        generic = cnorm * bnorm / (1 - q) ** 2
        used = min(upper, generic)
        assert used >= 0
        recomputed = {
            'exact_Stein_candidate_quadratic': center,
            'exact_residual1_entrywise_l1': a1,
            'exact_residual2_entrywise_l1': a2,
            'exact_error_correction': correction,
            'exact_Stein_upper_squared': upper,
            'exact_generic_upper_squared': generic,
            'selected_upper_squared': used,
        }
        for key, value in recomputed.items():
            assert value == rational(item[key]), key
        eta = rational(item['eta_rational'])
        assert eta >= 0 and eta ** 2 >= used
        assert F.from_float(item['eta_float']) >= eta
        assert item['eta_squared_ge_upper_exactly'] is True
        vinf = max(abs(x) for x in v[:, h])
        assert vinf == rational(item['v_infinity_exact'])
        aggregate += vinf * eta
    assert aggregate == rational(cert['uniform_logit_coefficient_A_exact'])
    for length in [32, 64, 96, 128]:
        key = str(length)
        delta = rational(cert['uniform_logit_bounds_exact_upward'][key])
        tv = rational(cert['uniform_TV_bounds_exact_upward'][key])
        assert delta >= 0 and delta ** 2 >= aggregate ** 2 * length
        assert tv == min(F(1), delta / 2)
        assert F.from_float(cert['uniform_logit_bounds_float'][key]) >= delta
        assert F.from_float(cert['uniform_TV_bounds_float'][key]) >= tv
    assert cert['all_stability_checks_passed'] is True
    assert cert['all_error_bounds_exactly_rounded_up'] is True
    return {'checkpoint': cert['checkpoint'], 'heads_verified': 4,
            'effective_parameters_match_checkpoint_exactly': True,
            'exact_stability_verified': True, 'exact_error_upper_bounds_verified': True,
            'exact_logit_coefficient_verified': True,
            'uniform_logit_and_TV_bounds_rounded_up_verified': True}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--require-complete', action='store_true')
    args = parser.parse_args()
    paths = sorted((ROOT / 'certificates').glob('*_multi_cert_seed*.json'))
    if args.require_complete:
        assert len(paths) == 10, len(paths)
    assert paths, 'No serialized certificates found'
    out = []
    for path in paths:
        cert = json.loads(path.read_text())
        out.append(verify(cert, ROOT / 'checkpoints'))
        print(json.dumps({'verified': path.name}), flush=True)
    # Negative control: changing a claimed exact bound must be rejected.
    bad = copy.deepcopy(json.loads(paths[0].read_text()))
    bad['heads'][0]['selected_upper_squared']['numerator'] = '-1'
    rejected = False
    try:
        verify(bad, ROOT / 'checkpoints')
    except AssertionError:
        rejected = True
    assert rejected
    result = {'runs': out, 'all_verified': True,
              'tampered_bound_rejected': rejected,
              'complete_ten_runs_required': args.require_complete,
              'verification_arithmetic': 'Python Fraction and independent Schur elimination; no certify.py imports',
              'scope': 'Exact real-arithmetic stored-coefficient certificate; no execution-roundoff or optimality claim',
              'verifier_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    destination = ROOT / 'results/certificate_verification.json'
    destination.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({'all_verified': True, 'runs': len(out), 'result': str(destination)}))


if __name__ == '__main__':
    main()
