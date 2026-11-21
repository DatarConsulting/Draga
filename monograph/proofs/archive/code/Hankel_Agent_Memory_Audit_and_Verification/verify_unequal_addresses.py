"""Exact finite tests of the audit's rational stopping-address construction.

The base codec counts letter-composition profiles and stores one integer ID.
The coupled codec adds the exact discrete mode at the oldest window endpoint.
It does not cache a tree, transition table, or subtree counts. The separate
small-tree oracle below is TEST-ONLY. These tests do not establish novelty,
asymptotic resource constants, neural accuracy, or a general formal proof.
"""
from fractions import Fraction as F
from itertools import product
from math import comb, log2, sqrt
from pathlib import Path
import json


def multinomial(ns):
    result, used = 1, 0
    for n in ns:
        result *= comb(used + n, n)
        used += n
    return result


def profiles(maxima):
    """O(d) counters; no pools of O(p) integers for Cartesian products."""
    ns = [0] * len(maxima)
    while True:
        yield tuple(ns)
        for j in range(len(ns) - 1, -1, -1):
            if ns[j] < maxima[j]:
                ns[j] += 1
                break
            ns[j] = 0
        else:
            return


class AddressCodec:
    def __init__(self, ratios, threshold):
        self.ratios = tuple(map(F, ratios))
        self.threshold = F(threshold)
        assert 0 < self.threshold < 1
        assert all(0 < r < 1 for r in self.ratios)

    def size(self, prefix_product=F(1)):
        """All nodes in the stopped subtree, counting its root."""
        if prefix_product <= self.threshold:
            return 1
        maxima = []
        for r in self.ratios:
            n, value = 0, prefix_product
            while value * r > self.threshold:
                n, value = n + 1, value * r
            maxima.append(n)
        internal = 0
        for ns in profiles(maxima):
            value = prefix_product
            for r, n in zip(self.ratios, ns):
                value *= r ** n
            if value > self.threshold:
                internal += multinomial(ns)
        return 1 + len(self.ratios) * internal

    def trim(self, newest_first):
        value, out = F(1), []
        for a in newest_first:
            out.append(a)
            value *= self.ratios[a]
            if value <= self.threshold:
                break
        return tuple(out)

    def rank(self, newest_first):
        index, value = 0, F(1)
        for a in newest_first:
            assert value > self.threshold
            index += 1
            for b in range(a):
                index += self.size(value * self.ratios[b])
            value *= self.ratios[a]
        return index

    def unrank(self, index):
        assert 0 <= index < self.size()
        value, out = F(1), []
        while index:
            assert value > self.threshold
            index -= 1
            for a, r in enumerate(self.ratios):
                count = self.size(value * r)
                if index < count:
                    out.append(a)
                    value *= r
                    break
                index -= count
            else:
                raise AssertionError('Invalid ordinal')
        return tuple(out)

    def update(self, index, letter):
        return self.rank(self.trim((letter,) + self.unrank(index)))


class CheckpointCodec:
    """Production representation: (base_mode, address_ordinal), with no log.

    mode_step depends ONLY on the old mode and letter. continuous_step uses
    the OLD mode as its second argument. Reference histories appear in tests
    only; neither update nor replay can access them.
    """
    def __init__(self, address, mode_step, continuous_step, x0, r0):
        self.address = address
        self.mode_step = mode_step
        self.continuous_step = continuous_step
        self.x0, self.r0 = x0, r0

    def initial(self):
        return self.r0, 0

    def update(self, state, letter):
        base, index = state
        word = (letter,) + self.address.unrank(index)
        retained = self.address.trim(word)
        for dropped in reversed(word[len(retained):]):
            base = self.mode_step(base, dropped)
        return base, self.address.rank(retained)

    def replay(self, state):
        mode, index = state
        value = self.x0
        for letter in reversed(self.address.unrank(index)):
            value, mode = (self.continuous_step(value, mode, letter),
                           self.mode_step(mode, letter))
        return mode, value


def verify_checkpoints():
    """Exhaustive comparison against full exact history, including resets."""
    ratios = (F(1,4), F(1,2))

    def mode_step(mode, letter):
        return 0 if letter == 0 else (mode + 1) % 3

    def continuous_step(value, mode, letter):
        # X=[0,1] is invariant: convex combination of x and (letter+mode)/4.
        r = ratios[letter]
        return r*value + (1-r)*F(letter+mode, 4)

    assert mode_step(0, 0) == mode_step(2, 0)  # Not invertible.
    for mode, letter in product(range(3), range(2)):
        assert 0 <= continuous_step(F(0), mode, letter)
        assert continuous_step(F(1), mode, letter) <= 1

    histories = updates = deletions = resets = multiple_drops = 0
    exact_startup_checks = 0
    mode_sensitive_witness = None
    for p in range(1, 5):
        address = AddressCodec(ratios, F(1, 2**p))
        codec = CheckpointCodec(address, mode_step, continuous_step, F(1,3), 2)
        for length in range(8):
            for word in product(range(2), repeat=length):
                state = codec.initial()
                actual, mode, modes = F(1,3), 2, [2]
                assert codec.replay(state) == (mode, actual)
                for t, letter in enumerate(word, 1):
                    old_retained = address.unrank(state[1])
                    actual, mode = (continuous_step(actual, mode, letter),
                                    mode_step(mode, letter))
                    modes.append(mode)
                    state = codec.update(state, letter)
                    retained = address.unrank(state[1])
                    assert retained == address.trim(tuple(reversed(word[:t])))
                    assert state[0] == modes[t-len(retained)]
                    replay_mode, rep = codec.replay(state)
                    assert replay_mode == mode
                    contraction = F(1)
                    for a in retained:
                        contraction *= ratios[a]
                    radius = contraction if contraction <= address.threshold else F(0)
                    assert abs(actual-rep) <= radius
                    if len(retained) == t:
                        assert actual == rep
                        exact_startup_checks += 1
                    else:
                        assert contraction <= address.threshold
                    dropped = 1 + len(old_retained) - len(retained)
                    deletions += dropped > 0
                    multiple_drops += dropped > 1
                    resets += letter == 0
                    if state[0] != mode and mode_sensitive_witness is None:
                        wrong_mode, wrong_value = mode, codec.x0
                        for a in reversed(retained):
                            wrong_value, wrong_mode = (
                                continuous_step(wrong_value, wrong_mode, a),
                                mode_step(wrong_mode, a))
                        if wrong_value != rep or wrong_mode != mode:
                            mode_sensitive_witness = dict(
                                precision=p, history=list(word[:t]),
                                base_mode=state[0], current_mode=mode,
                                correct_replay=str(rep), wrong_replay=str(wrong_value))
                    updates += 1
                histories += 1
    assert deletions > 0 and resets > 0 and multiple_drops > 0
    assert mode_sensitive_witness is not None
    return dict(all_passed=True, threshold_cases=4, histories=histories,
                exact_updates=updates, updates_dropping_letters=deletions,
                updates_dropping_multiple_letters=multiple_drops,
                reset_events=resets, exact_startup_checks=exact_startup_checks,
                current_mode_is_not_base_witness=mode_sensitive_witness)


def verify_guards_and_persistence():
    """Named exact regression examples for distinct invalid implications."""
    cases = []

    radius, numerical_error = F(4,256), F(2,256)
    exact_score, computed_score = F(3,256), F(5,256)
    assert abs(computed_score-exact_score) <= numerical_error
    assert abs(computed_score) > radius  # Incorrect test would accept.
    assert not abs(computed_score) > radius+numerical_error
    assert abs(exact_score)-radius == -F(1,256)
    assert exact_score-radius < 0 < exact_score+radius
    cases.append('numerical_score_without_error_budget_is_unsound')

    radius = F(1,8)
    center_score, lipschitz = radius**2, 2*radius
    exact_range = (radius**2, 2*radius**2)
    assert center_score-lipschitz*radius < 0 < exact_range[0]
    for numerator in range(-32, 33):
        x = radius*F(numerator, 32)
        assert exact_range[0] <= radius**2+x*x <= exact_range[1]
    cases.append('failed_lipschitz_test_does_not_imply_sign_ambiguity')

    # Check acceptance for positive margin and rejection for negative margin
    # at both extreme allowable score errors; also exhibit the boundary trap.
    A, margin = F(1,4), F(1,16)
    positive_exits = negative_exits = 0
    for q in range(1, 41):
        delta = F(1, 2**q)
        boundary_score = A+delta/2
        assert not abs(boundary_score) > A+delta
        assert not abs(boundary_score)+delta <= A
        for sign in (-1, 1):
            positive_score = A+margin+sign*delta
            negative_score = A-margin+sign*delta
            if 2*delta < margin:
                assert abs(positive_score) > A+delta
                assert abs(negative_score)+delta <= A
                positive_exits += 1
                negative_exits += 1
    assert positive_exits and negative_exits
    cases.append('zero_margin_can_prevent_both_finite_time_exits')

    # Genuine collision in the FULL stored state (one mode, one address).
    address = AddressCodec((F(1,2), F(1,2)), F(1,2))
    witness = []
    for word in ((0,0), (1,0)):
        value, index = F(0), 0
        for letter in word:
            value = value/2 + F(letter, 2)
            index = address.update(index, letter)
        witness.append(((0,index), value))
    assert witness[0][0] == witness[1][0]
    assert witness[0][1]-F(1,8) < 0 < witness[1][1]-F(1,8)
    cases.append('admissible_histories_with_same_memory_and_opposite_guards')

    value = F(1)
    for t in range(101):
        assert value == F(1,t+1)
        value = value/(1+value)
    cases.append('unit_asymptotic_rate_with_subexponential_forgetting')

    def multiply(a, b):
        return tuple(tuple(sum(a[i][k]*b[k][j] for k in range(2))
                           for j in range(2)) for i in range(2))

    P0, P1 = ((1,0),(0,0)), ((0,0),(0,1))
    zero = ((0,0),(0,0))
    assert multiply(P0,P0) == P0 and multiply(P1,P1) == P1
    assert multiply(P0,P1) == multiply(P1,P0) == zero
    cases.append('commuting_projectors_jsr_one_but_mixed_product_zero')
    return dict(all_passed=True, groups=len(cases), cases=cases,
                zero_margin_precision_checks=40,
                positive_margin_extreme_error_checks=positive_exits,
                negative_margin_extreme_error_checks=negative_exits,
                nonlinear_exact_iterates=101,
                scope='Finite exact examples; limits and the general guard lemma are proved in LaTeX.')


def oracle_nodes(codec):
    """Exponential oracle, only for bounded finite verification."""
    nodes = []
    def visit(word, value):
        nodes.append(word)
        if value > codec.threshold:
            for a, r in enumerate(codec.ratios):
                visit(word + (a,), value * r)
    visit((), F(1))
    return nodes


def main():
    records, state_count, transitions = [], 0, 0
    families = [((F(1,4), F(1,2)), 8),
                ((F(2,5), F(1,3)), 6),
                ((F(1,5), F(1,3), F(1,4)), 5)]
    for ratios, max_p in families:
        for p in range(1, max_p + 1):
            codec = AddressCodec(ratios, F(1, 2**p))
            oracle = oracle_nodes(codec)
            assert codec.size() == len(oracle)
            for i, word in enumerate(oracle):
                assert codec.rank(word) == i
                assert codec.unrank(i) == word
                for a in range(len(ratios)):
                    j = codec.update(i, a)
                    assert oracle[j] == codec.trim((a,) + word)
                    transitions += 1
            state_count += len(oracle)
            records.append(dict(ratios=list(map(str, ratios)), threshold=str(codec.threshold),
                                states=len(oracle), state_bits=(len(oracle)-1).bit_length()))

    # Exact trajectories for x -> x/4 and x -> x/2+1/2.
    # x0=1/4 belongs to the attractor and fixes neither branch.
    histories, updates = 0, 0
    for p in range(1, 7):
        codec = AddressCodec((F(1,4), F(1,2)), F(1, 2**p))
        for length in range(8):
            for word in product(range(2), repeat=length):
                actual, index = F(1,4), 0
                for a in word:
                    actual = actual / (4 if a == 0 else 2) + (0 if a == 0 else F(1,2))
                    index = codec.update(index, a)
                    retained = codec.unrank(index)
                    rep = F(1,4)
                    for b in reversed(retained):
                        rep = rep / (4 if b == 0 else 2) + (0 if b == 0 else F(1,2))
                    assert abs(actual-rep) <= codec.threshold
                    updates += 1
                retained = codec.unrank(index)
                assert retained == codec.trim(tuple(reversed(word)))
                rep = F(1,4)
                for b in reversed(retained):
                    rep = rep / (4 if b == 0 else 2) + (0 if b == 0 else F(1,2))
                assert abs(actual-rep) <= codec.threshold
                if len(retained) == len(word):
                    assert actual == rep
                histories += 1

    # Exact state counts at larger precisions, independently via Fibonacci.
    fib = [0, 1]
    for _ in range(130):
        fib.append(fib[-1] + fib[-2])
    examples = []
    for p in (8,16,32,64,128):
        count = 2*fib[p+2]-1
        examples.append(dict(threshold_exponent=p, states=count,
                             state_bits=(count-1).bit_length(),
                             fixed_window_startup_bits=p+1))
        if p <= 32:
            assert AddressCodec((F(1,4),F(1,2)),F(1,2**p)).size() == count

    result = dict(all_passed=True, finite_cases=len(records), states=state_count,
                  exact_transitions=transitions, histories=histories,
                  trajectory_updates=updates, cases=records, bit_examples=examples,
                  similarity_dimension=log2((1+sqrt(5))/2),
                  checkpoint_replay=verify_checkpoints(),
                  guard_and_persistence=verify_guards_and_persistence(),
                  scope='Exact finite verification and a reference implementation; no LLM or latency benchmark.')
    path = Path(__file__).with_name('unequal_address_results.json')
    path.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('cases','bit_examples')}))


if __name__ == '__main__':
    main()
