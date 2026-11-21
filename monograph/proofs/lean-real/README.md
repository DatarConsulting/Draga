# Exact persistent memory for two dyadic maps

A binary stream updates a real number by $T_0(x)=x/4$ or $T_1(x)=(1+x)/2$. After every finite history, including the empty history, the machine must output a number within $2^{-n-1}$ of the true state. The update receives the next letter and the persistent state. A decoder reads this state and may output any real number. The initial point and precision are fixed when the machine is chosen.

**Theorem.** Write $U_n$ for the union of the stopped leaf intervals defined below. For each $x_0\in U_n$, the minimum is $F_{n+2}$ states and the least fixed width is $\lceil\log_2 F_{n+2}\rceil$ persistent bits. Here $F_0=0$, $F_1=1$. The state lower bound holds throughout $[0,1]$. At $n=128$ and $x_0=0$, the minimum is 659,034,621,587,630,041,982,498,215 states and 90 bits. Working memory for arithmetic and the decoder's output buffer are separate resources in this model.

## Three steps of the proof

1. **Stop at a geometric scale.** Read an address from its newest letter. Charge two units for 0 and one for 1. Stop as soon as the cumulative cost reaches $n$. For an address $w$, let $I_w=T_w([0,1])$. Its width is $2^{-n}$ or $2^{-n-1}$. Splitting at the first letter gives $F_{n+2}$ leaves. Consecutive right endpoints differ by at least $3\cdot 2^{-n}/2$. [Leaf geometry and spacing](FibonacciMemory.lean#L137)
2. **Build indistinguishable-state contradictions from reachable values.** For each leaf, construct a chronological input of length $n+3$ by padding with ones and then appending the reversed address. Its final value lies close enough to that leaf's right endpoint that any two witnesses are more than $2^{-n}$ apart. A real readout shared by two witnesses would violate the triangle inequality. Their states must therefore be distinct. This argument applies to every deterministic machine meeting the contract. [Reachable witnesses and lower bound](FibonacciMemory.lean#L302)
3. **Attain the bound.** Initialize in any leaf containing $x_0$. On arrival of a symbol, prepend it to the current address and trim to the stopping threshold. The selected interval contains the true value by induction. Its midpoint is within the allowed error. A bounded ordinal numbers the leaves left to right. The codec unpacks the ordinal, updates the address and ranks the new address using Fibonacci counts. Its decoder returns an integer numerator and denominator. [Executable update and decoder](FibonacciMemory.lean#L696)

At $n=2$ the addresses in order are 0, 10 and 11. Their intervals are $[0,1/4]$, $[1/2,5/8]$ and $[3/4,1]$. With startup zero, ordinal 0 reads $1/8$, which meets the tolerance at time zero. Reading 1 changes the ordinal to 1 and the readout to $9/16$. Reading 0 next returns to ordinal 0 and reads the true value $1/8$ exactly.

## Correspondence with the book

| Result | Lean entry point | Conclusion |
|---|---|---|
| Lemma 3.13.1 | [minimum_cover](FibonacciMemory.lean#L494) | The least number of arbitrarily translated closed intervals of length $2^{-n}$ covering a nonempty invariant set inside $[0,1]$ is $F_{n+2}$. Invariance here means equality with the union of its two images. |
| Theorem 3.13.2 | [ninety_bits](FibonacciMemory.lean#L424) | The threshold at startup zero, with the exact concrete value at $n=128$. |
| Proposition 3.13.3 | [ordinal_accuracy](FibonacciMemory.lean#L747) | The integer update and rational decoder meet the real-valued contract at every finite history. The book gives a separate bit-complexity analysis of its algorithm; kernel verification here covers the codec's input-output behavior. |
| Theorems 3.3.1 and 3.13.4 | [minimum_states](FibonacciMemory.lean#L380) · [minimum_bits](FibonacciMemory.lean#L441) | Exact cardinality and fixed-width capacity for each admissible startup; the [packing bound](FibonacciMemory.lean#L368) also applies outside the leaf union. |
| Proposition 3.13.6 | [nondyadic_bounds](FibonacciMemory.lean#L608) | With $\theta=2^{-n}$, the same exact optimum holds for $\theta/2\le\epsilon<3\theta/4$. For $3\theta/4\le\epsilon<\theta$, the proved bounds are $F_{n+1}$ and $F_{n+2}$. |

## Reproduce or inspect a trace

Use Lean 4.19.0 with Lake, Git, curl and Python 3. From this directory run:

```bash
python ../checks/formal.py --real
python ../checks/fibonacci_codec.py
echo '2 0 10' | lake env lean --run Trace.lean
```

The first command restores the dependencies in the Lake manifest and fetches compiled objects for the imported mathlib modules on demand. It builds the project and checks the transitive axiom report for every declared theorem. The second command compares the runnable Lean code to an independent Python implementation using exact fractions. For thresholds 1 through 12, it exercises both edges of every state. For thresholds 32, 64 and 128, it checks constant, alternating and seeded random traces. The report counts the startup readout once per trace.

The trace driver takes a threshold, initial ordinal and chronological sequence on each input line. Use a hyphen for an empty sequence. It prints the initial state and each subsequent state as three integers: the ordinal, midpoint numerator and positive denominator. The example produces:

```text
0 1 8
1 9 16
0 1 8
```

[Machine-readable run and theorem-by-theorem axiom report](../checks/fibonacci-evidence.json) · [All formal declarations](Audit.lean) · [Full monograph proof index](../INDEX.md)
