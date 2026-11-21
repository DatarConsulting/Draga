import Std

/-!
# Finite arguments in the DRAGA monograph

Words are read from left to right. All continuation statements quantify over
arbitrary finite words, including the empty word. The file has no external
dependencies beyond Lean's standard library. The analytic assumptions of the
book enter only where they are explicitly supplied as hypotheses.
-/
namespace Draga

section Words
variable {A S Y : Type}

/-- The Hankel entry at a prefix and suffix is the coefficient of their concatenation. -/
def entry (f : List A → Y) (u v : List A) : Y := f (u ++ v)

/-- Moving the boundary across one letter preserves the entry. This gives the
coordinate identity behind the shift intertwining in Chapter 1. -/
theorem shiftEntry (f : List A → Y) (u v : List A) (a : A) :
    entry f (u ++ [a]) v = entry f u (a :: v) := by simp [entry, List.append_assoc]

/-- Two products differing only in the order of two letters agree for scalar
multiplicative coefficients, the algebraic obstruction used in §1.12. -/
theorem scalarSwap {R : Type} (op : R → R → R) (comm : ∀ x y, op x y = op y x)
    (f : List A → R) (hom : ∀ u v, f (u ++ v) = op (f u) (f v)) (a b : A) :
    f [a, b] = f [b, a] := by
  calc f [a, b] = op (f [a]) (f [b]) := hom [a] [b]
       _ = op (f [b]) (f [a]) := comm _ _
       _ = f [b, a] := (hom [b] [a]).symm

/-- Executing a word from a state uses only the transition and current letter. -/
def run (step : S → A → S) (s : S) (w : List A) : S := w.foldl step s

/-- Restarting a deterministic machine at a checkpoint preserves every subsequent state. -/
theorem replay (step : S → A → S) (s : S) (u v : List A) :
    run step s (u ++ v) = run step (run step s u) v := by simp [run, List.foldl_append]

/-- Exact equality of observable behaviour for all common continuations. -/
def equivalent (read : List A → Y) (u v : List A) : Prop :=
  ∀ w, read (u ++ w) = read (v ++ w)

/-- The future-behaviour equivalence in Proposition 3.10.2 is a right congruence. -/
theorem rightCongruence (read : List A → Y) {u v : List A}
    (h : equivalent read u v) (a : A) : equivalent read (u ++ [a]) (v ++ [a]) := by
  intro w
  simpa [List.append_assoc] using h (a :: w)

/-- A deterministic machine cannot merge histories with different future outputs.
The output type can itself be a function of the query. -/
theorem distinguish (step : S → A → S) (out : S → Y) (s : S) (u v w : List A)
    (diff : out (run step s (u ++ w)) ≠ out (run step s (v ++ w))) :
    run step s u ≠ run step s v := by
  intro equal
  apply diff
  simp only [replay, equal]

/-- The inductive core of Proposition 3.10.3. Cell membership is a relation,
so cells may overlap and the argument requires no inverse transition. -/
theorem invariantCover {I : Type} (step : S → A → S) (update : I → A → I)
    (cell : I → S → Prop)
    (closed : ∀ i x a, cell i x → cell (update i a) (step x a))
    (w : List A) (i : I) (x : S) (initial : cell i x) :
    cell (run update i w) (run step x w) := by
  induction w generalizing i x with
  | nil => exact initial
  | cons a w ih => exact ih (update i a) (step x a) (closed i x a initial)

/-- Apply a query-uniform readout certificate along every history. In the book
`good` asserts that the two answer laws are within the specified tolerance. -/
theorem certifiedReadout {I Q : Type} (step : S → A → S) (update : I → A → I)
    (cell : I → S → Prop) (good : I → S → Q → Prop)
    (closed : ∀ i x a, cell i x → cell (update i a) (step x a))
    (cert : ∀ i x, cell i x → ∀ q, good i x q)
    (w : List A) (i : I) (x : S) (initial : cell i x) (q : Q) :
    good (run update i w) (run step x w) q :=
  cert _ _ (invariantCover step update cell closed w i x initial) q
end Words


section Quotient
variable {A Y S I : Type}

/-- Observational equivalence is a setoid even with infinitely many histories
and arbitrary output types. A value of Y can encode the answers to all queries. -/
def futureSetoid (read : List A → Y) : Setoid (List A) where
  r := equivalent read
  iseqv := {
    refl := fun _ _ => rfl
    symm := fun h _ => (h _).symm
    trans := fun h k _ => (h _).trans (k _) }

/-- One transition of the canonical quotient machine of Proposition 3.10.2. -/
def quotientStep (read : List A → Y) (q : Quotient (futureSetoid read)) (a : A) :
    Quotient (futureSetoid read) :=
  Quotient.liftOn q (fun u => Quotient.mk (futureSetoid read) (u ++ [a]))
    (fun _ _ h => Quotient.sound (rightCongruence read h a))

/-- Empty continuation supplies the output on a class. -/
def quotientRead (read : List A → Y) (q : Quotient (futureSetoid read)) : Y :=
  Quotient.liftOn q read (fun _ _ h => by simpa using h [])

/-- The quotient transition realizes concatenation exactly. In particular,
starting at the class of the empty word realizes the given answer function. -/
theorem quotientRuns (read : List A → Y) (u w : List A) :
    run (quotientStep read) (Quotient.mk (futureSetoid read) u) w =
      Quotient.mk (futureSetoid read) (u ++ w) := by
  induction w generalizing u with
  | nil => simp [run]
  | cons a w ih =>
    change run (quotientStep read) (Quotient.mk (futureSetoid read) (u ++ [a])) w = _
    simpa [List.append_assoc] using ih (u ++ [a])

/-- Exact quotient readout for any history and any future. -/
theorem quotientAnswers (read : List A → Y) (u w : List A) :
    quotientRead read (run (quotientStep read) (Quotient.mk (futureSetoid read) u) w) =
      read (u ++ w) := by rw [quotientRuns]; rfl

/-- Any state reached in an exact deterministic machine determines a unique
future-equivalence class. Thus exact machines cannot identify distinct classes. -/
theorem minimalClasses (read : List A → Y) (step : S → A → S) (out : S → Y) (start : S)
    (exactOutput : ∀ w, out (run step start w) = read w) (u v : List A)
    (same : run step start u = run step start v) :
    Quotient.mk (futureSetoid read) u = Quotient.mk (futureSetoid read) v := by
  apply Quotient.sound
  intro w
  rw [← exactOutput, ← exactOutput, replay, replay, same]

/-- Packing as an injection. Two histories are separated if some future
admits no common approximate answer. A metric gap greater than twice the
uniform error supplies this hypothesis by the triangle inequality. -/
theorem packingInjection (step : S → A → S) (out : S → Y) (start : S)
    (read : List A → Y) (close : Y → Y → Prop) (hist : I → List A)
    (accurate : ∀ w, close (read w) (out (run step start w)))
    (separate : ∀ i j, i ≠ j → ∃ w, ∀ y,
      ¬ (close (read (hist i ++ w)) y ∧ close (read (hist j ++ w)) y)) :
    ∀ i j, run step start (hist i) = run step start (hist j) → i = j := by
  intro i j same
  apply Classical.byContradiction
  intro ne
  obtain ⟨w, witness⟩ := separate i j ne
  have h := accurate (hist i ++ w)
  have k := accurate (hist j ++ w)
  rw [replay, same] at h
  rw [replay] at k
  exact witness _ ⟨h,k⟩
end Quotient

section Maps
variable {K V X Y : Type}

/-- Absence is an explicit value. Writing none deletes the key. -/
inductive Edit (K V : Type)
  | write (key : K) (value : Option V)
  | skip

/-- Deterministic set, delete and no-op semantics from Proposition 3.15.9. -/
def editMap [DecidableEq K] (m : K → Option V) : Edit K V → K → Option V
  | .skip => m
  | .write k v => fun j => if j = k then v else m j

/-- Functions on active maps separate whenever a key query can separate them. -/
theorem mapsSeparated (m n : K → Option V) (different : m ≠ n) : ∃ k, m k ≠ n k := by
  apply Classical.byContradiction
  intro h
  apply different
  funext k
  apply Classical.byContradiction
  intro hk
  exact h ⟨k,hk⟩

/-- The factorization in Proposition 3.15.9 holds for arbitrary auxiliary
continuous state X: an update may consult the old map and the next edit.
No finiteness, invertibility, or contraction is required for this implication. -/
theorem activeMapFuture [DecidableEq K] (base : K → Option V)
    (update : ((K → Option V) × X) → Edit K V → ((K → Option V) × X))
    (answer : ((K → Option V) × X) → Y) (u v w : List (Edit K V)) (x : X)
    (same : run editMap base u = run editMap base v) :
    answer (run update (run editMap base u, x) w) =
      answer (run update (run editMap base v, x) w) := by rw [same]

/-- A one-step commuting diagram lifts to an arbitrary sequence of events.
Taking the concrete state to be a journal and the projection to be its active
map (together with any continuous state) proves exact forward simulation.
This is the substantive factorization assumption of Proposition 3.15.9. -/
theorem projectRun {Raw Small E : Type} (project : Raw → Small)
    (rawStep : Raw → E → Raw) (smallStep : Small → E → Small)
    (commute : ∀ x a, project (rawStep x a) = smallStep (project x) a)
    (x : Raw) (w : List E) :
    project (run rawStep x w) = run smallStep (project x) w := by
  induction w generalizing x with
  | nil => rfl
  | cons a w ih =>
    change project (run rawStep (rawStep x a) w) = run smallStep (smallStep (project x) a) w
    rw [ih, commute]

end Maps

section Budget
variable {I B : Type}

/-- Largest numerator at one block. The empty state list contributes zero. -/
def blockMax (weight : I → B → Nat) (b : B) : List I → Nat
  | [] => 0
  | i :: rest => max (weight i b) (blockMax weight b rest)

/-- Sum over the stated finite list, retaining multiplicities. -/
def total {T : Type} (f : T → Nat) : List T → Nat
  | [] => 0
  | t :: ts => f t + total f ts

/-- Summation preserves addition. -/
theorem totalAdd {T : Type} (f g : T → Nat) (ts : List T) :
    total (fun t => f t + g t) ts = total f ts + total g ts := by
  induction ts with
  | nil => rfl
  | cons t ts h => simp only [total, h]; omega

/-- A pointwise order on the listed elements survives summation. -/
theorem totalMono {T : Type} (f g : T → Nat) (ts : List T)
    (h : ∀ t ∈ ts, f t ≤ g t) : total f ts ≤ total g ts := by
  induction ts with
  | nil => exact Nat.le_refl 0
  | cons t ts ih =>
    exact Nat.add_le_add (h t (by simp)) (ih (fun j hj => h j (by simp [hj])))

/-- Core finite inequality of Theorem 2.7.6. For rational substochastic rows,
multiply all entries by a common denominator D. The sum of columnwise maxima
then costs at most one row budget per state. No distinctness or positivity
assumption is needed beyond the natural-number representation. -/
theorem blockBudget (weight : I → B → Nat) (blocks : List B) (states : List I)
    (D : Nat) (budget : ∀ i ∈ states, total (weight i) blocks ≤ D) :
    total (fun b => blockMax weight b states) blocks ≤ states.length * D := by
  induction states with
  | nil =>
    have hz : total (fun (_ : B) => 0) blocks = 0 := by
      induction blocks with
      | nil => rfl
      | cons b bs ih => simpa [total] using ih
    simpa [blockMax] using Nat.le_of_eq hz
  | cons i is ih =>
    have prev := ih (fun j hj => budget j (by simp [hj]))
    have row := budget i (by simp)
    have h := totalMono (fun b => blockMax weight b (i :: is))
      (fun b => weight i b + blockMax weight b is) blocks (by
        intro b _
        simp only [blockMax]
        exact Nat.max_le.mpr ⟨by omega, by omega⟩)
    rw [totalAdd] at h
    simp only [List.length_cons, Nat.succ_mul]
    omega
end Budget
section Stopping

/-- Binary stopping tree. At residual budget zero the current word is a leaf. -/
inductive StoppingTree
  | leaf
  | split (one two : StoppingTree)

/-- A letter of cost one or two is appended until the cumulative cost reaches
or exceeds the threshold. At threshold one either first letter stops. -/
def stopTree : Nat → StoppingTree
  | 0 => .leaf
  | 1 => .split .leaf .leaf
  | n + 2 => .split (stopTree (n+1)) (stopTree n)

/-- Count terminal words, without including internal prefixes. -/
def leaves : StoppingTree → Nat
  | .leaf => 1
  | .split a b => leaves a + leaves b

/-- Conventional normalization, with Fibonacci numbers zero and one at the origin. -/
def fib : Nat → Nat
  | 0 => 0
  | 1 => 1
  | n + 2 => fib (n+1) + fib n

/-- The exact number of stopped leaves in Lemma 3.13.1 for every threshold.
The geometric separation of endpoints is a separate part of that lemma. -/
theorem fibonacciLeaves (n : Nat) : leaves (stopTree n) = fib (n+2) := by
  suffices h : leaves (stopTree n) = fib (n+2) ∧
               leaves (stopTree (n+1)) = fib (n+3) from h.1
  induction n with
  | zero => decide
  | succ n ih =>
    refine ⟨ih.2, ?_⟩
    change leaves (stopTree (n+1)) + leaves (stopTree n) = fib (n+3) + fib (n+2)
    rw [ih.1, ih.2]

end Stopping

section Radix

/-- A fixed ordered domain of n keys with b possible values per key. Digit zero
can represent absence. The head is the least significant digit. -/
inductive Digits (b : Nat) : Nat → Type
  | nil : Digits b 0
  | cons {n} (digit : Fin b) (tail : Digits b n) : Digits b (n+1)

/-- Rank a persistent map in little-endian radix order. -/
def rankDigits {b n : Nat} : Digits b n → Nat
  | .nil => 0
  | .cons d t => d.val + b * rankDigits t

/-- Read exactly n base-b digits; the caller may prove a range bound to obtain
the inverse in both directions. -/
def unrankDigits (b : Nat) (pos : 0 < b) : (n : Nat) → Nat → Digits b n
  | 0, _ => .nil
  | n+1, k => .cons ⟨k % b, Nat.mod_lt k pos⟩ (unrankDigits b pos n (k / b))

/-- Every map identifier is in the claimed finite range. -/
theorem rankBound {b n : Nat} (m : Digits b n) : rankDigits m < b^n := by
  induction m with
  | nil => simp [rankDigits]
  | @cons n d t ih =>
    have scale := Nat.mul_le_mul_left b ih
    rw [Nat.mul_succ] at scale
    simp only [rankDigits, Nat.pow_succ]
    rw [Nat.mul_comm (b^n) b]
    have small := d.isLt
    omega

/-- The encoder is lossless for every map, including maps with high zero digits. -/
theorem decodeEncode {b n : Nat} (pos : 0 < b) (m : Digits b n) :
    unrankDigits b pos n (rankDigits m) = m := by
  induction m with
  | nil => rfl
  | cons d t ih =>
    simp [unrankDigits, rankDigits, Nat.add_mul_mod_self_left, Nat.add_mul_div_left _ _ pos,
          Nat.mod_eq_of_lt d.isLt, Nat.div_eq_of_lt d.isLt, ih]

/-- Every integer below the stated cardinality is the code of a map. -/
theorem encodeDecode (b : Nat) (pos : 0 < b) (n k : Nat) (within : k < b^n) :
    rankDigits (unrankDigits b pos n k) = k := by
  induction n generalizing k with
  | zero => simp at within; simp [unrankDigits, rankDigits, within]
  | succ n ih =>
    have lower : k / b < b^n := (Nat.div_lt_iff_lt_mul pos).mpr (by simpa [Nat.pow_succ] using within)
    simp only [unrankDigits, rankDigits, ih _ lower, Nat.mod_add_div]

/-- Finite pigeonhole principle on initial intervals. In the induction step,
remove the last image and close the gap above it. -/
theorem slotsLower (n m : Nat) (f : Nat → Nat)
    (bound : ∀ i, i < n → f i < m)
    (distinct : ∀ i j, i < n → j < n → f i = f j → i = j) : n ≤ m := by
  induction n generalizing m f with
  | zero => omega
  | succ n ih =>
    let z := f n
    let squeeze := fun y => if y < z then y else y-1
    have zsmall : z < m := bound n (by omega)
    have away (i : Nat) (hi : i < n) : f i ≠ z := by
      intro e
      have eq := distinct i n (by omega) (by omega) e
      omega
    have small (i : Nat) (hi : i < n) : squeeze (f i) < m-1 := by
      have := away i hi
      have := bound i (by omega)
      simp only [squeeze]
      split <;> omega
    have inject (i j : Nat) (hi : i < n) (hj : j < n) (eq : squeeze (f i) = squeeze (f j)) : i = j := by
      have := away i hi
      have := away j hj
      apply distinct i j (by omega) (by omega)
      by_cases a : f i < z <;> by_cases b : f j < z <;> simp [squeeze, a, b] at eq <;> omega
    have reduced := ih (m-1) (fun i => squeeze (f i)) small inject
    omega

/-- Any injective bounded identifier of maps needs at least b^n slots.
For bit strings of width w take capacity = 2^w. Together with rankBound,
this is the fixed-width optimality claim of Proposition 3.15.9 under exact
queries and full reachability. -/
theorem mapCapacity (b n capacity : Nat) (pos : 0 < b) (enc : Digits b n → Nat)
    (bounded : ∀ m, enc m < capacity)
    (faithful : ∀ m t, enc m = enc t → m = t) : b^n ≤ capacity := by
  apply slotsLower (b^n) capacity (fun i => enc (unrankDigits b pos n i))
  · intro i _; exact bounded _
  · intro i j hi hj same
    have e := congrArg rankDigits (faithful _ _ same)
    simpa only [encodeDecode b pos n i hi, encodeDecode b pos n j hj] using e

end Radix


/-- Pascal coefficients indexed on all integers. Out-of-range values are zero
by the recurrence, so convolution has no endpoint conventions to discharge. -/
def binomial : Nat → Int → Nat
  | 0, r => if r = 0 then 1 else 0
  | n+1, r => binomial n r + binomial n (r-1)

/-- Sum a reward over all binary words of length n, evaluated at their number
of ones. The two recursive branches add a zero or a one. -/
def pathSum : Nat → (Int → Nat) → Nat
  | 0, reward => reward 0
  | n+1, reward => pathSum n reward + pathSum n (fun j => reward (j+1))

/-- Splitting a binary word after a positions gives Vandermonde convolution. -/
theorem vandermonde (a b : Nat) (r : Int) :
    pathSum a (fun j => binomial b (r-j)) = binomial (a+b) r := by
  induction a generalizing r with
  | zero => simp [pathSum]
  | succ a ih =>
    have shift : (fun j : Int => binomial b (r-(j+1))) = (fun j => binomial b ((r-1)-j)) := by
      funext j; congr 1; omega
    simp only [pathSum, shift, ih, Nat.succ_add, binomial]

/-- Count how many times each index occurs in the path sum. -/
theorem coefficientMultiplicity (n : Nat) (r : Int) :
    pathSum n (fun j => if j = r then 1 else 0) = binomial n r := by
  have same : (fun j : Int => if j = r then 1 else 0) = (fun j => binomial 0 (r-j)) := by
    funext j
    by_cases h : j = r
    · simp [h, binomial]
    · have nz : r-j ≠ 0 := by omega
      simp [h, binomial, nz]
  rw [same]
  simpa using vandermonde n 0 r

/-- Interchanging zeros and ones reflects the Pascal row, at all integer indices. -/
theorem reflectBinomial (n : Nat) (r : Int) : binomial n ((n:Int)-r) = binomial n r := by
  induction n generalizing r with
  | zero =>
    by_cases h : r = 0
    · simp [h,binomial]
    · have hn : -r ≠ 0 := by omega
      simp [binomial,h,hn]
  | succ n ih =>
    have t : ((n+1:Nat):Int)-r = (n:Int)-(r-1) := by omega
    have u : ((n+1:Nat):Int)-r-1 = (n:Int)-r := by omega
    change binomial n (((n+1:Nat):Int)-r) + binomial n (((n+1:Nat):Int)-r-1) = _
    rw [u,t,ih,ih]
    exact Nat.add_comm _ _

/-- The path-sum version of the central-binomial identity in Theorem 1.8.1.
By coefficientMultiplicity each j contributes binomial(k,j) copies of the
same binomial(k,j), giving the familiar squared-coefficient sum. -/
theorem centralBinomial (k : Nat) :
    pathSum k (binomial k) = binomial (k+k) (k:Int) := by
  have row : binomial k = (fun j => binomial k ((k:Int)-j)) := by
    funext j; exact (reflectBinomial k j).symm
  calc pathSum k (binomial k) = pathSum k (fun j => binomial k ((k:Int)-j)) := congrArg (pathSum k) row
       _ = binomial (k+k) (k:Int) := vandermonde k k k
/-- The Rayleigh quotient in Theorem 1.8.1 uses t=q², 0≤t≤1. Write t=u/v, S=s/e and K=k/e with positive denominators. Clearing both
denominators gives the order certificate below. Its slack is (v-u)e(k-s). -/
theorem rayleighCross (u v e s k : Nat) (hu : u ≤ v) (hs : s ≤ k) :
    (v*e + u*k)*(e+s) ≤ (v*e + u*s)*(e+k) := by
  obtain ⟨a,rfl⟩ := Nat.exists_eq_add_of_le hu
  obtain ⟨b,rfl⟩ := Nat.exists_eq_add_of_le hs
  have identity : ((u+a)*e+u*s)*(e+(s+b)) =
      ((u+a)*e+u*(s+b))*(e+s) + a*e*b := by
    simp only [Nat.add_mul,Nat.mul_add]
    ac_rfl
  omega

end Draga
