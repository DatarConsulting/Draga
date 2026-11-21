import Draga
import Mathlib.Data.Real.Archimedean
import Mathlib.Data.List.Pairwise
import Mathlib.Data.Fintype.Card
import Mathlib.Data.Nat.Log
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.Ring
import Mathlib.Tactic.NormNum

/-!
# Exact memory for the two dyadic contractions

A history is chronological. An address is newest letter first. The maps on
the real unit interval are x ↦ x/4 and x ↦ (1+x)/2. All time horizons below
are finite, with no upper bound. A memory state is numbered by a value of
Fin m and its decoder can be any function into the real numbers.
-/
namespace FibonacciMemory

abbrev Word := List Bool

/-- A false letter costs two dyadic scales and a true letter costs one. -/
def cost : Word → Nat
  | [] => 0
  | false :: w => 2 + cost w
  | true :: w => 1 + cost w

/-- The leaves in left-to-right order, stopping as soon as the threshold is met. -/
def words : Nat → List Word
  | 0 => [[]]
  | 1 => [[false],[true]]
  | n+2 => (words n).map (false :: ·) ++ (words (n+1)).map (true :: ·)

/-- Binary threshold, expressed as a length in the real unit interval. -/
noncomputable def theta (n : Nat) : ℝ := (1/2 : ℝ)^n

/-- One physical update. -/
noncomputable def step : Bool → ℝ → ℝ
  | false, x => x/4
  | true, x => (1+x)/2

/-- Evaluate a newest-first address on its remaining tail value. -/
noncomputable def value (w : Word) (x : ℝ) := w.foldr step x
noncomputable def leftEnd (w : Word) := value w 0
noncomputable def rightEnd (w : Word) := value w 1
noncomputable def width (w : Word) := rightEnd w - leftEnd w

/-- The complete startup domain covered by the finite leaf machine. -/
def leafUnion (n : Nat) (x : ℝ) : Prop :=
  ∃ w ∈ words n, leftEnd w ≤ x ∧ x ≤ rightEnd w

theorem theta_pos (n : Nat) : 0 < theta n := pow_pos (by norm_num) n

theorem theta_step (n : Nat) : theta (n+1) = theta n / 2 := by
  simp [theta, pow_succ, div_eq_mul_inv]

theorem theta_le_one (n : Nat) : theta n ≤ 1 := by
  induction n with
  | zero => norm_num [theta]
  | succ n ih => rw [theta_step]; linarith [theta_pos n]

theorem value_cons (a : Bool) (w : Word) (x : ℝ) : value (a::w) x = step a (value w x) := rfl

theorem value_append (v w : Word) (x : ℝ) : value (v++w) x = value v (value w x) := by
  simp [value, List.foldr_append]

theorem step_mono (a : Bool) : Monotone (step a) := by
  intro x y h; cases a <;> simp only [step] <;> linarith

theorem value_mono (w : Word) : Monotone (value w) := by
  induction w with
  | nil => exact fun _ _ h => h
  | cons a w ih => exact (step_mono a).comp ih

theorem value_unit (w : Word) {x : ℝ} (hx : 0 ≤ x ∧ x ≤ 1) :
    0 ≤ value w x ∧ value w x ≤ 1 := by
  induction w with
  | nil => exact hx
  | cons a w ih =>
    rw [value_cons]; cases a <;> simp only [step] <;> constructor <;> linarith [ih.1,ih.2]

/-- Affine interpolation also covers the two boundary values of the tail. -/
theorem interpolation (w : Word) (x : ℝ) :
    value w x = leftEnd w + width w * x := by
  induction w with
  | nil => simp [value,leftEnd,rightEnd,width]
  | cons a w ih =>
    unfold width leftEnd rightEnd at *
    cases a <;> simp only [value_cons,step] <;> rw [ih] <;> ring

theorem word_count (n : Nat) : (words n).length = Draga.fib (n+2) := by
  induction n using Nat.twoStepInduction with
  | zero => rfl
  | one => rfl
  | more n h k =>
    simp only [words,List.length_append,List.length_map,h,k]
    exact Nat.add_comm _ _

/-- Every stopped word has enough cost and uses at most n actual letters. -/
theorem leaf_cost_length (n : Nat) : ∀ w ∈ words n, n ≤ cost w ∧ w.length ≤ n := by
  induction n using Nat.twoStepInduction with
  | zero => simp [words,cost]
  | one => simp [words,cost]
  | more n h k =>
    intro w hw
    simp only [words,List.mem_append,List.mem_map] at hw
    rcases hw with ⟨u,hu,rfl⟩ | ⟨u,hu,rfl⟩
    · obtain ⟨hc,hl⟩ := h u hu
      simp only [cost,List.length_cons]; omega
    · obtain ⟨hc,hl⟩ := k u hu
      simp only [cost,List.length_cons]; omega

theorem leaf_width (n : Nat) : ∀ w ∈ words n, width w = theta n ∨ width w = theta n/2 := by
  induction n using Nat.twoStepInduction with
  | zero => simp [words,width,leftEnd,rightEnd,value,theta]
  | one => norm_num [words,width,leftEnd,rightEnd,value,step,theta]
  | more n h k =>
    intro w hw
    simp only [words,List.mem_append,List.mem_map] at hw
    rcases hw with ⟨u,hu,rfl⟩ | ⟨u,hu,rfl⟩
    · have d : width (false::u) = width u / 4 := by
        simp only [width,leftEnd,rightEnd,value_cons,step]; ring
      rw [d,theta_step,theta_step]
      rcases h u hu with e|e
      · left; linarith
      · right; linarith
    · have d : width (true::u) = width u / 2 := by
        simp only [width,leftEnd,rightEnd,value_cons,step]; ring
      rw [d,theta_step]
      rcases k u hu with e|e
      · left; linarith
      · right; linarith

/-- Any two ordered leaf endpoints, including adjacent ones, are at least
three half-thresholds apart. The cross-branch estimate uses the right
leaf's positive width in addition to the gap between the two children. -/
theorem endpoint_spacing (n : Nat) :
    (words n).Pairwise (fun u v => rightEnd u + 3*theta n/2 ≤ rightEnd v) := by
  induction n using Nat.twoStepInduction with
  | zero => simp [words]
  | one => norm_num [words,rightEnd,value,step,theta]
  | more n h k =>
    rw [words,List.pairwise_append]
    refine ⟨List.pairwise_map.mpr (h.imp ?_),List.pairwise_map.mpr (k.imp ?_),?_⟩
    · intro u v d
      change rightEnd u / 4 + 3*theta (n+2)/2 ≤ rightEnd v / 4
      rw [theta_step,theta_step]; linarith
    · intro u v d
      change (1+rightEnd u)/2 + 3*theta (n+2)/2 ≤ (1+rightEnd v)/2
      rw [theta_step]; linarith
    · intro a ha b hb
      obtain ⟨u,_,rfl⟩ := List.mem_map.mp ha
      obtain ⟨v,hv,rfl⟩ := List.mem_map.mp hb
      change rightEnd u / 4 + 3*theta (n+2)/2 ≤ (1+rightEnd v)/2
      have hu : rightEnd u ≤ 1 := (value_unit u (x:=1) (by norm_num)).2
      have hl : 0 ≤ leftEnd v := (value_unit v (x:=0) (by norm_num)).1
      have hw : theta (n+1)/2 ≤ width v := by
        rcases leaf_width _ _ hv with q|q <;> linarith [theta_pos (n+1)]
      unfold width at hw
      rw [theta_step] at hw
      rw [theta_step,theta_step]
      linarith [theta_le_one n]

/-- Keep the shortest newest-first prefix whose cost reaches n. An insufficient
word is returned in full. Once a leaf is stored, that case cannot occur. -/
def trim : Nat → Word → Word
  | 0, _ => []
  | _+1, [] => []
  | n+1, true::w => true :: trim n w
  | n+1, false::w => false :: trim (n-1) w
termination_by _ w => w.length

theorem trim_prefix (w : Word) (n : Nat) : ∃ tail, w = trim n w ++ tail := by
  induction w generalizing n with
  | nil => cases n <;> exact ⟨[], by simp [trim]⟩
  | cons a w ih =>
    cases n with
    | zero => exact ⟨a::w,by simp [trim]⟩
    | succ n =>
      cases a with
      | false =>
        obtain ⟨t,ht⟩ := ih (n-1)
        exact ⟨t, by simpa [trim] using congrArg (false::·) ht⟩
      | true =>
        obtain ⟨t,ht⟩ := ih n
        exact ⟨t, by simpa [trim] using congrArg (true::·) ht⟩

theorem trim_leaf (n : Nat) : ∀ w, n ≤ cost w → trim n w ∈ words n := by
  induction n using Nat.twoStepInduction with
  | zero => simp [words,trim]
  | one =>
    intro w h
    cases w with
    | nil => simp [cost] at h
    | cons a w => cases a <;> simp [words,trim]
  | more n h k =>
    intro w hw
    cases w with
    | nil => simp [cost] at hw
    | cons a w =>
      cases a with
      | false =>
        simp only [cost] at hw
        have main := h w (by omega)
        exact List.mem_append_left _ (List.mem_map.mpr ⟨_,main,by simp [trim]⟩)
      | true =>
        simp only [cost] at hw
        have main := k w (by omega)
        exact List.mem_append_right _ (List.mem_map.mpr ⟨_,main,by simp [trim]⟩)

theorem trim_update_leaf {n : Nat} {w : Word} (hw : w ∈ words n) (a : Bool) :
    trim n (a::w) ∈ words n := by
  apply trim_leaf
  have h := (leaf_cost_length n w hw).1
  cases a <;> simp only [cost] <;> omega

/-- One forward update preserves interval containment, without inverting a map. -/
theorem trim_invariant (n : Nat) (a : Bool) (w : Word) (x : ℝ)
    (hx : leftEnd w ≤ x ∧ x ≤ rightEnd w) :
    leftEnd (trim n (a::w)) ≤ step a x ∧ step a x ≤ rightEnd (trim n (a::w)) := by
  let v := trim n (a::w)
  obtain ⟨tail,ht⟩ := trim_prefix (a::w) n
  change leftEnd v ≤ step a x ∧ step a x ≤ rightEnd v
  have inside : leftEnd (a::w) ≤ step a x ∧ step a x ≤ rightEnd (a::w) :=
    ⟨step_mono a hx.1,step_mono a hx.2⟩
  change a::w = v++tail at ht
  rw [ht] at inside
  unfold leftEnd rightEnd at *
  simp only [value_append] at inside
  exact ⟨le_trans (value_mono v (value_unit tail (x:=0) (by norm_num)).1) inside.1,
         le_trans inside.2 (value_mono v (value_unit tail (x:=1) (by norm_num)).2)⟩

/-- An autonomous machine using m labelled states. Accuracy includes the empty
history. Decoder values are unrestricted real numbers. -/
def FitsAt (ε x₀ : ℝ) (m : Nat) : Prop :=
  ∃ (start : Fin m) (update : Fin m → Bool → Fin m) (read : Fin m → ℝ),
    ∀ h : Word, |read (h.foldl update start) - value h.reverse x₀| ≤ ε

/-- Specialization to half a dyadic leaf width. -/
def Fits (n : Nat) (x₀ : ℝ) (m : Nat) := FitsAt (theta n/2) x₀ m

local instance (priority := 2000) : BEq Word := instBEqOfDecidableEq

abbrev Index (n : Nat) := Fin (words n).length
def address {n : Nat} (i : Index n) : Word := (words n).get i

/-- The index is executable by searching the canonical, ordered leaf list. -/
def indexOf {n : Nat} (w : Word) (hw : w ∈ words n) : Index n :=
  ⟨(words n).idxOf w,List.idxOf_lt_length hw⟩

theorem index_spec {n : Nat} (w : Word) (hw : w ∈ words n) :
    address (indexOf w hw) = w := by
  exact List.getElem_idxOf (List.idxOf_lt_length hw)

def updateLeaf (n : Nat) (i : Index n) (a : Bool) : Index n :=
  indexOf (trim n (a::address i)) (trim_update_leaf (List.get_mem _ _) a)

theorem update_address (n : Nat) (i : Index n) (a : Bool) :
    address (updateLeaf n i a) = trim n (a::address i) := index_spec _ _

/-- Any implementation whose address update commutes with trimming preserves
the leaf interval for every chronological history. This argument is shared
by the list-index realization and the executable table-free realization. -/
theorem history_interval {S : Type} (n : Nat) (addr : S → Word)
    (update : S → Bool → S)
    (commute : ∀ i a, addr (update i a) = trim n (a :: addr i))
    (h : Word) (i : S) (x : ℝ) (hx : leftEnd (addr i) ≤ x ∧ x ≤ rightEnd (addr i)) :
    leftEnd (addr (h.foldl update i)) ≤ value h.reverse x ∧
      value h.reverse x ≤ rightEnd (addr (h.foldl update i)) := by
  induction h generalizing i x with
  | nil => exact hx
  | cons a h ih =>
    have hz := trim_invariant n a (addr i) x hx
    rw [← commute i a] at hz
    simpa [List.reverse_cons,value_append,value] using ih (update i a) (step a x) hz

/-- The midpoint decoder meets the prescribed tolerance for any history from
a permitted starting interval. No property of the unknown trajectory is used
to choose future memory updates. -/
theorem machine_upper (n : Nat) (x₀ : ℝ) (startup : leafUnion n x₀) :
    Fits n x₀ (Draga.fib (n+2)) := by
  rw [← word_count]
  obtain ⟨w,hw,hlo,hhi⟩ := startup
  let start := indexOf w hw
  let read := fun i : Index n => (leftEnd (address i)+rightEnd (address i))/2
  refine ⟨start,updateLeaf n,read,?_⟩
  intro h
  have ini : leftEnd (address start) ≤ x₀ ∧ x₀ ≤ rightEnd (address start) := by
    simpa [start,index_spec] using And.intro hlo hhi
  obtain ⟨l,r⟩ := history_interval n address (updateLeaf n) (update_address n) h start x₀ ini
  have bound := leaf_width n _ (List.get_mem (words n) (h.foldl (updateLeaf n) start))
  change width (address (h.foldl (updateLeaf n) start)) = theta n ∨
    width (address (h.foldl (updateLeaf n) start)) = theta n/2 at bound
  have d : width (address (h.foldl (updateLeaf n) start)) ≤ theta n := by
    rcases bound with t|t <;> linarith [theta_pos n]
  unfold width at d
  dsimp only [read]
  exact abs_le.mpr ⟨by linarith,by linarith⟩

/-- Three final warm-up ones ensure proximity to 1. Padding before them puts
every witness at the same time n+3, independent of the number of letters. -/
def witnessHistory (n : Nat) (w : Word) : Word :=
  (List.replicate (n-w.length) true ++ [true,true,true]) ++ w.reverse

theorem witness_time {n : Nat} {w : Word} (hw : w ∈ words n) :
    (witnessHistory n w).length = n+3 := by
  have d := (leaf_cost_length n w hw).2
  simp [witnessHistory]; omega

/-- This inequality concerns actual finite histories, including starts outside
the leaf union. The attractor endpoint need not itself be reached in finite time. -/
theorem witness_near {n : Nat} {w : Word} (hw : w ∈ words n) {x₀ : ℝ}
    (hx : 0 ≤ x₀ ∧ x₀ ≤ 1) :
    rightEnd w - theta n/8 ≤ value (witnessHistory n w).reverse x₀ ∧
      value (witnessHistory n w).reverse x₀ ≤ rightEnd w := by
  let z := value (List.replicate (n-w.length) true) x₀
  let y := step true (step true (step true z))
  have eq : value (witnessHistory n w).reverse x₀ = value w y := by
    simp [witnessHistory,List.reverse_append,value_append,y,z,value,List.foldr_append]
  have hz : 0 ≤ z ∧ z ≤ 1 := value_unit _ hx
  have hy : 7/8 ≤ y ∧ y ≤ 1 := by
    dsimp only [y,step]; constructor <;> linarith [hz.1,hz.2]
  have hbound := leaf_width n w hw
  have w0 : 0 ≤ width w := by rcases hbound with e|e <;> linarith [theta_pos n]
  have w1 : width w ≤ theta n := by rcases hbound with e|e <;> linarith [theta_pos n]
  rw [eq]
  refine ⟨?_, value_mono w hy.2⟩
  have s := mul_nonneg w0 (sub_nonneg.mpr hy.1)
  rw [interpolation]
  have endpoint : rightEnd w = leftEnd w + width w := by unfold width; ring
  rw [endpoint]; nlinarith

theorem witness_gap (n : Nat) {i j : Index n} (order : i < j) {x₀ : ℝ}
    (hx : 0 ≤ x₀ ∧ x₀ ≤ 1) :
    value (witnessHistory n (address i)).reverse x₀ + theta n <
      value (witnessHistory n (address j)).reverse x₀ := by
  have gap := List.pairwise_iff_get.mp (endpoint_spacing n) i j order
  change rightEnd (address i) + 3*theta n/2 ≤ rightEnd (address j) at gap
  have a := (witness_near (w:=address i) (List.get_mem _ i) hx).2
  have b := (witness_near (w:=address j) (List.get_mem _ j) hx).1
  linarith [theta_pos n]

/-- A family of actually reachable and strictly separated outputs injects into
the states of any uniformly accurate machine, for any real error tolerance. -/
theorem finite_packing {n m : Nat} {x₀ ε : ℝ} (hist : Index n → Word)
    (gap : ∀ i j : Index n, i < j → value (hist i).reverse x₀ + 2*ε < value (hist j).reverse x₀)
    (machine : FitsAt ε x₀ m) : Draga.fib (n+2) ≤ m := by
  obtain ⟨start,update,read,accurate⟩ := machine
  let encode := fun i : Index n => (hist i).foldl update start
  have distinct {i j : Index n} (lt : i < j) : encode i ≠ encode j := by
    intro same
    have a := (abs_le.mp (accurate (hist i))).2
    have b := (abs_le.mp (accurate (hist j))).1
    change read (encode i) - _ ≤ _ at a
    change _ ≤ read (encode j) - _ at b
    rw [same] at a
    linarith [gap i j lt]
  have inj : Function.Injective encode := by
    intro i j e
    rcases Nat.lt_trichotomy i.val j.val with lt|eq|gt
    · exact False.elim (distinct lt e)
    · exact Fin.ext eq
    · exact False.elim (distinct gt e.symm)
  have bound := Fintype.card_le_of_injective encode inj
  simpa [Index,word_count] using bound

/-- The dyadic lower bound holds at every starting point of the unit interval. -/
theorem machine_lower {n m : Nat} {x₀ : ℝ} (hx : 0 ≤ x₀ ∧ x₀ ≤ 1)
    (machine : Fits n x₀ m) : Draga.fib (n+2) ≤ m := by
  apply finite_packing (fun i : Index n => witnessHistory n (address i)) _ machine
  intro i j hij
  linarith [witness_gap n hij hx]

theorem union_unit {n : Nat} {x : ℝ} (hx : leafUnion n x) : 0 ≤ x ∧ x ≤ 1 := by
  obtain ⟨w,_,lo,hi⟩ := hx
  exact ⟨le_trans (value_unit w (x:=0) (by norm_num)).1 lo,
         le_trans hi (value_unit w (x:=1) (by norm_num)).2⟩

/-- Theorem 3.13.4, sharp state count at every permitted startup point. -/
theorem minimum_states {n : Nat} {x₀ : ℝ} (hx : leafUnion n x₀) :
    Fits n x₀ (Draga.fib (n+2)) ∧ ∀ m, Fits n x₀ m → Draga.fib (n+2) ≤ m :=
  ⟨machine_upper n x₀ hx,fun _ => machine_lower (union_unit hx)⟩

theorem zero_start (n : Nat) : leafUnion n 0 := by
  induction n using Nat.twoStepInduction with
  | zero => exact ⟨[],by simp [words],by norm_num [leftEnd,rightEnd,value]⟩
  | one => exact ⟨[false],by simp [words],by norm_num [leftEnd,rightEnd,value,step]⟩
  | more n ih _ =>
    obtain ⟨w,hw,lo,hi⟩ := ih
    refine ⟨false::w,List.mem_append_left _ (List.mem_map.mpr ⟨w,hw,rfl⟩),?_⟩
    change leftEnd w/4 ≤ 0 ∧ 0 ≤ rightEnd w/4
    constructor <;> linarith

/-- Theorem 3.13.2 at its prescribed startup point, including time-zero accuracy. -/
theorem zero_minimum (n : Nat) :
    Fits n 0 (Draga.fib (n+2)) ∧ ∀ m, Fits n 0 m → Draga.fib (n+2) ≤ m :=
  minimum_states (zero_start n)

/-- A fixed-width register can label at most 2^b states. This formulation
also allows unused bit patterns. -/
def BitFits (n : Nat) (x : ℝ) (b : Nat) : Prop :=
  ∃ m, Fits n x m ∧ m ≤ 2^b

theorem bit_threshold {n : Nat} {x : ℝ} (hx : leafUnion n x) (b : Nat) :
    BitFits n x b ↔ Draga.fib (n+2) ≤ 2^b := by
  constructor
  · rintro ⟨m,fit,cap⟩
    exact le_trans ((minimum_states hx).2 m fit) cap
  · intro cap
    exact ⟨_,(minimum_states hx).1,cap⟩

/-- Compute adjacent Fibonacci numbers with one recursive evaluation per index. -/
def fibPair : Nat → Nat × Nat
  | 0 => (0,1)
  | n+1 => let z := fibPair n; (z.2,z.1+z.2)

theorem fibPair_spec (n : Nat) : fibPair n = (Draga.fib n,Draga.fib (n+1)) := by
  induction n with
  | zero => rfl
  | succ n h => simp [fibPair,h,Draga.fib,Nat.add_comm]

/-- Concrete finite minimax guarantee. Both the impossibility below 90 bits
and attainment at 90 use precisely the same real-valued readout contract. -/
theorem ninety_bits (b : Nat) : BitFits 128 0 b ↔ 90 ≤ b := by
  rw [bit_threshold (zero_start 128) b]
  have f : Draga.fib (128+2) = 659034621587630041982498215 := by
    rw [← show (fibPair 130).1 = Draga.fib (128+2) by rw [fibPair_spec]]
    decide
  rw [f]
  constructor
  · intro cap
    by_contra h
    have upper := Nat.pow_le_pow_right (n:=2) (by decide) (show b ≤ 89 by omega)
    norm_num at upper
    omega
  · intro h
    exact le_trans (by decide : 659034621587630041982498215 ≤ 2^90)
      (Nat.pow_le_pow_right (by decide) h)

/-- The least width is the base-two ceiling logarithm, including powers of two. -/
theorem minimum_bits {n : Nat} {x : ℝ} (hx : leafUnion n x) (b : Nat) :
    BitFits n x b ↔ Nat.clog 2 (Draga.fib (n+2)) ≤ b :=
  (bit_threshold hx b).trans (Nat.le_pow_iff_clog_le (by decide))

/-- A finite cover by arbitrarily translated, closed real intervals of length
theta n. Left endpoints have no grid restriction. -/
def Covers (K : Set ℝ) (n m : Nat) : Prop :=
  ∃ origin : Fin m → ℝ, ∀ x ∈ K, ∃ i, origin i ≤ x ∧ x ≤ origin i + theta n

/-- The hypotheses are precisely nonemptiness, inclusion in the unit interval
and K = tau_0(K) ∪ tau_1(K), written without set-image coercions. The usual
compact attractor satisfies these. Neither closedness nor uniqueness of K is
needed for the covering-number conclusion. -/
structure InvariantSet (K : Set ℝ) : Prop where
  nonempty : ∃ x, x ∈ K
  unit : ∀ x ∈ K, 0 ≤ x ∧ x ≤ 1
  forward : ∀ a x, x ∈ K → step a x ∈ K
  backward : ∀ x ∈ K, ∃ a z, z ∈ K ∧ x = step a z

theorem invariant_union {K : Set ℝ} (stable : InvariantSet K) (n : Nat) :
    ∀ x ∈ K, leafUnion n x := by
  induction n using Nat.twoStepInduction with
  | zero =>
    intro x hx
    exact ⟨[],by simp [words],stable.unit x hx⟩
  | one =>
    intro x hx
    obtain ⟨a,z,hz,rfl⟩ := stable.backward x hx
    have h := stable.unit z hz
    refine ⟨[a],by cases a <;> simp [words],?_⟩
    exact ⟨step_mono a h.1,step_mono a h.2⟩
  | more n h k =>
    intro x hx
    obtain ⟨a,z,hz,rfl⟩ := stable.backward x hx
    cases a with
    | false =>
      obtain ⟨w,hw,lo,hi⟩ := h z hz
      exact ⟨false::w,List.mem_append_left _ (List.mem_map.mpr ⟨w,hw,rfl⟩),
             step_mono false lo,step_mono false hi⟩
    | true =>
      obtain ⟨w,hw,lo,hi⟩ := k z hz
      exact ⟨true::w,List.mem_append_right _ (List.mem_map.mpr ⟨w,hw,rfl⟩),
             step_mono true lo,step_mono true hi⟩

theorem invariant_value {K : Set ℝ} (stable : InvariantSet K) (w : Word)
    {x : ℝ} (hx : x ∈ K) : value w x ∈ K := by
  induction w with
  | nil => exact hx
  | cons a w ih => exact stable.forward a _ ih

/-- Full covering-number conclusion of Lemma 3.13.1. The lower bound uses
actual points in K obtained from a common starting point. It is therefore
valid even for invariant sets which omit an extremal endpoint. -/
theorem minimum_cover {K : Set ℝ} (stable : InvariantSet K) (n : Nat) :
    Covers K n (Draga.fib (n+2)) ∧ ∀ m, Covers K n m → Draga.fib (n+2) ≤ m := by
  constructor
  · rw [← word_count]
    refine ⟨fun i => leftEnd (address i),?_⟩
    intro x hx
    obtain ⟨w,hw,lo,hi⟩ := invariant_union stable n x hx
    refine ⟨indexOf w hw,?_⟩
    change leftEnd (address (indexOf w hw)) ≤ x ∧ x ≤ leftEnd (address (indexOf w hw)) + theta n
    rw [index_spec]
    have wb := leaf_width n w hw
    unfold width at wb
    refine ⟨lo,?_⟩
    rcases wb with e|e <;> linarith [theta_pos n]
  · rintro m ⟨origin,coverage⟩
    obtain ⟨x,hx⟩ := stable.nonempty
    let y := fun i : Index n => value (witnessHistory n (address i)).reverse x
    have hits (i : Index n) := coverage (y i) (invariant_value stable _ hx)
    let label := fun i : Index n => Classical.choose (hits i)
    have bounds (i : Index n) : origin (label i) ≤ y i ∧ y i ≤ origin (label i) + theta n :=
      Classical.choose_spec (hits i)
    have different {i j : Index n} (lt : i < j) : label i ≠ label j := by
      intro e
      have a := (bounds i).1
      rw [e] at a
      have b := (bounds j).2
      have gap : y i + theta n < y j := witness_gap n lt (stable.unit x hx)
      linarith
    have inj : Function.Injective label := by
      intro i j e
      rcases Nat.lt_trichotomy i.val j.val with lt|eq|gt
      · exact False.elim (different lt e)
      · exact Fin.ext eq
      · exact False.elim (different gt e.symm)
    have bound := Fintype.card_le_of_injective label inj
    simpa [Index,word_count] using bound

/-- A word is a leaf precisely when it reaches the threshold and trimming
retains the whole word. This links the recursive enumeration to the update. -/
theorem stopped_word (n : Nat) (w : Word) :
    w ∈ words n ↔ n ≤ cost w ∧ trim n w = w := by
  constructor
  · intro hw
    refine ⟨(leaf_cost_length _ _ hw).1,?_⟩
    induction n using Nat.twoStepInduction generalizing w with
    | zero => simp [words] at hw; simp [hw,trim]
    | one => simp [words] at hw; rcases hw with rfl|rfl <;> simp [trim]
    | more n h k =>
      simp only [words,List.mem_append,List.mem_map] at hw
      rcases hw with ⟨u,hu,rfl⟩ | ⟨u,hu,rfl⟩
      · simp [trim,h u hu]
      · simp [trim,k u hu]
  · rintro ⟨hc,ht⟩
    simpa [ht] using trim_leaf n w hc

theorem warm_value (j : Nat) (z : ℝ) :
    value (List.replicate j true) z = 1-theta j*(1-z) := by
  induction j with
  | zero => simp [theta,value]
  | succ j h => rw [List.replicate_succ,value_cons,h,theta_step]; unfold step; ring

/-- The free parameter j is needed when twice the allowed error is close to
the endpoint gap. Padding first keeps all histories at a common time. -/
def padded (n j : Nat) (w : Word) : Word :=
  (List.replicate (n-w.length) true ++ List.replicate j true) ++ w.reverse

theorem padded_time {n : Nat} {w : Word} (hw : w ∈ words n) (j : Nat) :
    (padded n j w).length = n+j := by
  have l := (leaf_cost_length n w hw).2
  simp [padded]; omega

theorem padded_near (n j : Nat) (w : Word) {x : ℝ} (hx : 0 ≤ x ∧ x ≤ 1) :
    rightEnd w-theta j ≤ value (padded n j w).reverse x ∧
      value (padded n j w).reverse x ≤ rightEnd w := by
  let z := value (List.replicate (n-w.length) true) x
  let y := value (List.replicate j true) z
  have eq : value (padded n j w).reverse x = value w y := by
    simp only [padded,List.reverse_append,List.reverse_reverse,List.reverse_replicate,value_append,y,z]
  have hz : 0 ≤ z ∧ z ≤ 1 := value_unit _ hx
  have hy : 0 ≤ y ∧ y ≤ 1 := value_unit _ hz
  have near : 1-theta j ≤ y := by
    dsimp only [y]; rw [warm_value]
    nlinarith [mul_nonneg (le_of_lt (theta_pos j)) hz.1]
  have lo : 0 ≤ leftEnd w := (value_unit w (x:=0) (by norm_num)).1
  have hi : rightEnd w ≤ 1 := (value_unit w (x:=1) (by norm_num)).2
  have wd : width w ≤ 1 := by unfold width; linarith
  have delta := mul_nonneg (sub_nonneg.mpr wd) (sub_nonneg.mpr hy.2)
  have e : rightEnd w = leftEnd w + width w := by unfold width; ring
  rw [eq]
  refine ⟨?_,value_mono w hy.2⟩
  rw [e,interpolation]; nlinarith

/-- Strictness at three quarters is essential: as the error approaches the
boundary, the finite warm-up needed to separate the witnesses can increase. -/
theorem below_three_quarters {n m : Nat} {ε x : ℝ} (hx : 0 ≤ x ∧ x ≤ 1)
    (threshold : ε < 3*theta n/4) (machine : FitsAt ε x m) : Draga.fib (n+2) ≤ m := by
  obtain ⟨j,hj⟩ := exists_pow_lt_of_lt_one (x:=3*theta n/2-2*ε) (y:=(1/2:ℝ))
    (by linarith) (by norm_num)
  change theta j < _ at hj
  apply finite_packing (fun i : Index n => padded n j (address i)) _ machine
  intro i k hik
  have g := List.pairwise_iff_get.mp (endpoint_spacing n) i k hik
  change rightEnd (address i) + 3*theta n/2 ≤ rightEnd (address k) at g
  have a := (padded_near n j (address i) hx).2
  have b := (padded_near n j (address k) hx).1
  linarith

theorem relax_accuracy {ε δ x : ℝ} {m : Nat} (le : ε ≤ δ) (machine : FitsAt ε x m) :
    FitsAt δ x m := by
  obtain ⟨start,update,read,h⟩ := machine
  exact ⟨start,update,read,fun w => le_trans (h w) le⟩

/-- Proposition 3.13.6. Write ell=n+1 so the coarser Fibonacci index is
unambiguous at ell=1. This also holds for startup in the whole leaf union. -/
theorem nondyadic_bounds {n : Nat} {ε x : ℝ} (hx : leafUnion (n+1) x)
    (lo : theta (n+1)/2 ≤ ε) (hi : ε < theta (n+1)) :
    FitsAt ε x (Draga.fib (n+3)) ∧
    (ε < 3*theta (n+1)/4 → ∀ m, FitsAt ε x m → Draga.fib (n+3) ≤ m) ∧
    (3*theta (n+1)/4 ≤ ε → ∀ m, FitsAt ε x m → Draga.fib (n+2) ≤ m) := by
  refine ⟨relax_accuracy lo (machine_upper _ _ hx),?_,?_⟩
  · intro gap m fit
    exact below_three_quarters (n:=n+1) (union_unit hx) gap fit
  · intro _ m fit
    apply below_three_quarters (union_unit hx) _ fit
    rw [theta_step] at hi
    linarith [theta_pos n]


/-- Executable subtree size using a rolling Fibonacci pair. -/
def count (n : Nat) : Nat := (fibPair (n+2)).1
theorem count_correct (n : Nat) : count n = Draga.fib (n+2) := by simp [count,fibPair_spec]
theorem count_split (n : Nat) : count (n+2) = count n + count (n+1) := by
  simp only [count_correct]
  exact Nat.add_comm _ _

/-- Lexicographic rank. One branch decision and at most one subtree-count
addition per letter. The recursive computation never constructs the leaf list. -/
def rankWord : Nat → Word → Nat
  | 0, _ => 0
  | 1, true::_ => 1
  | 1, _ => 0
  | _+2, [] => 0
  | n+2, false::w => rankWord n w
  | n+2, true::w => count n + rankWord (n+1) w

/-- Inverse on the initial interval of valid ordinals. -/
def unrankWord : Nat → Nat → Word
  | 0, _ => []
  | 1, i => if i = 0 then [false] else [true]
  | n+2, i => if i < count n then false :: unrankWord n i
              else true :: unrankWord (n+1) (i-count n)

/-- Both membership of the unranked leaf and recovery of its ordinal. -/
theorem ordinal_decode (n : Nat) : ∀ i < count n, unrankWord n i ∈ words n ∧
    rankWord n (unrankWord n i) = i := by
  induction n using Nat.twoStepInduction with
  | zero =>
    intro i hi
    change i < 1 at hi
    have : i=0 := by omega
    subst i; decide
  | one =>
    intro i hi
    change i < 2 at hi
    have : i=0 ∨ i=1 := by omega
    rcases this with rfl|rfl <;> decide
  | more n h k =>
    intro i hi
    rw [count_split] at hi
    by_cases lt : i < count n
    · obtain ⟨mem,rank⟩ := h i lt
      rw [unrankWord,if_pos lt]
      exact ⟨List.mem_append_left _ (List.mem_map.mpr ⟨_,mem,rfl⟩),rank⟩
    · have sub : i-count n < count (n+1) := by omega
      obtain ⟨mem,rank⟩ := k _ sub
      rw [unrankWord,if_neg lt]
      refine ⟨List.mem_append_right _ (List.mem_map.mpr ⟨_,mem,rfl⟩),?_⟩
      change count n + rankWord (n+1) _ = i
      rw [rank]; omega

/-- Every enumerated leaf has an in-range, lossless persistent ordinal. -/
theorem ordinal_encode (n : Nat) : ∀ w ∈ words n, rankWord n w < count n ∧
    unrankWord n (rankWord n w) = w := by
  induction n using Nat.twoStepInduction with
  | zero => intro w hw; simp [words] at hw; subst w; decide
  | one => intro w hw; simp [words] at hw; rcases hw with rfl|rfl <;> decide
  | more n h k =>
    intro w hw
    simp only [words,List.mem_append,List.mem_map] at hw
    rw [count_split]
    rcases hw with ⟨u,hu,rfl⟩ | ⟨u,hu,rfl⟩
    · obtain ⟨lt,eq⟩ := h _ hu
      change rankWord n u < _ ∧ unrankWord (n+2) (rankWord n u) = _
      exact ⟨by omega,by rw [unrankWord,if_pos lt,eq]⟩
    · obtain ⟨lt,eq⟩ := k _ hu
      change count n + rankWord (n+1) u < _ ∧ unrankWord (n+2) (count n + rankWord (n+1) u) = _
      refine ⟨by omega,?_⟩
      have ge : ¬ count n + rankWord (n+1) u < count n := by omega
      rw [unrankWord,if_neg ge,Nat.add_sub_cancel_left,eq]

/-- A total, table-free update on finite ordinals. Its type certifies the
persistent range. Arithmetic scratch space is a distinct resource. -/
def ordinalStep (n : Nat) (i : Fin (count n)) (a : Bool) : Fin (count n) :=
  let w := trim n (a::unrankWord n i.val)
  ⟨rankWord n w, (ordinal_encode n w (trim_update_leaf (ordinal_decode n i.val i.isLt).1 a)).1⟩

theorem ordinal_step_correct (n : Nat) (i : Fin (count n)) (a : Bool) :
    unrankWord n (ordinalStep n i a).val = trim n (a::unrankWord n i.val) :=
  (ordinal_encode n _ (trim_update_leaf (ordinal_decode n i.val i.isLt).1 a)).2

/-- Store a dyadic cylinder as (integer offset, denominator). Its affine map
is x ↦ (offset + x) / denominator. Both values are computed exactly. -/
def cylinder : Word → Nat × Nat
  | [] => (0,1)
  | false :: w => let z := cylinder w; (z.1,4*z.2)
  | true :: w => let z := cylinder w; (z.1+z.2,2*z.2)

theorem cylinder_scale (w : Word) : (cylinder w).2 = 2 ^ cost w := by
  induction w with
  | nil => rfl
  | cons a w ih => cases a <;> simp [cylinder,cost,pow_add,ih]

theorem cylinder_value (w : Word) (x : ℝ) :
    ((cylinder w).2 : ℝ) * value w x = (cylinder w).1 + x := by
  induction w with
  | nil => simp [cylinder,value]
  | cons a w ih =>
    cases a <;> simp only [cylinder,Nat.cast_mul,Nat.cast_add,value_cons,step] <;>
      push_cast <;> linarith

/-- Decoder returns an exact rational pair, including when the threshold is
overshot by one cost unit. No rounding of the readout is assumed. -/
def midpointCode (n : Nat) (i : Fin (count n)) : Nat × Nat :=
  let z := cylinder (unrankWord n i.val)
  (2*z.1+1,2*z.2)

theorem midpoint_exact (n : Nat) (i : Fin (count n)) :
    ((midpointCode n i).1 : ℝ) / (midpointCode n i).2 =
      (leftEnd (unrankWord n i.val)+rightEnd (unrankWord n i.val))/2 := by
  let w := unrankWord n i.val
  have den : 0 < (cylinder w).2 := by rw [cylinder_scale]; exact Nat.pow_pos (by decide)
  have dz : (0:ℝ) < (cylinder w).2 := by exact_mod_cast den
  have l := cylinder_value w 0
  have r := cylinder_value w 1
  change ((2*(cylinder w).1+1:Nat):ℝ) / ((2*(cylinder w).2:Nat):ℝ) =
    (value w 0+value w 1)/2
  push_cast
  apply (div_eq_iff (by linarith : (2:ℝ)*(cylinder w).2 ≠ 0)).2
  linarith

/-- All-history accuracy of the actual integer encoder, updater and rational
decoder. Initialization is by any leaf containing the prescribed x, including
endpoints. The base contract allows the same tolerance at time zero. -/
theorem ordinal_accuracy {n : Nat} (i : Fin (count n)) {x : ℝ}
    (hx : leftEnd (unrankWord n i.val) ≤ x ∧ x ≤ rightEnd (unrankWord n i.val))
    (h : Word) :
    let result := midpointCode n (h.foldl (ordinalStep n) i)
    |(result.1 : ℝ) / result.2 - value h.reverse x| ≤ theta n/2 := by
  let finish := h.foldl (ordinalStep n) i
  let w := unrankWord n finish.val
  obtain ⟨l,r⟩ := history_interval n (fun s : Fin (count n) => unrankWord n s.val)
    (ordinalStep n) (ordinal_step_correct n) h i x hx
  change leftEnd w ≤ value h.reverse x at l
  change value h.reverse x ≤ rightEnd w at r
  have hw : w ∈ words n := (ordinal_decode n finish.val finish.isLt).1
  have bound := leaf_width n w hw
  unfold width at bound
  change |((midpointCode n finish).1 : ℝ) / (midpointCode n finish).2 - _| ≤ _
  rw [midpoint_exact]
  change |(leftEnd w+rightEnd w)/2 - _| ≤ _
  rcases bound with eq|eq <;> apply abs_le.mpr <;> constructor <;> linarith [theta_pos n]

end FibonacciMemory
