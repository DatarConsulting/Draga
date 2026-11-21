import FibonacciMemory

/-! Command-line driver for the verified ordinal codec.
Each input line is "threshold initial_ordinal chronological_bits". Use "-"
for the empty history. One output line is emitted at startup and after each
transition: ordinal, exact midpoint numerator, exact midpoint denominator. -/
open FibonacciMemory

def printState (n : Nat) (state : Fin (count n)) : IO Unit := do
  let (p,q) := midpointCode n state
  IO.println s!"{state.val} {p} {q}"

def runLine (line : String) : IO Unit := do
  let [threshold, ordinal, bits] := line.splitOn " "
    | throw <| IO.userError "Expected threshold, ordinal and binary word"
  let some n := threshold.toNat? | throw <| IO.userError "Invalid threshold"
  let some i := ordinal.toNat? | throw <| IO.userError "Invalid ordinal"
  if hi : i < count n then
    let mut state : Fin (count n) := ⟨i,hi⟩
    printState n state
    if bits != "-" then
      for bit in bits.toList do
        let a ← match bit with
          | '0' => pure false
          | '1' => pure true
          | _ => throw <| IO.userError "Invalid letter"
        state := ordinalStep n state a
        printState n state
  else throw <| IO.userError "Ordinal exceeds the leaf count"

def main : IO Unit := do
  let input ← IO.getStdin
  repeat
    let line ← input.getLine
    if line.isEmpty then break
    runLine line.trim
