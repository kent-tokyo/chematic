import Kekule.Search
/-!
# Differential checker

Reads cases exported by the Rust harness (`formal/kekule/harness`) on stdin and
judges every Rust `kekulize` result against the Lean specification:

* `checkKekule` is proved equal to the specification (`checkKekule_iff`);
* `search` is proved to find an assignment iff one exists (`search_none_iff`).

Input, one block per case:

    c <id>
    e a0 b0 a1 b1 ...      aromatic bonds, in Rust bond order
    m i j k ...            must-match atoms
    x i ...                atoms a Rust heuristic may leave unmatched
    r ok k0 k1 ...         Rust Ok: indices (into the `e` list) of Double bonds
    r err                  Rust Err
-/
open Kekule

structure Case where
  id : String := ""
  edges : Array Edge := #[]
  must : List Nat := []
  excl : List Nat := []
  result : Option (List Nat) := none

def nums (ws : List String) : List Nat := ws.filterMap String.toNat?

def pairs : List Nat → List Edge
  | a :: b :: rest => (a, b) :: pairs rest
  | _ => []

inductive Verdict
  | pass | passRelaxed | passErr | heuristicDeviation
  | invalidOutput | completenessBug | badInput
  deriving BEq, Repr

def Verdict.isBug : Verdict → Bool
  | .invalidOutput | .completenessBug | .badInput => true
  | _ => false

def judge (c : Case) : Verdict := Id.run do
  let S := c.must
  if !(decide S.Nodup) then return .badInput
  let E := restrict c.edges.toList S
  let feasible := (search E S).isSome
  match c.result with
  | none => return if feasible then .completenessBug else .passErr
  | some idxs =>
    if idxs.any (· ≥ c.edges.size) then return .badInput
    let D := idxs.map (fun i => c.edges[i]!)
    if checkKekule E S D then return .pass
    let relaxed := S.filter (fun v => !c.excl.contains v)
    if checkKekule E relaxed D then
      return if feasible then .heuristicDeviation else .passRelaxed
    return .invalidOutput

def parseLine (c : Case) (line : String) : Case :=
  match line.trimAscii.toString.splitOn " " with
  | "c" :: id :: _ => { c with id := id }
  | "e" :: rest => { c with edges := (pairs (nums rest)).toArray }
  | "m" :: rest => { c with must := nums rest }
  | "x" :: rest => { c with excl := nums rest }
  | "r" :: "ok" :: rest => { c with result := some (nums rest) }
  | "r" :: "err" :: _ => { c with result := none }
  | _ => c

partial def loop (h : IO.FS.Stream) (cur : Case) (counts : Array Nat) : IO (Array Nat) := do
  let line ← h.getLine
  if line.isEmpty then return counts
  if line.trimAscii.toString == "end" then
    let v := judge cur
    let k := match v with
      | .pass => 0 | .passRelaxed => 1 | .passErr => 2 | .heuristicDeviation => 3
      | .invalidOutput => 4 | .completenessBug => 5 | .badInput => 6
    if v.isBug || v == .heuristicDeviation then
      IO.println s!"{repr v} {cur.id}"
    loop h {} (counts.modify k (· + 1))
  else
    loop h (parseLine cur line) counts

def main : IO UInt32 := do
  let counts ← loop (← IO.getStdin) {} #[0, 0, 0, 0, 0, 0, 0]
  let names := #["pass", "pass_relaxed", "pass_err", "heuristic_deviation",
    "INVALID_OUTPUT", "COMPLETENESS_BUG", "BAD_INPUT"]
  for i in [0:names.size] do
    IO.println s!"{names[i]!}: {counts[i]!}"
  return if counts[4]! + counts[5]! + counts[6]! > 0 then 1 else 0
