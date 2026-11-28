/-!
# Kekulé assignment: specification

An aromatic system is a graph whose vertices are atom indices and whose edges
are the aromatic bonds. Some atoms *must* receive a double bond (the Rust
`atom_must_be_matched` rule decides which; that chemistry is an input here,
not part of the model).

A Kekulé assignment for a must-set `S` is a set `D` of edges such that

* every edge of `D` is an aromatic bond of the graph,
* no two edges of `D` share an atom (each atom gets at most one double bond),
* every atom of `S` is touched by some edge of `D`.

This file states that specification as a `Prop` and as a `Bool` checker, and
proves the two agree, so the executable checker in `Main.lean` checks exactly
the specification.
-/

namespace Kekule

/-- An undirected edge between two atom indices. -/
abbrev Edge := Nat × Nat

/-- `e` touches atom `v`. -/
def touches (e : Edge) (v : Nat) : Bool := e.1 == v || e.2 == v

/-- Edges `e` and `f` share no atom. -/
def disj (e f : Edge) : Bool := !(touches f e.1) && !(touches f e.2)

/-- The endpoint of `e` that is not `s` (for an edge touching `s`). -/
def other (e : Edge) (s : Nat) : Nat := if e.1 = s then e.2 else e.1

/-- The specification, as a proposition. -/
structure IsKekule (E : List Edge) (S : List Nat) (D : List Edge) : Prop where
  subset : ∀ e ∈ D, e ∈ E
  matching : D.Pairwise (fun e f => disj e f = true)
  covers : ∀ v ∈ S, ∃ e ∈ D, touches e v = true

/-- The same specification as an executable check. -/
def checkKekule (E : List Edge) (S : List Nat) (D : List Edge) : Bool :=
  D.all (fun e => E.contains e) &&
  pairwiseB D &&
  S.all (fun v => D.any (fun e => touches e v))
where
  pairwiseB : List Edge → Bool
    | [] => true
    | e :: rest => rest.all (fun f => disj e f) && pairwiseB rest

theorem pairwiseB_iff (D : List Edge) :
    checkKekule.pairwiseB D = true ↔ D.Pairwise (fun e f => disj e f = true) := by
  induction D with
  | nil => simp [checkKekule.pairwiseB]
  | cons e rest ih => simp [checkKekule.pairwiseB, ih, List.pairwise_cons]

/-- The executable checker decides exactly the specification. -/
theorem checkKekule_iff (E : List Edge) (S : List Nat) (D : List Edge) :
    checkKekule E S D = true ↔ IsKekule E S D := by
  constructor
  · intro h
    simp only [checkKekule, Bool.and_eq_true, List.all_eq_true, List.any_eq_true,
      List.contains_iff_mem, pairwiseB_iff] at h
    exact ⟨h.1.1, h.1.2, h.2⟩
  · intro ⟨hs, hm, hc⟩
    simp only [checkKekule, Bool.and_eq_true, List.all_eq_true, List.any_eq_true,
      List.contains_iff_mem, pairwiseB_iff]
    exact ⟨⟨hs, hm⟩, hc⟩

theorem disj_self (e : Edge) : disj e e = false := by
  simp [disj, touches]

theorem disj_comm (e f : Edge) : disj e f = disj f e := by
  obtain ⟨a, b⟩ := e
  obtain ⟨c, d⟩ := f
  simp only [disj, touches]
  rw [Bool.eq_iff_iff]
  simp only [Bool.and_eq_true, Bool.not_eq_true', Bool.or_eq_false_iff, beq_eq_false_iff_ne]
  omega

end Kekule
