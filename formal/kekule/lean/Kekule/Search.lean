import Kekule.Spec
/-!
# Reference model: exhaustive search for a Kekulé assignment

`search E S` looks for a Kekulé assignment of the must-set `S` using only edges
of `E`. It is deliberately simple (backtracking, exponential in the worst case)
so that its correctness is easy to prove; it is a reference to test the Rust
augmenting-path + blossom matcher against, not a replacement for it.

Proved here:

* `search_sound`: anything `search` returns is a valid Kekulé assignment.
* `search_complete`: if any Kekulé assignment exists (with edges inside `S`),
  `search` finds one.

Together: `search E S` returns `none` **iff** no Kekulé assignment exists, so a
Rust `Err` on an input where `search` returns `some` is a real completeness bug.
-/

namespace Kekule

def search (E : List Edge) : List Nat → Option (List Edge)
  | [] => some []
  | s :: rest =>
    (E.filter (fun e => touches e s)).findSome? fun e =>
      (search (E.filter (fun f => disj e f)) (rest.erase (other e s))).map (e :: ·)
termination_by S => S.length
decreasing_by
  have := List.length_erase_le (l := rest) (a := other e s)
  simp; omega

theorem touches_other {e : Edge} {s : Nat} (h : touches e s = true) :
    touches e (other e s) = true := by
  obtain ⟨a, b⟩ := e
  simp [touches, other] at *
  by_cases ha : a = s <;> simp_all

theorem touches_iff {e : Edge} {s v : Nat} (h : touches e s = true)
    (hv : touches e v = true) : v = s ∨ v = other e s := by
  obtain ⟨a, b⟩ := e
  simp only [touches, other, Bool.or_eq_true, beq_iff_eq] at *
  by_cases ha : a = s <;> simp only [ha, if_true, if_false] <;> omega

/-- An edge disjoint from `e` avoids both endpoints `s` and `other e s` of `e`. -/
theorem disj_avoid {e f : Edge} {s v : Nat} (hes : touches e s = true)
    (h : disj e f = true) (hv : touches f v = true) : v ≠ s ∧ v ≠ other e s := by
  obtain ⟨a, b⟩ := e
  obtain ⟨c, d⟩ := f
  simp only [disj, touches, other, Bool.and_eq_true, Bool.not_eq_true',
    Bool.or_eq_false_iff, beq_eq_false_iff_ne, Bool.or_eq_true, beq_iff_eq] at *
  by_cases ha : a = s <;> simp only [ha, if_true, if_false] <;> omega

theorem findSome?_some {α β} {f : α → Option β} {l : List α} {b : β}
    (h : l.findSome? f = some b) : ∃ a ∈ l, f a = some b := by
  induction l with
  | nil => simp at h
  | cons x xs ih =>
    simp only [List.findSome?_cons] at h
    split at h
    · rename_i b' hb; cases h; exact ⟨x, by simp, hb⟩
    · obtain ⟨a, ha, hfa⟩ := ih h; exact ⟨a, by simp [ha], hfa⟩

theorem findSome?_isSome {α β} {f : α → Option β} {l : List α} {a : α}
    (ha : a ∈ l) (hf : (f a).isSome) : (l.findSome? f).isSome := by
  induction l with
  | nil => simp at ha
  | cons x xs ih =>
    simp only [List.findSome?_cons]
    split
    · simp
    · rename_i hx
      rcases List.mem_cons.mp ha with rfl | h
      · simp [hx] at hf
      · exact ih h

/-- Soundness: a result of `search` is a Kekulé assignment. -/
theorem search_sound (E : List Edge) (S : List Nat) (D : List Edge)
    (h : search E S = some D) : IsKekule E S D := by
  match S with
  | [] => simp [search] at h; subst h; exact ⟨by simp, by simp, by simp⟩
  | s :: rest =>
    rw [search] at h
    obtain ⟨e, he, hres⟩ := findSome?_some h
    rw [List.mem_filter] at he
    obtain ⟨heE, hes⟩ := he
    cases hrec : search (E.filter (fun f => disj e f)) (rest.erase (other e s)) with
    | none => simp [hrec] at hres
    | some D' =>
      simp [hrec] at hres; subst hres
      have ih := search_sound _ _ D' hrec
      refine ⟨?_, ?_, ?_⟩
      · intro f hf
        rcases List.mem_cons.mp hf with rfl | hf
        · exact heE
        · exact (List.mem_filter.mp (ih.subset f hf)).1
      · refine List.Pairwise.cons ?_ ih.matching
        intro f hf; exact (List.mem_filter.mp (ih.subset f hf)).2
      · intro v hv
        rcases List.mem_cons.mp hv with rfl | hv
        · exact ⟨e, by simp, hes⟩
        · by_cases hvt : v = other e s
          · exact ⟨e, by simp, hvt ▸ touches_other hes⟩
          · obtain ⟨f, hf, hfv⟩ := ih.covers v ((List.mem_erase_of_ne hvt).mpr hv)
            exact ⟨f, List.mem_cons_of_mem _ hf, hfv⟩
termination_by S.length
decreasing_by
  have := List.length_erase_le (l := rest) (a := other e s)
  simp; omega

/-- Well-formed input: no self-loops, and every edge lies inside `S`. -/
structure WF (E : List Edge) (S : List Nat) : Prop where
  nodup : S.Nodup
  noLoop : ∀ e ∈ E, e.1 ≠ e.2
  inside : ∀ e ∈ E, e.1 ∈ S ∧ e.2 ∈ S

/-- Completeness: if a Kekulé assignment exists, `search` finds one. -/
theorem search_complete (E : List Edge) (S : List Nat) (D : List Edge)
    (wf : WF E S) (hD : IsKekule E S D) : (search E S).isSome := by
  match S with
  | [] => simp [search]
  | s :: rest =>
    rw [search]
    obtain ⟨e, heD, hes⟩ := hD.covers s (by simp)
    have heE := hD.subset e heD
    refine findSome?_isSome (a := e) (List.mem_filter.mpr ⟨heE, hes⟩) ?_
    let t := other e s
    let E' := E.filter (fun f => disj e f)
    let D' := D.filter (fun f => f != e)
    have hsr : s ∉ rest := (List.nodup_cons.mp wf.nodup).1
    have hnd : rest.Nodup := (List.nodup_cons.mp wf.nodup).2
    -- every other edge of D is disjoint from e
    have hdisj : ∀ f ∈ D, f ≠ e → disj e f = true := by
      have hR : D.Pairwise (fun a b => a ≠ b → disj a b = true) :=
        hD.matching.imp (fun {a b} h (_ : a ≠ b) => h)
      have hflip : D.Pairwise (flip (fun a b => a ≠ b → disj a b = true)) :=
        hD.matching.imp (fun {a b} h (_ : b ≠ a) => by rw [disj_comm]; exact h)
      have hall := List.Pairwise.forall_of_forall_of_flip
        (fun x _ h => absurd rfl h) hR hflip
      intro f hf hne
      exact hall heD hf (Ne.symm hne)
    have herase : ∀ v, v ∈ rest.erase t ↔ v ∈ rest ∧ v ≠ t := by
      intro v
      rw [hnd.erase_eq_filter, List.mem_filter]
      simp
    have wf' : WF E' (rest.erase t) := by
      refine ⟨hnd.erase t, fun f hf => wf.noLoop f (List.mem_filter.mp hf).1, ?_⟩
      intro f hf
      obtain ⟨hfE, hfe⟩ := List.mem_filter.mp hf
      obtain ⟨h1, h2⟩ := wf.inside f hfE
      have a1 := disj_avoid hes hfe (v := f.1) (by simp [touches])
      have a2 := disj_avoid hes hfe (v := f.2) (by simp [touches])
      simp only [List.mem_cons] at h1 h2
      refine ⟨(herase _).mpr ⟨?_, a1.2⟩, (herase _).mpr ⟨?_, a2.2⟩⟩
      · rcases h1 with h | h
        · exact absurd h a1.1
        · exact h
      · rcases h2 with h | h
        · exact absurd h a2.1
        · exact h
    have hD' : IsKekule E' (rest.erase t) D' := by
      refine ⟨?_, hD.matching.filter _, ?_⟩
      · intro f hf
        obtain ⟨hfD, hne⟩ := List.mem_filter.mp hf
        have hne : f ≠ e := by simpa using hne
        exact List.mem_filter.mpr ⟨hD.subset f hfD, hdisj f hfD hne⟩
      · intro v hv
        obtain ⟨hvr, hvt⟩ := (herase v).mp hv
        have hvs : v ≠ s := fun h => hsr (h ▸ hvr)
        obtain ⟨f, hfD, hfv⟩ := hD.covers v (List.mem_cons_of_mem _ hvr)
        refine ⟨f, List.mem_filter.mpr ⟨hfD, ?_⟩, hfv⟩
        have : f ≠ e := by
          rintro rfl
          rcases touches_iff hes hfv with h | h
          · exact hvs h
          · exact hvt h
        simpa using this
    have ih := search_complete E' (rest.erase t) D' wf' hD'
    simp only [Option.isSome_map]
    exact ih
termination_by S.length
decreasing_by
  have := List.length_erase_le (l := rest) (a := other e s)
  simp; omega

end Kekule

namespace Kekule

/-- Keep only the aromatic edges whose two atoms are distinct and both in `S`
(the Rust matcher only ever puts a double bond between two must-match atoms). -/
def restrict (edges : List Edge) (S : List Nat) : List Edge :=
  edges.filter (fun e => e.1 != e.2 && S.contains e.1 && S.contains e.2)

theorem restrict_wf (edges : List Edge) (S : List Nat) (h : S.Nodup) :
    WF (restrict edges S) S := by
  refine ⟨h, ?_, ?_⟩
  · intro e he
    have := (List.mem_filter.mp he).2
    simp only [Bool.and_eq_true, bne_iff_ne] at this
    exact this.1.1
  · intro e he
    have := (List.mem_filter.mp he).2
    simp only [Bool.and_eq_true, List.contains_iff_mem] at this
    exact ⟨this.1.2, this.2⟩

/-- `search` decides existence: it fails exactly when no Kekulé assignment
exists. This is the fact the differential checker relies on when it flags a
Rust `Err` as a bug. -/
theorem search_none_iff (edges : List Edge) (S : List Nat) (h : S.Nodup) :
    search (restrict edges S) S = none ↔ ¬ ∃ D, IsKekule (restrict edges S) S D := by
  constructor
  · intro hn ⟨D, hD⟩
    have := search_complete _ _ D (restrict_wf edges S h) hD
    simp [hn] at this
  · intro hne
    cases hs : search (restrict edges S) S with
    | none => rfl
    | some D => exact absurd ⟨D, search_sound _ _ D hs⟩ hne

end Kekule
