//! QueryMolecule: typed representation of a SMARTS pattern graph.
//!
//! A `QueryMolecule` is a graph whose nodes carry `AtomQuery` conditions and
//! whose edges carry `BondQuery` conditions. The VF2 matcher walks this graph
//! and checks conditions against a target `Molecule`.

/// A single primitive condition on an atom.
#[derive(Debug, Clone, PartialEq)]
pub enum AtomPrimitive {
    /// `[#N]` — matches by atomic number.
    AtomicNum(u8),
    /// `C`, `N`, `O`, … — matches by element symbol.
    Symbol(String),
    /// `[a]` (true) or `[A]` (false) — aromatic / aliphatic.
    Aromatic(bool),
    /// `[+1]`, `[-1]`, … — formal charge.
    Charge(i8),
    /// `[H2]`, `[H]` — total hydrogen count (explicit + implicit).
    HCount(u8),
    /// `[h2]`, `[h]` — implicit hydrogen count only (not explicit H).
    ImplicitHCount(u8),
    /// `[D3]` — heavy-atom degree (number of bonds).
    Degree(u8),
    /// `[R]` (true) or `[!R]` (false) — ring membership.
    RingMembership(bool),
    /// `[k5]` — atom belongs to *any* perceived ring of exactly N members. This is
    /// distinct from `[rN]`/[`MinRingSize`] — confirmed against real RDKit: for a
    /// fusion atom shared between a 5-ring and a 6-ring, RDKit's `[k6]` matches but
    /// `[r6]` does not (see `docs/rdkit_compat.md`'s "SMARTS-R0"/"SMARTS-R1" entries).
    RingSize(u8),
    /// `[r5]` — the *smallest* perceived ring containing this atom has exactly N
    /// members (RDKit's actual `[rN]` semantics — not the same predicate as `[kN]`/
    /// [`RingSize`], despite chematic previously aliasing them to the same query).
    MinRingSize(u8),
    /// `*` — wildcard; matches any atom.
    Wildcard,
    /// `$(smarts)` — recursive SMARTS: the atom must be the root of a match for `smarts`.
    Recursive(Box<QueryMolecule>),
    /// `[vN]` — total valence: sum of explicit bond orders plus implicit H count.
    Valence(u8),
    /// `[xN]` — ring-bond count: number of bonds where both endpoints share a SSSR ring.
    RingBondCount(u8),
    /// `[^N]` — hybridization: 1 = sp, 2 = sp2, 3 = sp3.
    Hybridization(u8),
    /// `[XN]` — total connectivity: heavy-atom degree + implicit H count.
    TotalConnectivity(u8),
    /// `[RN]` — ring membership count: the number of SSSR rings containing this atom.
    /// `[R0]` = not in any ring; `[R1]` = in exactly 1 ring; etc.
    RingCount(u8),
    /// `[13C]`, `[2H]` — specific isotope (mass number).
    ///
    /// Whether this constraint is enforced depends on `MatchConfig::use_isotopes`.
    Isotope(u16),
    /// `[@]` (counterclockwise, value 1) or `[@@]` (clockwise, value 2) —
    /// tetrahedral chirality specification.
    ///
    /// Whether this constraint is enforced depends on `MatchConfig::use_chirality`.
    Chirality(u8),
}

/// Logical combination of atom primitives.
#[derive(Debug, Clone, PartialEq)]
pub enum AtomQuery {
    Primitive(AtomPrimitive),
    /// High-precedence AND (`&`) or juxtaposition inside brackets.
    And(Box<AtomQuery>, Box<AtomQuery>),
    /// OR (`,`).
    Or(Box<AtomQuery>, Box<AtomQuery>),
    /// NOT (`!`).
    Not(Box<AtomQuery>),
}

/// A single primitive condition on a bond.
#[derive(Debug, Clone, PartialEq)]
pub enum BondPrimitive {
    /// `-` single bond (also matches `/` and `\` stereo bonds).
    Single,
    /// `=` double bond.
    Double,
    /// `#` triple bond.
    Triple,
    /// `:` aromatic bond.
    Aromatic,
    /// `~` any bond.
    Any,
    /// `@` ring bond (both endpoints share at least one ring).
    Ring,
    /// `/` cis-like geometric direction (for E/Z double bonds).
    Up,
    /// `\` trans-like geometric direction (for E/Z double bonds).
    Down,
}

/// Logical combination of bond primitives.
#[derive(Debug, Clone, PartialEq)]
pub enum BondQuery {
    Primitive(BondPrimitive),
    And(Box<BondQuery>, Box<BondQuery>),
    Or(Box<BondQuery>, Box<BondQuery>),
    Not(Box<BondQuery>),
    /// Implicit bond: no bond specification between two atoms in SMARTS — matches any bond (`~`).
    Any,
}

/// A node in a `QueryMolecule` graph.
#[derive(Debug, Clone, PartialEq)]
pub struct QueryAtom {
    pub query: AtomQuery,
    /// Atom map number from `:N` syntax (e.g. `[O:3]`).
    /// Stored as metadata only; never used as a matching criterion.
    pub atom_map: Option<u16>,
}

/// An edge in a `QueryMolecule` graph.
#[derive(Debug, Clone, PartialEq)]
pub struct QueryBond {
    /// Index into `QueryMolecule::atoms` for one endpoint.
    pub atom1: usize,
    /// Index into `QueryMolecule::atoms` for the other endpoint.
    pub atom2: usize,
    pub query: BondQuery,
}

/// A query molecule built from a SMARTS string.
///
/// Stores atoms, bonds, and a per-atom adjacency list for fast neighbour lookup
/// during VF2 matching.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct QueryMolecule {
    pub atoms: Vec<QueryAtom>,
    pub bonds: Vec<QueryBond>,
    /// `adj[i]` = `Vec<(bond_idx, neighbour_atom_idx)>`
    pub adj: Vec<Vec<(usize, usize)>>,
}

impl QueryMolecule {
    /// Create an empty `QueryMolecule`.
    pub fn new() -> Self {
        Self::default()
    }

    /// Add an atom with the given query condition. Returns the atom index.
    pub fn add_atom(&mut self, query: AtomQuery) -> usize {
        let idx = self.atoms.len();
        self.atoms.push(QueryAtom {
            query,
            atom_map: None,
        });
        self.adj.push(vec![]);
        idx
    }

    /// Add an atom with an optional atom map number (`:N` from SMARTS/SMIRKS).
    pub fn add_atom_with_map(&mut self, query: AtomQuery, atom_map: Option<u16>) -> usize {
        let idx = self.atoms.len();
        self.atoms.push(QueryAtom { query, atom_map });
        self.adj.push(vec![]);
        idx
    }

    /// Add an undirected bond between atoms `a` and `b` with the given query condition.
    pub fn add_bond(&mut self, a: usize, b: usize, query: BondQuery) {
        let bidx = self.bonds.len();
        self.bonds.push(QueryBond {
            atom1: a,
            atom2: b,
            query,
        });
        self.adj[a].push((bidx, b));
        self.adj[b].push((bidx, a));
    }

    /// Number of atoms in this query molecule.
    pub fn atom_count(&self) -> usize {
        self.atoms.len()
    }
}
