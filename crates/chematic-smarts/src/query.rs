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
    /// `[H2]`, `[H]` — implicit hydrogen count.
    HCount(u8),
    /// `[D3]` — heavy-atom degree (number of bonds).
    Degree(u8),
    /// `[R]` (true) or `[!R]` (false) — ring membership.
    RingMembership(bool),
    /// `[r5]` — atom is in a ring of exactly N members.
    RingSize(u8),
    /// `*` — wildcard; matches any atom.
    Wildcard,
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
#[derive(Debug, Clone)]
pub struct QueryAtom {
    pub query: AtomQuery,
}

/// An edge in a `QueryMolecule` graph.
#[derive(Debug, Clone)]
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
#[derive(Debug, Clone, Default)]
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
        self.atoms.push(QueryAtom { query });
        self.adj.push(vec![]);
        idx
    }

    /// Add an undirected bond between atoms `a` and `b` with the given query condition.
    pub fn add_bond(&mut self, a: usize, b: usize, query: BondQuery) {
        let bidx = self.bonds.len();
        self.bonds.push(QueryBond { atom1: a, atom2: b, query });
        self.adj[a].push((bidx, b));
        self.adj[b].push((bidx, a));
    }
}
