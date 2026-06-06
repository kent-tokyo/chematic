/* tslint:disable */
/* eslint-disable */

/**
 * A conformer ensemble: one molecule geometry with multiple 3D coordinate sets.
 *
 * Create with `new(smiles)`, then add conformers with `add_generated_conformer`
 * or `add_minimized_conformer`.  Retrieve coordinates as PDB strings via
 * `get_conformer_pdb(idx)`.  Compare conformers with `conformer_rmsd`.
 */
export class ConformerHandle {
    free(): void;
    [Symbol.dispose](): void;
    /**
     * Generate a new 3D conformer using distance-geometry and add it to the ensemble.
     *
     * Returns the index of the newly added conformer.
     */
    add_generated_conformer(): number;
    /**
     * Generate a new 3D conformer, run force-field minimization, and add it.
     *
     * Returns the index of the newly added conformer.
     */
    add_minimized_conformer(): number;
    /**
     * Number of conformers currently stored.
     */
    conformer_count(): number;
    /**
     * Kabsch-aligned RMSD (Å) between conformers `a` and `b`.
     *
     * Returns `NaN` if either index is out of range.
     */
    conformer_rmsd(a: number, b: number): number;
    /**
     * Un-aligned (translation + rotation NOT removed) RMSD (Å) between conformers `a` and `b`.
     *
     * Returns `NaN` if either index is out of range.
     */
    conformer_rmsd_no_align(a: number, b: number): number;
    /**
     * Return conformer `idx` as a PDB string, or `null` if `idx` is out of range.
     */
    get_conformer_pdb(idx: number): string | undefined;
    /**
     * The ensemble's molecule as a `MolHandle`.
     */
    mol(): MolHandle;
    /**
     * Create a new empty ensemble for the molecule given by `smiles`.
     *
     * Returns a JS error on SMILES parse failure.
     */
    constructor(smiles: string);
    /**
     * Remove conformer `idx` and return `true`, or `false` if `idx` is out of range.
     */
    remove_conformer(idx: number): boolean;
}

/**
 * Style options for [`MolHandle::depict_svg_opts`].
 *
 * Construct with `new DepictOptions()`, then call setters:
 * ```js
 * const opts = new DepictOptions();
 * opts.set_background("transparent");
 * opts.set_dark(true);
 * opts.set_width(240);
 * opts.set_height(240);
 * ```
 */
export class DepictOptions {
    free(): void;
    [Symbol.dispose](): void;
    constructor();
    /**
     * Set a per-atom color override (CSS color string).  Calling multiple times
     * for the same `idx` uses the last value.  The atom is highlighted even if
     * not in `set_highlight_atoms`.
     */
    set_atom_color(idx: number, color: string): void;
    set_atom_ids(v: boolean): void;
    set_background(bg: string): void;
    set_dark(dark: boolean): void;
    set_height(h: number): void;
    set_highlight_atoms(atoms: Uint32Array): void;
    set_highlight_bonds(bonds: Uint32Array): void;
    set_highlight_color(color: string): void;
    set_kekulize(v: boolean): void;
    set_padding(p: number): void;
    set_show_atom_indices(v: boolean): void;
    set_width(w: number): void;
}

/**
 * A handle to a parsed molecule.  Owns the molecule behind an `Rc` so that
 * it can be cheaply cloned on the JS side without copying atom/bond data.
 */
export class MolHandle {
    private constructor();
    free(): void;
    [Symbol.dispose](): void;
    /**
     * Number of aromatic rings (all ring atoms aromatic).
     */
    aromatic_ring_count(): number;
    /**
     * Number of heavy atoms (explicit atoms in the graph; does not count implicit H).
     */
    atom_count(): number;
    /**
     * Bertz complexity index (BertzCT).
     */
    bertz_ct(): number;
    /**
     * Number of bonds.
     */
    bond_count(): number;
    /**
     * Canonical SMILES string.
     */
    canonical_smiles(): string;
    /**
     * Kier–Hall χ0 molecular connectivity index.
     */
    chi0(): number;
    /**
     * Kier–Hall χ0v valence-weighted connectivity index.
     */
    chi0v(): number;
    /**
     * Kier–Hall χ1 molecular connectivity index.
     */
    chi1(): number;
    /**
     * Kier–Hall χ1v valence-weighted connectivity index.
     */
    chi1v(): number;
    /**
     * Kier–Hall χ2 molecular connectivity index.
     */
    chi2(): number;
    /**
     * Kier–Hall χ2v valence-weighted connectivity index.
     */
    chi2v(): number;
    /**
     * Kier–Hall χ3 molecular connectivity index.
     */
    chi3(): number;
    /**
     * Kier–Hall χ3v valence-weighted connectivity index.
     */
    chi3v(): number;
    /**
     * Kier–Hall χ4 molecular connectivity index.
     */
    chi4(): number;
    /**
     * Kier–Hall χ4v valence-weighted connectivity index.
     */
    chi4v(): number;
    /**
     * 2D SVG depiction of the molecule (CPK coloring).
     */
    depict_svg(): string;
    /**
     * 2D SVG depiction with style options.
     */
    depict_svg_opts(opts: DepictOptions): string;
    /**
     * Returns `true` if the molecule passes Egan's absorption criteria
     * (TPSA ≤ 131.6 Å² and LogP ≤ 5.88).
     */
    egan_passes(): boolean;
    /**
     * Monoisotopic (exact) mass.
     */
    exact_mass(): number;
    /**
     * Sum of formal charges.
     */
    formal_charge_sum(): number;
    /**
     * Molecular formula string (Hill notation: C first, H second, then alphabetical).
     */
    formula(): string;
    /**
     * Fraction of sp3 carbons (Fsp3).
     */
    fsp3(): number;
    /**
     * Returns `true` if the molecule passes Ghose's drug-likeness filter
     * (MW 160–480, LogP −0.4–5.6, HeavyAtoms 20–70, MR 40–130).
     */
    ghose_passes(): boolean;
    /**
     * Number of hydrogen bond acceptors (Lipinski: all N and O atoms).
     */
    hba_count(): number;
    /**
     * Number of hydrogen bond donors (N-H or O-H groups).
     */
    hbd_count(): number;
    /**
     * Number of non-hydrogen heavy atoms.
     */
    heavy_atom_count(): number;
    /**
     * Hall–Kier κ1 shape index.
     */
    kappa1(): number;
    /**
     * Hall–Kier κ2 shape index.
     */
    kappa2(): number;
    /**
     * Hall–Kier κ3 shape index.
     */
    kappa3(): number;
    /**
     * Labute approximate surface area (Å²).
     */
    labute_asa(): number;
    /**
     * Returns `true` if the molecule satisfies Lipinski's Rule of Five.
     */
    lipinski_passes(): boolean;
    /**
     * Crippen–Wildman octanol/water partition coefficient (LogP).
     */
    logp_crippen(): number;
    /**
     * Maximum EState index across all heavy atoms.
     */
    max_estate(): number;
    /**
     * Minimum EState index across all heavy atoms.
     */
    min_estate(): number;
    /**
     * Wildman–Crippen molar refractivity (MR).
     */
    molar_refractivity(): number;
    /**
     * Average molecular weight (Da).
     */
    molecular_weight(): number;
    /**
     * Morgan count fingerprint as a JSON object string (`{"<hash>": count, …}`).
     *
     * `radius` controls the ECFP radius (2 = ECFP4-equivalent).
     */
    morgan_fp_counts_json(radius: number): string;
    /**
     * Number of non-aromatic rings containing at least one heteroatom.
     */
    num_aliphatic_heterocycles(): number;
    /**
     * Count of aliphatic (non-aromatic) rings in the SSSR.
     */
    num_aliphatic_rings(): number;
    /**
     * Number of aromatic rings containing at least one heteroatom (N, O, S, …).
     */
    num_aromatic_heterocycles(): number;
    /**
     * Number of bridgehead atoms (shared by ≥2 rings with ≥3 ring bonds).
     */
    num_bridgehead_atoms(): number;
    /**
     * Number of heteroatoms (non-C, non-H heavy atoms).
     */
    num_heteroatoms(): number;
    /**
     * Number of fully saturated rings containing at least one heteroatom.
     */
    num_saturated_heterocycles(): number;
    /**
     * Count of fully saturated rings in the SSSR.
     */
    num_saturated_rings(): number;
    /**
     * Number of spiro atoms (sole shared atom between exactly 2 rings).
     */
    num_spiro_atoms(): number;
    /**
     * Number of assigned stereocenters (R/S).
     */
    num_stereocenters(): number;
    /**
     * Count of tetrahedral stereocenters with unspecified configuration.
     */
    num_unspecified_stereocenters(): number;
    /**
     * Returns `true` if the molecule has no PAINS structural alerts.
     */
    pains_passes(): boolean;
    /**
     * Quantitative Estimate of Drug-likeness (QED); range [0, 1].
     */
    qed(): number;
    /**
     * Returns `true` if the molecule passes the REOS (Rapid Elimination Of Swill) filter.
     */
    reos_passes(): boolean;
    /**
     * Total number of rings (SSSR count).
     */
    ring_count(): number;
    /**
     * Number of rotatable bonds.
     */
    rotatable_bond_count(): number;
    /**
     * Sum of EState indices over all heavy atoms.
     */
    sum_estate(): number;
    /**
     * Topological polar surface area (Å²).
     */
    tpsa(): number;
    /**
     * Returns `true` if the molecule passes Veber's oral bioavailability criteria
     * (TPSA ≤ 140 Å² and rotatable bonds ≤ 10).
     */
    veber_passes(): boolean;
    /**
     * Wiener topological index (sum of all pairwise shortest-path distances).
     */
    wiener_index(): number;
}

/**
 * Return a copy of the molecule with all implicit hydrogens converted to explicit H atoms.
 */
export function add_hydrogens(mol: MolHandle): MolHandle;

/**
 * AtomPair fingerprint as a bit-packed byte vector (256 bytes = 2048 bits).
 */
export function atom_pair_bitvec(mol: MolHandle): Uint8Array;

/**
 * Check whether a reaction SMILES is atom-balanced.
 *
 * Returns JSON: `{ "balanced": true|false, "diff": ["C: 1 reactant vs 2 product", ...] }`
 * Returns `"error:<msg>"` on parse failure.
 */
export function balance_check_json(reaction_smiles: string): string;

/**
 * Number of BRICS fragments produced by fragmenting the molecule.
 *
 * Returns 1 if no BRICS-breakable bonds exist (whole molecule is one fragment).
 */
export function brics_fragment_count(mol: MolHandle): number;

/**
 * BRICS fragment SMILES as a JSON array.
 *
 * Applies the BRICS fragmentation rules and returns the canonical SMILES of
 * every resulting fragment.  Returns `[]` for molecules with no BRICS-breakable
 * bonds (e.g. benzene).
 *
 * The count of fragments equals `brics_fragment_count`.
 */
export function brics_fragments_json(mol: MolHandle): string;

/**
 * Cluster molecules by structural similarity (Butina algorithm, ECFP4 Tanimoto).
 *
 * `smiles_json` — a JSON array of SMILES strings.
 * `cutoff` — Tanimoto similarity threshold (0.0–1.0); molecules within this
 *   distance of a cluster centre are assigned to that cluster.
 * Returns a JSON array of clusters, each cluster being an array of 0-based input indices.
 * Returns a JS error if any SMILES fails to parse.
 */
export function butina_cluster_ecfp4_json(smiles_json: string, cutoff: number): string;

/**
 * Canonical tautomer of `mol`.
 *
 * Applies a rule-based tautomer normalisation and returns the canonical form
 * as a new `MolHandle`.
 */
export function canonical_tautomer(mol: MolHandle): MolHandle;

/**
 * Parse all molecular fragments from a CDXML string.
 *
 * Returns a JSON array of SMILES strings, one per fragment:
 * `["CC","c1ccccc1"]`
 *
 * Stereochemistry (wedge/dash bonds) is read from the `Display` attribute
 * of bond elements.
 */
export function cdxml_to_smiles_json(cdxml: string): string;

/**
 * CIP stereo assignments as a JSON array of `{atomIdx, cipCode}` objects.
 *
 * `cipCode` is one of `"R"`, `"S"`, `"E"`, or `"Z"`.
 * Returns `[]` for molecules with no specified stereocenters.
 */
export function cip_assignments_json(mol: MolHandle): string;

/**
 * Compute direct Coulomb energy for a molecule with Gasteiger partial charges.
 *
 * Returns JSON object: `{ "coulomb_energy": E, "unit": "kcal/mol" }`
 *
 * # Arguments
 * * `mol` - Molecule to evaluate
 *
 * # Example (JavaScript)
 * ```js
 * const mol = parse_smiles("CCO");
 * const result = coulomb_energy_json(mol);
 * // { "coulomb_energy": -12.34, "unit": "kcal/mol" }
 * ```
 */
export function coulomb_energy_json(mol: MolHandle): string;

/**
 * Return the CPK color (CSS hex string) for the given element symbol.
 *
 * Returns `"#000000"` (black) for carbon and unknown elements.
 */
export function cpk_color(element_symbol: string): string;

/**
 * Compute structured depiction data for `mol` as a JSON object.
 *
 * Returns:
 * ```json
 * {
 *   "atoms": [
 *     {"idx": 0, "element": "C", "x": 1.5, "y": 0.0, "charge": 0,
 *      "label": null, "color": "#000000"},
 *     ...
 *   ],
 *   "bonds": [
 *     {"idx": 0, "atom1": 0, "atom2": 1, "kind": "Single"},
 *     ...
 *   ]
 * }
 * ```
 *
 * `label` is `null` for carbon atoms in skeletal structures (label suppressed).
 * `kind` is one of `"Single"`, `"Double"`, `"Triple"`, `"Aromatic"`, `"Up"`, `"Down"`.
 */
export function depict_data_json(mol: MolHandle): string;

/**
 * Compute structured depiction data using caller-supplied 2D coordinates.
 *
 * `coords_json` — JSON array of `[x, y]` pairs, one per atom in order.
 *
 * Returns the same JSON format as `depict_data_json`.
 */
export function depict_data_with_coords_json(mol: MolHandle, coords_json: string): string;

/**
 * Render a reaction SMILES string (e.g. `"CC(=O)O.CCO>>CC(=O)OCC.O"`) as a
 * single SVG showing reactants → products with `+` separators.
 *
 * Returns a self-contained SVG string.  Returns a JS error on invalid input.
 */
export function depict_reaction_svg(rxn_smiles: string): string;

/**
 * Render a grid SVG from newline-separated SMILES (one per line).
 *
 * Lines that fail to parse are silently skipped.
 * `cols` controls the number of columns (each cell is 200×200 px).
 */
export function depict_svg_grid(smiles_block: string, cols: number): string;

/**
 * Render a molecule grid with SMARTS-based atom highlighting.
 *
 * `smiles_block` — newline-separated SMILES strings (same format as `depict_svg_grid`).
 * `cols` — number of grid columns.
 * `match_smarts` — SMARTS pattern; matched atoms in each molecule are highlighted.
 *   Pass an empty string `""` to render without any highlighting.
 *
 * Invalid SMILES are rendered as empty cells; SMARTS parse failure returns an
 * unhighlighted grid (the SMARTS is silently ignored).
 */
export function depict_svg_grid_highlighted(smiles_block: string, cols: number, match_smarts: string): string;

/**
 * Detect named functional groups in `mol`.
 *
 * Returns a JSON array of `{"name":"hydroxyl","atoms":[3]}` objects.
 * Multiple matches of the same group (e.g. two hydroxyl groups) each appear
 * as a separate entry.  Overlapping groups (carboxylic acid → "carboxyl" +
 * "hydroxyl" + "carbonyl") are all returned.
 */
export function detect_functional_groups(mol: MolHandle): string;

/**
 * Dice similarity between `a` and `b` using ECFP4 fingerprints.
 */
export function dice_ecfp4(a: MolHandle, b: MolHandle): number;

/**
 * Dice similarity between `a` and `b` using ECFP6 fingerprints.
 */
export function dice_ecfp6(a: MolHandle, b: MolHandle): number;

/**
 * Dice similarity between `a` and `b` using MACCS 166-bit fingerprints.
 */
export function dice_maccs(a: MolHandle, b: MolHandle): number;

/**
 * Compute the ECFP4 fingerprint as a bit-packed byte vector (256 bytes = 2048 bits).
 */
export function ecfp4_bitvec(mol: MolHandle): Uint8Array;

/**
 * ECFP6 (radius-3) fingerprint as a bit-packed byte vector (256 bytes = 2048 bits).
 */
export function ecfp6_bitvec(mol: MolHandle): Uint8Array;

/**
 * Compute a fingerprint bit-vector with configurable ECFP radius and bit width.
 *
 * `radius` — Morgan radius (1 = ECFP2, 2 = ECFP4, 3 = ECFP6).
 * `nbits` — bit width; must be one of 256, 512, 1024, or 2048.
 *   Returns a `Uint8Array` of `nbits/8` bytes.
 *
 * The hash modulo is applied at fingerprint-generation time (`id % nbits`),
 * so no post-processing fold is needed.
 */
export function ecfp_bitvec_custom(mol: MolHandle, radius: number, nbits: number): Uint8Array;

/**
 * Enumerate all stereoisomers arising from unspecified tetrahedral stereocenters.
 *
 * Only considers carbon stereocenters without explicit `@`/`@@` annotation.
 * Already-specified centers and E/Z double-bond geometry are unchanged.
 * Returns a JSON array of canonical SMILES strings.
 *
 * At most 2^6 = 64 combinations are enumerated; if more than 6 unspecified
 * centers are present this function returns a JS error to avoid combinatorial
 * explosion.
 */
export function enumerate_stereo_isomers_json(mol: MolHandle): string;

/**
 * All enumerated tautomers of `mol` as a JSON array of canonical SMILES strings.
 *
 * Example return value: `["Oc1cccc2ccccc12","O=C1C=CC=Cc2ccccc21"]`
 */
export function enumerate_tautomers_json(mol: MolHandle): string;

/**
 * Per-atom EState values as a JSON array of f64.
 *
 * Indices match `mol.atoms()` order.  Hydrogen atoms get 0.0.
 */
export function estate_indices_json(mol: MolHandle): string;

/**
 * FCFP4 (pharmacophore, radius-2) fingerprint as a bit-packed byte vector (256 bytes).
 */
export function fcfp4_bitvec(mol: MolHandle): Uint8Array;

/**
 * FCFP6 (pharmacophore, radius-3) fingerprint as a bit-packed byte vector (256 bytes).
 */
export function fcfp6_bitvec(mol: MolHandle): Uint8Array;

/**
 * Analyze a reaction SMILES and return the reaction center as JSON.
 *
 * JSON schema: `{ broken: [[a1,a2],...], formed: [[a1,a2],...], changed: [a,...] }`
 * where atom indices are 0-based within the first reactant molecule.
 * Returns an error string prefixed with `"error:"` on failure.
 */
export function find_reaction_center_json(reaction_smiles: string): string;

/**
 * Gasteiger-Marsili PEOE partial charges as a JSON array of f64.
 */
export function gasteiger_charges_json(mol: MolHandle): string;

/**
 * Generate energy-minimized 3D coordinates and return a PDB string.
 *
 * Runs distance-geometry placement followed by gradient-descent force-field
 * minimization.  Geometry quality is better than `generate_3d_pdb` for
 * flexible molecules; the force field is approximate (not MMFF94/UFF).
 */
export function generate_3d_minimized_pdb(mol: MolHandle): string;

/**
 * Generate 3D coordinates for the molecule and return a PDB string.
 *
 * Coordinates are generated using distance-geometry placement with ring templates.
 * Returns heavy-atom PDB (HETATM records, no explicit H).
 */
export function generate_3d_pdb(mol: MolHandle): string;

/**
 * Generic (atom-type-erased) Murcko scaffold of `mol`.
 *
 * All atoms become carbon and all bonds become single bonds, giving the pure
 * graph topology of the scaffold.
 */
export function generic_murcko_scaffold(mol: MolHandle): MolHandle;

/**
 * Return information about a single atom as a JSON object.
 *
 * `idx` is the 0-based atom index (matching `atoms()` order).
 * Returns `"null"` if `idx` is out of range.
 *
 * Fields: `element` (symbol), `hybridization` ("sp"/"sp2"/"sp3"),
 * `charge` (formal charge integer), `isAromatic` (bool),
 * `totalHydrogens` (explicit + implicit H count, integer).
 * sp3d/sp3d2 (hypervalent P/S) are not distinguished from sp3/sp2.
 */
export function get_atom_info(mol: MolHandle, idx: number): string;

/**
 * Return bond information as a JSON object, looked up by the two bonded atom indices.
 *
 * Useful when you know the atom indices from SMARTS matching or `data-atom-idx` SVG
 * attributes but not the bond index.  Returns `"null"` if no bond exists between them.
 *
 * Fields: same as `get_bond_info` plus `bondIdx` (u32).
 */
export function get_bond_between(mol: MolHandle, atom1: number, atom2: number): string;

/**
 * Return bond information as a JSON object, looked up by bond index.
 *
 * `idx` is the 0-based bond index (order matches `mol.bonds()` iteration).
 * Returns `"null"` if `idx` is out of range.
 *
 * Fields: `bondOrder` (1.0/1.5/2.0/3.0), `isAromatic` (bool),
 * `isInRing` (bool), `atomFrom` (u32), `atomTo` (u32).
 */
export function get_bond_info(mol: MolHandle, idx: number): string;

/**
 * All scalar molecular descriptors as a single JSON object.
 *
 * Keys use camelCase and match the individual `MolHandle` method names.
 * Drug-likeness rule outcomes are included as boolean fields.
 */
export function get_descriptors_json(mol: MolHandle): string;

/**
 * Identify functional groups. Returns a JSON array of objects:
 * `[{"atoms":[0,2,3],"type":"C,N,O"}, …]`
 */
export function identify_functional_groups(mol: MolHandle): string;

/**
 * Returns `true` if the SMILES string can be parsed without error.
 */
export function is_valid_smiles(s: string): boolean;

/**
 * Per-atom Labute approximate surface area contributions as a JSON array of f64.
 *
 * Non-finite values (single-atom molecules etc.) are emitted as JSON `null`.
 */
export function labute_asa_per_atom_json(mol: MolHandle): string;

/**
 * Return the largest fragment of `mol` (salt/solvent stripping).
 *
 * For single-component molecules returns a copy of the same molecule.
 */
export function largest_fragment(mol: MolHandle): MolHandle;

/**
 * Per-atom Crippen LogP contributions as a JSON array of f64.
 *
 * Index `i` corresponds to atom `i` in `mol.atoms()` order.
 */
export function logp_per_atom_json(mol: MolHandle): string;

/**
 * MACCS 166-bit structural keys fingerprint as a byte array (21 bytes, LSB-first).
 *
 * Bit `i` (0-indexed) corresponds to MACCS key `i+1`.
 */
export function maccs_bitvec(mol: MolHandle): Uint8Array;

/**
 * Find all SMARTS matches in a molecule given only SMILES strings.
 *
 * Convenience wrapper around `smarts_match_atoms` that accepts raw SMILES
 * instead of a `MolHandle`.  Returns the same JSON format: `[[0,1],[3,4]]`.
 * Returns a JS error on SMILES or SMARTS parse failure.
 */
export function match_smarts_smiles(smiles: string, smarts: string): string;

/**
 * Select `n` maximally-diverse molecules (MaxMin algorithm, ECFP4 Tanimoto).
 *
 * `smiles_json` — a JSON array of SMILES strings, e.g. `["CC","c1ccccc1","CCO"]`.
 * Returns a JSON array of 0-based indices into the input array.
 * Returns a JS error if any SMILES fails to parse (indices would otherwise shift).
 */
export function maxmin_picks_ecfp4_json(smiles_json: string, n: number): string;

/**
 * Maximum Common Substructure of a set of molecules, returned as a canonical SMILES string.
 *
 * `smiles_json` — a JSON array of at least 2 SMILES strings.
 * Returns the MCS SMILES, or `"null"` when no common substructure was found.
 * Returns a JS error on SMILES parse failure.
 */
export function mcs_smiles_json(smiles_json: string): string;

/**
 * Find matched molecular pairs in a set of molecules as JSON.
 *
 * `smiles_json` — JSON array of SMILES strings to analyze.
 *
 * Returns a JSON array of matched pairs:
 * ```json
 * [
 *   {
 *     "mol_a": "CC(=O)Oc1ccccc1",
 *     "mol_b": "CC(=O)Nc1ccccc1",
 *     "core": "c1ccccc1[*]",
 *     "fragment_a": "[*]OC(C)=O",
 *     "fragment_b": "[*]NC(C)=O"
 *   }
 * ]
 * ```
 *
 * Each pair represents molecules that share a common core scaffold but differ
 * by exactly one structural fragment at a single BRICS-breakable bond cut.
 *
 * Returns a JS error if any SMILES fails to parse.
 */
export function mmp_pairs_json(smiles_json: string): string;

/**
 * Parse a Tripos MOL2 string and return SMILES.
 *
 * Returns `"error:<msg>"` on failure.
 */
export function mol2_to_smiles(mol2_str: string): string;

/**
 * Parse a MOL V2000 string and return 2D coordinates as a JSON array.
 *
 * Returns `[[x0,y0],[x1,y1],...]` in atom-insertion order.
 * Coordinates are in Ångström as stored in the MOL file.
 */
export function mol_block_coords_json(mol_block: string): string;

/**
 * Serialize a SMILES string directly to a MOL V2000 block with 2D coordinates.
 *
 * Returns a JS error on SMILES parse failure.
 */
export function mol_block_from_smiles(smiles: string): string;

/**
 * Parse a ChemDraw XML (CDXML) string into a `MolHandle`.
 *
 * Only the first molecular fragment in the document is returned.
 * Returns a JS error if the document cannot be parsed.
 */
export function mol_from_cdxml(cdxml: string): MolHandle;

/**
 * Parse a CML string into a `MolHandle`.
 *
 * Returns a JS error if the CML is invalid (unknown element, bad bond, etc.).
 */
export function mol_from_cml(cml: string): MolHandle;

/**
 * Parse a PDB file and return a `MolHandle` (topology only; coordinates are discarded).
 *
 * Uses CONECT records for connectivity if present; otherwise infers bonds from
 * atom distances (the same heuristic as the internal `pdb_to_molecule` function).
 */
export function mol_from_pdb(pdb: string): MolHandle;

/**
 * Parse a MOL V2000 block and return a `MolHandle`.
 *
 * Returns a JS error string on parse failure.
 */
export function mol_from_sdf_block(block: string): MolHandle;

/**
 * Parse a MOL V3000 block and return a `MolHandle`.
 *
 * Returns a JS error string on parse failure.
 */
export function mol_from_v3000_block(block: string): MolHandle;

/**
 * Parse an XYZ file and return a `MolHandle` (topology only; coordinates are discarded).
 *
 * Returns a JS error on parse failure.
 */
export function mol_from_xyz(xyz: string): MolHandle;

/**
 * Return the index that would be assigned to an atom appended to `mol`.
 */
export function mol_next_atom_idx(mol: MolHandle): number;

/**
 * Return a new `MolHandle` with one atom appended.
 *
 * The second return value is the new atom's index (as a JS number).
 * Use `with_atom_added_idx` to retrieve the index.
 */
export function mol_with_atom_added(mol: MolHandle, element_symbol: string): MolHandle;

/**
 * Return a new `MolHandle` with the formal charge of atom `idx` changed.
 *
 * Returns a JS error if `idx` is out of range.
 */
export function mol_with_atom_charge(mol: MolHandle, idx: number, charge: number): MolHandle;

/**
 * Return a new `MolHandle` with the element of atom `idx` changed.
 *
 * `element_symbol` — periodic-table symbol, e.g. `"N"`, `"O"`, `"Cl"`.
 * Returns a JS error if `idx` is out of range or the symbol is unknown.
 */
export function mol_with_atom_element(mol: MolHandle, idx: number, element_symbol: string): MolHandle;

/**
 * Return a new `MolHandle` with atom `idx` and all its bonds removed.
 *
 * Atom indices above `idx` shift down by 1.  Returns a JS error if `idx`
 * is out of range.
 */
export function mol_with_atom_removed(mol: MolHandle, idx: number): MolHandle;

/**
 * Return a new `MolHandle` with one bond added between `a` and `b`.
 *
 * `order` — 1 = single, 2 = double, 3 = triple.
 * Returns a JS error if the bond already exists or `a == b`.
 */
export function mol_with_bond_added(mol: MolHandle, a: number, b: number, order: number): MolHandle;

/**
 * Return a new `MolHandle` with bond `idx` removed.
 *
 * Atom indices are unchanged; bond indices above `idx` shift down.
 * Returns a JS error if `idx` is out of range.
 */
export function mol_with_bond_removed(mol: MolHandle, idx: number): MolHandle;

/**
 * Per-atom molar refractivity contributions as a JSON array of f64.
 */
export function mr_per_atom_json(mol: MolHandle): string;

/**
 * Murcko scaffold of `mol` — the ring system plus linkers, side-chains removed.
 *
 * Returns a new `MolHandle`.  For acyclic molecules returns an empty molecule.
 */
export function murcko_scaffold(mol: MolHandle): MolHandle;

/**
 * Find the k nearest neighbours of a query SMILES in a list of db SMILES.
 *
 * `db_smiles_json`: JSON array of SMILES strings, e.g. `["CC","c1ccccc1"]`.
 * Returns JSON: `[{"index":0,"tanimoto":0.95},...]` sorted by descending Tanimoto.
 * Returns `"error:<msg>"` on parse failure.
 */
export function nearest_neighbors_json(query_smiles: string, db_smiles_json: string, k: number): string;

/**
 * Neutralize formal charges on `mol` by proton addition/removal.
 *
 * Returns a new `MolHandle` with all formal charges set to zero where possible.
 */
export function neutralize_charges(mol: MolHandle): MolHandle;

/**
 * Parse and re-serialise a reaction SMILES string, returning the normalised form.
 *
 * Useful for validating reaction SMILES and obtaining a canonical representation.
 * Returns a JS error on parse failure.
 */
export function normalize_reaction_smiles(rxn_smiles: string): string;

/**
 * PAINS structural alert names matched by `mol` as a JSON array.
 *
 * Returns `[]` when no alerts fire, or e.g. `["ene_six_het_A(483)"]`.
 * Use alongside `pains_passes()` to know *which* alerts triggered.
 */
export function pains_matches_json(mol: MolHandle): string;

/**
 * Parse a SMILES string into a `MolHandle`.
 *
 * Returns a JS error string on parse failure.
 */
export function parse_smiles(s: string): MolHandle;

/**
 * PEOE_VSA descriptors (14 bins) as a JSON array.
 */
export function peoe_vsa_json(mol: MolHandle): string;

/**
 * Return a copy of the molecule with all explicit hydrogen atoms removed.
 */
export function remove_hydrogens(mol: MolHandle): MolHandle;

/**
 * Decompose a set of molecules against a core SMARTS, returning R-group SMILES.
 *
 * `smiles_json` — JSON array of SMILES strings.
 * `core_smarts` — SMARTS pattern with `*` (wildcard) atoms marking R-group
 *   attachment points.  For example `c1ccc(*)cc1` for para-substituted benzene.
 *
 * Returns a JSON array with one entry per input molecule:
 * ```json
 * [
 *   {"matched":true, "r1":"C"},
 *   {"matched":true, "r1":"CC"},
 *   {"matched":false}
 * ]
 * ```
 * R-group keys are `"r1"`, `"r2"`, … in the order the `*` atoms appear in
 * the SMARTS pattern.  A molecule that does not contain the core gets
 * `"matched": false` and no R-group keys.
 *
 * Returns a JS error if the SMARTS fails to parse or any SMILES is invalid.
 */
export function rgroup_decompose_json(smiles_json: string, core_smarts: string): string;

/**
 * Run molecular dynamics simulation and return trajectory as JSON.
 *
 * Returns JSON object with trajectory frames: `{ "frames": [{ "step": N, "potential": E, "kinetic": K, "temp": T }, …] }`
 * Uses NVT ensemble (Berendsen thermostat) at 300 K by default.
 * Note: Limited to molecules with ~50 atoms or fewer for practical WASM performance.
 */
export function run_md_json(mol: MolHandle, steps: number, temp_k: number): string;

/**
 * Apply a SMIRKS reaction template and return product SMILES as a JSON string.
 *
 * `reactants_smiles`: pipe-separated SMILES, one per reactant slot in the SMIRKS.
 * Returns a JSON array of arrays: `[["product_smi", …], …]`.
 * Returns a JS error on parse failure or arity mismatch.
 */
export function run_reactants(smirks: string, reactants_smiles: string): string;

/**
 * Synthetic Accessibility Score (1 = easy, 10 = hard).
 */
export function sa_score(mol: MolHandle): number;

/**
 * Serialize multiple molecules with properties to an SDF string.
 *
 * # Arguments
 * * `smiles_json` — JSON array of SMILES strings, e.g. `["CC(=O)O","c1ccccc1"]`
 * * `names_json`  — JSON array of molecule names (same length as `smiles_json`)
 * * `props_json`  — JSON array where each element encodes one molecule's SD data fields
 *   as `"key1\tvalue1\nkey2\tvalue2"` (tab-separated key/value, `\n`-separated pairs;
 *   pass `""` for a molecule with no properties)
 *
 * Returns the SDF string, or a JS error if any SMILES fails to parse or the
 * arrays have mismatched lengths.
 *
 * The `\n` and `\t` sequences in `props_json` are JSON-escaped — they are
 * decoded to the actual characters before SDF formatting.
 */
export function sdf_from_records_json(smiles_json: string, names_json: string, props_json: string): string;

/**
 * Parse an SDF string and return a JSON array of record objects.
 *
 * Each record has the shape:
 * ```json
 * {"smiles":"CC(=O)O","name":"aspirin","properties":{"MW":"180.2","Activity":"high"}}
 * ```
 *
 * Invalid records are represented as `null`.  SD data fields are included in
 * `properties`; multi-line values are joined with `\n`.
 */
export function sdf_to_records_json(sdf: string): string;

/**
 * Parse an SDF string and return a JSON array of canonical SMILES strings.
 *
 * Invalid records are represented as `null` in the array.
 */
export function sdf_to_smiles_json(sdf: string): string;

/**
 * 3D shape descriptors as a JSON object.
 *
 * Keys: `pmi1`, `pmi2`, `pmi3`, `npr1`, `npr2`, `asphericity`, `eccentricity`,
 * `radiusOfGyration`, `planeOfBestFit`.  Non-finite values (e.g. single-atom
 * molecules where pmi3 = 0) are serialised as JSON `null`.
 */
export function shape_descriptors_json(mol: MolHandle): string;

/**
 * SlogP_VSA descriptors (12 bins) as a JSON array.
 */
export function slogp_vsa_json(mol: MolHandle): string;

/**
 * Find all substructure matches of a SMARTS pattern in `mol`.
 *
 * Returns JSON array of arrays of atom indices (sorted, 0-based).
 * Example: `[[0,1,2],[3,4,5]]` — two matches.
 * Returns `"[]"` if no match. Returns a JS error on invalid SMARTS.
 */
export function smarts_match_atoms(smarts: string, mol: MolHandle): string;

/**
 * Serialise a JSON array of SMILES to an SDF string.
 *
 * Generates 2D coordinates for each molecule.  Property data can be
 * included by using `sdf_from_records_json` instead.
 */
export function smiles_array_to_sdf(smiles_json: string): string;

/**
 * Convert a SMILES to a minimal Tripos MOL2 string (no 3D coordinates).
 *
 * Returns `"error:<msg>"` on parse failure.
 */
export function smiles_to_mol2(smiles: string): string;

/**
 * Render a highlighted SVG from a SMILES string in one call.
 *
 * `atoms` — 0-based atom indices to highlight (Uint32Array in JS).
 * `bonds` — 0-based bond indices to highlight (Uint32Array in JS).
 * `color` — CSS color for highlights (e.g. `"#ef4444"`); empty string uses default yellow.
 *
 * Returns a JS error on SMILES parse failure.
 */
export function smiles_to_svg_highlighted(smiles: string, atoms: Uint32Array, bonds: Uint32Array, color: string): string;

/**
 * SMR_VSA descriptors (10 bins) as a JSON array.
 */
export function smr_vsa_json(mol: MolHandle): string;

/**
 * Smallest Set of Smallest Rings (SSSR) as a JSON array of atom-index arrays.
 *
 * Example return value for naphthalene:
 * `[[0,1,2,3,4,5],[5,6,7,8,9,4]]`
 */
export function sssr_rings_json(mol: MolHandle): string;

/**
 * Standardize a SMILES string and return the canonical SMILES of the result.
 *
 * Applies: largest fragment extraction → charge neutralization.
 * Returns `"error:<msg>"` on parse failure.
 */
export function standardize_smiles(smiles: string): string;

export function start(): void;

/**
 * Tanimoto similarity between two molecules using AtomPair fingerprints.
 */
export function tanimoto_atom_pair(a: MolHandle, b: MolHandle): number;

/**
 * Tanimoto similarity between two molecules using ECFP4 fingerprints.
 */
export function tanimoto_ecfp4(a: MolHandle, b: MolHandle): number;

/**
 * Tanimoto similarity between `a` and `b` using ECFP6 fingerprints.
 */
export function tanimoto_ecfp6(a: MolHandle, b: MolHandle): number;

/**
 * Tanimoto similarity between two molecules using FCFP4 fingerprints (pharmacophore-based).
 */
export function tanimoto_fcfp4(a: MolHandle, b: MolHandle): number;

/**
 * Tanimoto similarity between `a` and `b` using FCFP6 (radius-3 pharmacophore) fingerprints.
 */
export function tanimoto_fcfp6(a: MolHandle, b: MolHandle): number;

/**
 * Tanimoto similarity between `a` and `b` using MACCS 166-bit fingerprints.
 */
export function tanimoto_maccs(a: MolHandle, b: MolHandle): number;

/**
 * Tanimoto similarity between two molecules given only SMILES strings (ECFP4).
 *
 * Returns a JS error on parse failure.
 */
export function tanimoto_smiles(smiles1: string, smiles2: string): number;

/**
 * Tanimoto similarity between two molecules using topological path fingerprints.
 */
export function tanimoto_topo_path(a: MolHandle, b: MolHandle): number;

/**
 * Tanimoto similarity between two molecules using Topological Torsion fingerprints.
 */
export function tanimoto_torsion(a: MolHandle, b: MolHandle): number;

/**
 * Serialise a `MolHandle` to a CML string with 2D coordinates.
 *
 * Coordinates are generated using the same 2D layout engine as `to_mol_block`.
 */
export function to_cml(mol: MolHandle): string;

/**
 * Serialize a molecule to a MOL V2000 block with 2D coordinates.
 *
 * Atom positions are computed via the same layout engine used for SVG depiction
 * and converted to Ångström units (`1.5 Å` per bond).
 */
export function to_mol_block(mol: MolHandle): string;

/**
 * Serialise a `MolHandle` to MOL V3000 format with 2D coordinates.
 */
export function to_mol_v3000_block(mol: MolHandle): string;

/**
 * Serialize a molecule to XYZ format.
 *
 * 3D coordinates are generated via distance-geometry placement.
 */
export function to_xyz(mol: MolHandle): string;

/**
 * Torsion fingerprint as a bit-packed byte vector (256 bytes = 2048 bits).
 */
export function torsion_bitvec(mol: MolHandle): Uint8Array;

/**
 * Non-canonical SMILES for `mol`.
 *
 * Unlike `canonical_smiles`, the output depends on the internal atom ordering
 * and is not normalised.  Useful when round-trip fidelity (preserving atom
 * order) matters more than a canonical form.
 */
export function write_smiles(mol: MolHandle): string;

export type InitInput = RequestInfo | URL | Response | BufferSource | WebAssembly.Module;

export interface InitOutput {
    readonly memory: WebAssembly.Memory;
    readonly __wbg_conformerhandle_free: (a: number, b: number) => void;
    readonly __wbg_depictoptions_free: (a: number, b: number) => void;
    readonly __wbg_molhandle_free: (a: number, b: number) => void;
    readonly add_hydrogens: (a: number) => number;
    readonly atom_pair_bitvec: (a: number) => [number, number];
    readonly balance_check_json: (a: number, b: number) => [number, number];
    readonly brics_fragment_count: (a: number) => number;
    readonly brics_fragments_json: (a: number) => [number, number];
    readonly butina_cluster_ecfp4_json: (a: number, b: number, c: number) => [number, number, number, number];
    readonly canonical_tautomer: (a: number) => number;
    readonly cdxml_to_smiles_json: (a: number, b: number) => [number, number, number, number];
    readonly cip_assignments_json: (a: number) => [number, number];
    readonly conformerhandle_add_generated_conformer: (a: number) => number;
    readonly conformerhandle_add_minimized_conformer: (a: number) => number;
    readonly conformerhandle_conformer_count: (a: number) => number;
    readonly conformerhandle_conformer_rmsd: (a: number, b: number, c: number) => number;
    readonly conformerhandle_conformer_rmsd_no_align: (a: number, b: number, c: number) => number;
    readonly conformerhandle_get_conformer_pdb: (a: number, b: number) => [number, number];
    readonly conformerhandle_mol: (a: number) => number;
    readonly conformerhandle_new: (a: number, b: number) => [number, number, number];
    readonly conformerhandle_remove_conformer: (a: number, b: number) => number;
    readonly coulomb_energy_json: (a: number) => [number, number];
    readonly cpk_color: (a: number, b: number) => [number, number];
    readonly depict_data_json: (a: number) => [number, number];
    readonly depict_data_with_coords_json: (a: number, b: number, c: number) => [number, number];
    readonly depict_reaction_svg: (a: number, b: number) => [number, number, number, number];
    readonly depict_svg_grid: (a: number, b: number, c: number) => [number, number];
    readonly depict_svg_grid_highlighted: (a: number, b: number, c: number, d: number, e: number) => [number, number];
    readonly depictoptions_new: () => number;
    readonly depictoptions_set_atom_color: (a: number, b: number, c: number, d: number) => void;
    readonly depictoptions_set_atom_ids: (a: number, b: number) => void;
    readonly depictoptions_set_background: (a: number, b: number, c: number) => void;
    readonly depictoptions_set_dark: (a: number, b: number) => void;
    readonly depictoptions_set_height: (a: number, b: number) => void;
    readonly depictoptions_set_highlight_atoms: (a: number, b: number, c: number) => void;
    readonly depictoptions_set_highlight_bonds: (a: number, b: number, c: number) => void;
    readonly depictoptions_set_highlight_color: (a: number, b: number, c: number) => void;
    readonly depictoptions_set_kekulize: (a: number, b: number) => void;
    readonly depictoptions_set_padding: (a: number, b: number) => void;
    readonly depictoptions_set_show_atom_indices: (a: number, b: number) => void;
    readonly depictoptions_set_width: (a: number, b: number) => void;
    readonly detect_functional_groups: (a: number) => [number, number];
    readonly dice_ecfp4: (a: number, b: number) => number;
    readonly dice_ecfp6: (a: number, b: number) => number;
    readonly dice_maccs: (a: number, b: number) => number;
    readonly ecfp4_bitvec: (a: number) => [number, number];
    readonly ecfp6_bitvec: (a: number) => [number, number];
    readonly ecfp_bitvec_custom: (a: number, b: number, c: number) => [number, number];
    readonly enumerate_stereo_isomers_json: (a: number) => [number, number, number, number];
    readonly enumerate_tautomers_json: (a: number) => [number, number];
    readonly estate_indices_json: (a: number) => [number, number];
    readonly fcfp4_bitvec: (a: number) => [number, number];
    readonly fcfp6_bitvec: (a: number) => [number, number];
    readonly find_reaction_center_json: (a: number, b: number) => [number, number];
    readonly gasteiger_charges_json: (a: number) => [number, number];
    readonly generate_3d_minimized_pdb: (a: number) => [number, number];
    readonly generate_3d_pdb: (a: number) => [number, number];
    readonly generic_murcko_scaffold: (a: number) => number;
    readonly get_atom_info: (a: number, b: number) => [number, number];
    readonly get_bond_between: (a: number, b: number, c: number) => [number, number];
    readonly get_bond_info: (a: number, b: number) => [number, number];
    readonly get_descriptors_json: (a: number) => [number, number];
    readonly identify_functional_groups: (a: number) => [number, number];
    readonly is_valid_smiles: (a: number, b: number) => number;
    readonly labute_asa_per_atom_json: (a: number) => [number, number];
    readonly largest_fragment: (a: number) => number;
    readonly logp_per_atom_json: (a: number) => [number, number];
    readonly maccs_bitvec: (a: number) => [number, number];
    readonly match_smarts_smiles: (a: number, b: number, c: number, d: number) => [number, number, number, number];
    readonly maxmin_picks_ecfp4_json: (a: number, b: number, c: number) => [number, number, number, number];
    readonly mcs_smiles_json: (a: number, b: number) => [number, number, number, number];
    readonly mmp_pairs_json: (a: number, b: number) => [number, number, number, number];
    readonly mol2_to_smiles: (a: number, b: number) => [number, number];
    readonly mol_block_coords_json: (a: number, b: number) => [number, number, number, number];
    readonly mol_block_from_smiles: (a: number, b: number) => [number, number, number, number];
    readonly mol_from_cdxml: (a: number, b: number) => [number, number, number];
    readonly mol_from_cml: (a: number, b: number) => [number, number, number];
    readonly mol_from_pdb: (a: number, b: number) => number;
    readonly mol_from_sdf_block: (a: number, b: number) => [number, number, number];
    readonly mol_from_v3000_block: (a: number, b: number) => [number, number, number];
    readonly mol_from_xyz: (a: number, b: number) => [number, number, number];
    readonly mol_next_atom_idx: (a: number) => number;
    readonly mol_with_atom_added: (a: number, b: number, c: number) => [number, number, number];
    readonly mol_with_atom_charge: (a: number, b: number, c: number) => [number, number, number];
    readonly mol_with_atom_element: (a: number, b: number, c: number, d: number) => [number, number, number];
    readonly mol_with_atom_removed: (a: number, b: number) => [number, number, number];
    readonly mol_with_bond_added: (a: number, b: number, c: number, d: number) => [number, number, number];
    readonly mol_with_bond_removed: (a: number, b: number) => [number, number, number];
    readonly molhandle_aromatic_ring_count: (a: number) => number;
    readonly molhandle_bertz_ct: (a: number) => number;
    readonly molhandle_bond_count: (a: number) => number;
    readonly molhandle_canonical_smiles: (a: number) => [number, number];
    readonly molhandle_chi0: (a: number) => number;
    readonly molhandle_chi0v: (a: number) => number;
    readonly molhandle_chi1: (a: number) => number;
    readonly molhandle_chi1v: (a: number) => number;
    readonly molhandle_chi2: (a: number) => number;
    readonly molhandle_chi2v: (a: number) => number;
    readonly molhandle_chi3: (a: number) => number;
    readonly molhandle_chi3v: (a: number) => number;
    readonly molhandle_chi4: (a: number) => number;
    readonly molhandle_chi4v: (a: number) => number;
    readonly molhandle_depict_svg: (a: number) => [number, number];
    readonly molhandle_depict_svg_opts: (a: number, b: number) => [number, number];
    readonly molhandle_egan_passes: (a: number) => number;
    readonly molhandle_exact_mass: (a: number) => number;
    readonly molhandle_formal_charge_sum: (a: number) => number;
    readonly molhandle_formula: (a: number) => [number, number];
    readonly molhandle_fsp3: (a: number) => number;
    readonly molhandle_ghose_passes: (a: number) => number;
    readonly molhandle_hba_count: (a: number) => number;
    readonly molhandle_hbd_count: (a: number) => number;
    readonly molhandle_heavy_atom_count: (a: number) => number;
    readonly molhandle_kappa1: (a: number) => number;
    readonly molhandle_kappa2: (a: number) => number;
    readonly molhandle_kappa3: (a: number) => number;
    readonly molhandle_labute_asa: (a: number) => number;
    readonly molhandle_lipinski_passes: (a: number) => number;
    readonly molhandle_logp_crippen: (a: number) => number;
    readonly molhandle_max_estate: (a: number) => number;
    readonly molhandle_min_estate: (a: number) => number;
    readonly molhandle_molar_refractivity: (a: number) => number;
    readonly molhandle_molecular_weight: (a: number) => number;
    readonly molhandle_morgan_fp_counts_json: (a: number, b: number) => [number, number];
    readonly molhandle_num_aliphatic_heterocycles: (a: number) => number;
    readonly molhandle_num_aliphatic_rings: (a: number) => number;
    readonly molhandle_num_aromatic_heterocycles: (a: number) => number;
    readonly molhandle_num_bridgehead_atoms: (a: number) => number;
    readonly molhandle_num_heteroatoms: (a: number) => number;
    readonly molhandle_num_saturated_heterocycles: (a: number) => number;
    readonly molhandle_num_saturated_rings: (a: number) => number;
    readonly molhandle_num_spiro_atoms: (a: number) => number;
    readonly molhandle_num_stereocenters: (a: number) => number;
    readonly molhandle_num_unspecified_stereocenters: (a: number) => number;
    readonly molhandle_pains_passes: (a: number) => number;
    readonly molhandle_qed: (a: number) => number;
    readonly molhandle_reos_passes: (a: number) => number;
    readonly molhandle_ring_count: (a: number) => number;
    readonly molhandle_rotatable_bond_count: (a: number) => number;
    readonly molhandle_sum_estate: (a: number) => number;
    readonly molhandle_tpsa: (a: number) => number;
    readonly molhandle_veber_passes: (a: number) => number;
    readonly molhandle_wiener_index: (a: number) => number;
    readonly mr_per_atom_json: (a: number) => [number, number];
    readonly murcko_scaffold: (a: number) => number;
    readonly nearest_neighbors_json: (a: number, b: number, c: number, d: number, e: number) => [number, number];
    readonly neutralize_charges: (a: number) => number;
    readonly normalize_reaction_smiles: (a: number, b: number) => [number, number, number, number];
    readonly pains_matches_json: (a: number) => [number, number];
    readonly parse_smiles: (a: number, b: number) => [number, number, number];
    readonly peoe_vsa_json: (a: number) => [number, number];
    readonly remove_hydrogens: (a: number) => number;
    readonly rgroup_decompose_json: (a: number, b: number, c: number, d: number) => [number, number, number, number];
    readonly run_md_json: (a: number, b: number, c: number) => [number, number];
    readonly run_reactants: (a: number, b: number, c: number, d: number) => [number, number, number, number];
    readonly sa_score: (a: number) => number;
    readonly sdf_from_records_json: (a: number, b: number, c: number, d: number, e: number, f: number) => [number, number, number, number];
    readonly sdf_to_records_json: (a: number, b: number) => [number, number];
    readonly sdf_to_smiles_json: (a: number, b: number) => [number, number];
    readonly shape_descriptors_json: (a: number) => [number, number];
    readonly slogp_vsa_json: (a: number) => [number, number];
    readonly smarts_match_atoms: (a: number, b: number, c: number) => [number, number, number, number];
    readonly smiles_array_to_sdf: (a: number, b: number) => [number, number, number, number];
    readonly smiles_to_mol2: (a: number, b: number) => [number, number];
    readonly smiles_to_svg_highlighted: (a: number, b: number, c: number, d: number, e: number, f: number, g: number, h: number) => [number, number, number, number];
    readonly smr_vsa_json: (a: number) => [number, number];
    readonly sssr_rings_json: (a: number) => [number, number];
    readonly standardize_smiles: (a: number, b: number) => [number, number];
    readonly tanimoto_atom_pair: (a: number, b: number) => number;
    readonly tanimoto_ecfp4: (a: number, b: number) => number;
    readonly tanimoto_ecfp6: (a: number, b: number) => number;
    readonly tanimoto_fcfp4: (a: number, b: number) => number;
    readonly tanimoto_fcfp6: (a: number, b: number) => number;
    readonly tanimoto_maccs: (a: number, b: number) => number;
    readonly tanimoto_smiles: (a: number, b: number, c: number, d: number) => [number, number, number];
    readonly tanimoto_topo_path: (a: number, b: number) => number;
    readonly tanimoto_torsion: (a: number, b: number) => number;
    readonly to_cml: (a: number) => [number, number];
    readonly to_mol_block: (a: number) => [number, number];
    readonly to_mol_v3000_block: (a: number) => [number, number];
    readonly to_xyz: (a: number) => [number, number];
    readonly torsion_bitvec: (a: number) => [number, number];
    readonly write_smiles: (a: number) => [number, number];
    readonly start: () => void;
    readonly molhandle_atom_count: (a: number) => number;
    readonly __wbindgen_free: (a: number, b: number, c: number) => void;
    readonly __wbindgen_malloc: (a: number, b: number) => number;
    readonly __wbindgen_realloc: (a: number, b: number, c: number, d: number) => number;
    readonly __wbindgen_externrefs: WebAssembly.Table;
    readonly __externref_table_dealloc: (a: number) => void;
    readonly __wbindgen_start: () => void;
}

export type SyncInitInput = BufferSource | WebAssembly.Module;

/**
 * Instantiates the given `module`, which can either be bytes or
 * a precompiled `WebAssembly.Module`.
 *
 * @param {{ module: SyncInitInput }} module - Passing `SyncInitInput` directly is deprecated.
 *
 * @returns {InitOutput}
 */
export function initSync(module: { module: SyncInitInput } | SyncInitInput): InitOutput;

/**
 * If `module_or_path` is {RequestInfo} or {URL}, makes a request and
 * for everything else, calls `WebAssembly.instantiate` directly.
 *
 * @param {{ module_or_path: InitInput | Promise<InitInput> }} module_or_path - Passing `InitInput` directly is deprecated.
 *
 * @returns {Promise<InitOutput>}
 */
export default function __wbg_init (module_or_path?: { module_or_path: InitInput | Promise<InitInput> } | InitInput | Promise<InitInput>): Promise<InitOutput>;
