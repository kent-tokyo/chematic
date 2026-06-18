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
     * Cluster conformers by Kabsch-aligned RMSD and return a JSON object
     * describing which conformers to keep.
     *
     * Uses greedy leader-linkage: conformers are visited in index order; each
     * is compared against existing cluster representatives. If the RMSD to any
     * representative is < `rms_threshold`, the conformer is discarded; otherwise
     * it starts a new cluster and is kept.
     *
     * Returns `{"kept_indices":[0,3,7,...],"removed_count":5}` on success.
     */
    cluster_conformers_json(rms_threshold: number): string;
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
 * MinHash LSH index: insert MHFP fingerprints and query by approximate similarity.
 *
 * ```js
 * const idx = new MhfpLshHandle(128);
 * const i0 = idx.add_smiles("c1ccccc1");    // benzene → index 0
 * const i1 = idx.add_smiles("Cc1ccccc1");   // toluene → index 1
 * const hits = JSON.parse(idx.query_json("c1ccccc1", 0.5));
 * // hits: [{index:0,similarity:1.0}, {index:1,similarity:0.xxx}]
 * ```
 */
export class MhfpLshHandle {
    free(): void;
    [Symbol.dispose](): void;
    /**
     * Add a molecule by SMILES; returns its 0-based index in the index.
     */
    add_smiles(smiles: string): number;
    /**
     * True if the index contains no molecules.
     */
    is_empty(): boolean;
    /**
     * Number of molecules in the index.
     */
    len(): number;
    /**
     * Create a new LSH index for MHFP fingerprints with `num_hashes` hash lanes.
     * Default band decomposition: 16 bands × (num_hashes / 16) rows.
     * `num_hashes` must be a multiple of 16 (e.g. 128).
     */
    constructor(num_hashes: number);
    /**
     * Query by SMILES for all entries with similarity ≥ threshold.
     *
     * Returns a JSON array `[{"index":N,"similarity":0.xxx},...]` sorted by
     * descending similarity.  Empty array `[]` when nothing qualifies.
     */
    query_json(query_smiles: string, threshold: number): string;
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
     * Assign CIP (R/S/E/Z) stereocenters and return JSON.
     *
     * Format: `{"centers":[{"atom":0,"code":"R"},{"atom":3,"code":"E"}]}`
     */
    assign_cip_json(): string;
    /**
     * Number of heavy atoms (explicit atoms in the graph; does not count implicit H).
     */
    atom_count(): number;
    /**
     * Returns true when TPSA < 90 Å², MW < 400, HBD ≤ 3.
     */
    bbb_passes(): boolean;
    /**
     * Clark (2000) blood-brain barrier logBB score.
     * logBB > −1.0 = likely CNS penetrant.
     */
    bbb_score(): number;
    /**
     * Bertz complexity index (BertzCT).
     */
    bertz_ct(): number;
    /**
     * Number of bonds.
     */
    bond_count(): number;
    /**
     * Palm (1997) Caco-2 intestinal permeability (logPCaco2).
     * > −5.5 = high permeability.
     */
    caco2_permeability(): number;
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
     * CYP3A4 metabolic inhibition risk score (0.0–1.0).
     */
    cyp3a4_inhibition_risk(): number;
    /**
     * 2D PNG depiction (rasterized from SVG).
     * Returns PNG data as base64-encoded string for embedding in HTML/JS.
     */
    depict_png(): Uint8Array;
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
     * hERG cardiac toxicity risk score (0.0–1.0).
     */
    herg_risk_score(): number;
    /**
     * Isotope distribution as JSON.
     *
     * Returns `[{"mass":100.0,"abundance":0.9},...]` sorted by mass.
     * `resolution`: m/z bin width in Da (e.g. `0.1` for nominal, `0.01` for high-res).
     */
    isotope_distribution_json(resolution: number): string;
    /**
     * Generate IUPAC systematic name for the molecule.
     *
     * Returns the name string on success, or an empty string when the
     * structure is outside the supported naming scope (complex polycyclics,
     * multi-functional groups, etc.).
     */
    iupac_name(): string;
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
     * LogD (distribution coefficient) at a specific pH.
     *
     * Accounts for ionization state: neutral molecules return LogP unchanged,
     * ionizable molecules are adjusted by log(neutral_fraction).
     */
    logd_at_ph(ph: number): number;
    /**
     * LogD profile across a pH range as JSON.
     *
     * Returns `[{"ph":0.0,"logd":2.5}, ...]` with `steps` evenly-spaced pH points.
     */
    logd_profile_json(ph_start: number, ph_end: number, steps: number): string;
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
     * Most acidic pKa in the molecule, or NaN if no acidic site.
     */
    pka_acid_value(): number;
    /**
     * Most basic pKa in the molecule, or NaN if no basic site.
     */
    pka_base_value(): number;
    /**
     * Quantitative Estimate of Drug-likeness (QED); range [0, 1].
     */
    qed(): number;
    /**
     * Randić connectivity index (χ₀).
     *
     * χ₀ = Σ 1/√(d_i × d_j) over all bonds, where d is heavy-atom degree.
     */
    randic_index(): number;
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
     * InChI string representation of the molecule.
     */
    to_inchi(): string;
    /**
     * InChIKey (27-character identifier) for the molecule.
     */
    to_inchikey(): string;
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
    /**
     * Zagreb index M1: Σ d_i² over all heavy atoms.
     */
    zagreb_index_m1(): number;
}

/**
 * Return a copy of the molecule with all implicit hydrogens converted to explicit H atoms.
 */
export function add_hydrogens(mol: MolHandle): MolHandle;

/**
 * Compute a full ADMET property profile for a molecule.
 *
 * Returns a JSON object with fields:
 * `bbb_score`, `bbb_passes`, `caco2`, `herg_risk`, `cyp3a4_risk`,
 * `pka_acid` (null if absent), `pka_base` (null if absent),
 * `esol`, `logd74`, `mw`, `logp`, `tpsa`, `hbd`, `hba`, `rotatable_bonds`
 *
 * Returns `{"error":"..."}` on parse failure.
 */
export function admet_profile_json(smiles: string): string;

/**
 * AtomPair fingerprint as a bit-packed byte vector (256 bytes = 2048 bits).
 */
export function atom_pair_bitvec(mol: MolHandle): Uint8Array;

/**
 * AutoCorr2D descriptor (7 values: topological distance lags 1-7).
 */
export function autocorr_2d_json(mol: MolHandle): string;

/**
 * AutoCorr3D descriptor (8 values: Euclidean distance bins 1-8 Å).
 * Requires 3D coordinates (generated automatically).
 */
export function autocorr_3d_json(mol: MolHandle): string;

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
 * Compute the canonical tautomer with specific atoms blocked from H-transfer.
 *
 * `blocked_atom_indices_json`: JSON array of 0-based atom indices, e.g. `[0, 3]`.
 * Any tautomer move whose donor, bridge, or acceptor is in the blocked set is suppressed.
 *
 * Returns canonical SMILES of the result, or `{"error":"..."}` on failure.
 * Out-of-range indices are silently ignored (no effect).
 */
export function canonical_tautomer_with_blocked_atoms_json(mol: MolHandle, blocked_atom_indices_json: string): string;

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
 * Compare multiple SMILES strings (up to 256 by default).
 * Accepts a delimiter-separated list (e.g., newline or comma).
 *
 * # Example (JS)
 * ```javascript
 * const smilesList = "c1ccccc1\nCc1ccccc1\nCCc1ccccc1";
 * const json = module.compare_molecules_batch_json(smilesList, "\n");
 * const comparison = JSON.parse(json);
 * ```
 */
export function compare_molecules_batch_json(smiles_batch: string, delimiter: string): string;

/**
 * Compare two or more SMILES strings (JSON string output).
 * Returns the JSON representation of a `MoleculeComparison` struct.
 *
 * # Example (JS)
 * ```javascript
 * const json = module.compare_molecules_json("c1ccccc1", "Cc1ccccc1");
 * const comparison = JSON.parse(json);
 * console.log(comparison.pairwise[0].similarities.ecfp4_tanimoto);
 * ```
 */
export function compare_molecules_json(smiles1: string, smiles2: string): string;

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
 * Infer bond connectivity and bond orders from an XYZ-format string.
 *
 * Explicit hydrogen atoms must be present in the XYZ for reliable bond-order
 * assignment (without H, carbonyl C=O cannot be distinguished from C-O).
 *
 * Returns JSON on success: `{"smiles":"CCO","atom_count":3,"bond_count":2}`.
 * `atom_count` and `bond_count` refer to the heavy-atom skeleton (H removed).
 *
 * Returns JSON on error: `{"error":"molecule has 450 atoms; maximum is 300"}`.
 *
 * Safe: never freezes. All internal loops are O(n²). Capped at 300 atoms.
 */
export function determine_bonds_from_xyz_json(xyz_str: string): string;

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
 * Like `ecfp4_bitvec` but with explicit chirality control.
 *
 * When `use_chirality=true`, tetrahedral stereochemistry is included in the
 * initial atom hash, making enantiomers have different fingerprints.
 * When `false` (default), chirality is ignored.
 */
export function ecfp4_bitvec_with_chirality(mol: MolHandle, use_chirality: boolean): Uint8Array;

/**
 * ECFP6 (radius-3) fingerprint as a bit-packed byte vector (256 bytes = 2048 bits).
 */
export function ecfp6_bitvec(mol: MolHandle): Uint8Array;

/**
 * Like `ecfp6_bitvec` but with explicit chirality control.
 *
 * When `use_chirality=true`, tetrahedral stereochemistry is included in the
 * initial atom hash, making enantiomers have different fingerprints.
 * When `false` (default), chirality is ignored.
 */
export function ecfp6_bitvec_with_chirality(mol: MolHandle, use_chirality: boolean): Uint8Array;

/**
 * Compute a fingerprint bit-vector with configurable ECFP radius and bit width.
 *
 * `radius` — Morgan radius (1 = ECFP2, 2 = ECFP4, 3 = ECFP6).
 * `nbits` — bit width; must be one of 256, 512, 1024, or 2048.
 *   Returns a `Uint8Array` of `nbits/8` bytes.
 *
 * The hash modulo is applied at fingerprint-generation time (`id % nbits`),
 * so no post-processing fold is needed.
 * Compute a custom ECFP (Extended Connectivity FingerPrint) with specified radius and bit count.
 *
 * When `use_chirality=true`, tetrahedral stereochemistry is included in the initial
 * atom hash. When `false` (default), chirality is ignored.
 */
export function ecfp_bitvec_custom(mol: MolHandle, radius: number, nbits: number, use_chirality: boolean): Uint8Array;

/**
 * Enumerate a combinatorial library from a SMIRKS template and two fragment sets.
 *
 * Generates all products by combining every scaffold with every building block.
 * Input format: `scaffolds_smiles` and `building_blocks_smiles` are pipe-delimited
 * SMILES strings (e.g., `"c1ccccc1|Cc1ccccc1"`).
 *
 * Returns JSON array of product SMILES strings.
 * Example: `enumerate_library_2way("[C:1][Cl].[C:2][NH2]>>[C:1]N[C:2]", "c1ccccc1|Cc1ccccc1", "NCc1ccccc1|NCC")`
 */
export function enumerate_library_2way(template: string, scaffolds_smiles: string, building_blocks_smiles: string): string;

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
 * Compute ERG-style 315-element float histogram fingerprint.
 * Returns JSON: {"len":315,"values":[f64,...]} or {"error":"..."}.
 * Format: 21 pharmacophore-feature-pair × 15 distance bins with Gaussian fuzzing.
 * See `chematic_fp::erg_vec` for details.
 */
export function erg_vec_json(mol: MolHandle): string;

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
 * Generate 3D coordinates using ETKDG and minimize with DREIDING force field.
 */
export function generate_3d_etkdg_minimized_pdb(mol: MolHandle): string;

/**
 * Generate 3D coordinates using ETKDG (torsion angle preferences) and return PDB block.
 * ETKDG produces higher-quality conformations than rule-based DG by applying
 * experimental torsion angle preferences to common structural patterns.
 */
export function generate_3d_etkdg_pdb(mol: MolHandle): string;

/**
 * Generate 3D coordinates from SMILES (raw distance geometry, no minimization).
 * Returns PDB format string with atoms positioned in 3D space.
 *
 * # Example (JS)
 * ```javascript
 * const pdbStr = module.generate_3d_from_smiles("c1ccccc1");
 * console.log(pdbStr);  // PDB file content
 * ```
 */
export function generate_3d_from_smiles(smiles: string): string;

/**
 * Generate energy-minimized 3D coordinates and return a PDB string.
 *
 * Runs distance-geometry placement followed by gradient-descent force-field
 * minimization.  Geometry quality is better than `generate_3d_pdb` for
 * flexible molecules; the force field is approximate (not MMFF94/UFF).
 */
export function generate_3d_minimized_pdb(mol: MolHandle): string;

/**
 * Generate 3D coordinates and minimize from SMILES string.
 * Pipeline: distance geometry → DREIDING minimization.
 * Better geometry quality than raw DG; suitable for graphics.
 *
 * # Example (JS)
 * ```javascript
 * const pdbStr = module.generate_3d_optimized_pdb("c1ccccc1");
 * console.log(pdbStr);  // PDB file with optimized geometry
 * ```
 */
export function generate_3d_optimized_pdb(smiles: string): string;

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
 * Get bond length in Ångströms between two atoms from a SMILES string.
 * Returns -1.0 if parsing fails or atom indices are out of range.
 *
 * # Arguments
 * - `smiles`: SMILES string
 * - `a`: first atom index
 * - `b`: second atom index
 *
 * # Example
 * ```javascript
 * const len = get_bond_length_json("CC", 0, 1);  // C-C single bond ≈ 1.54 Å
 * ```
 */
export function get_bond_length_json(smiles: string, a: number, b: number): number;

/**
 * All scalar molecular descriptors as a single JSON object.
 *
 * Keys use camelCase and match the individual `MolHandle` method names.
 * Drug-likeness rule outcomes are included as boolean fields.
 */
export function get_descriptors_json(mol: MolHandle): string;

/**
 * Get dihedral angle A—B—C—D in degrees from a SMILES string.
 * Returns null (JSON null) if any atom index is out of range or atoms are collinear.
 *
 * # Arguments
 * - `smiles`: SMILES string
 * - `a`, `b`, `c`, `d`: atom indices
 *
 * # Example
 * ```javascript
 * const dihedral = get_dihedral_json("CCCC", 0, 1, 2, 3);  // A-B-C-D
 * ```
 */
export function get_dihedral_json(smiles: string, a: number, b: number, c: number, d: number): any;

/**
 * Compute GETAWAY descriptors (GEometric, Topologic And wAveleT descriptors) from 3D coordinates.
 * Returns JSON array of 9 values: [G1, G2, G3, D1, D2, D3, T, V, A]
 * where G* = geometric autocorrelations (lag-1,2,3), D* = topologic distances,
 * T = total pairwise distance, V = bounding-box volume, A = anisotropy ratio.
 */
export function getaway_descriptors_json(mol: MolHandle): string;

/**
 * Identify functional groups. Returns a JSON array of objects:
 * `[{"atoms":[0,2,3],"type":"C,N,O"}, …]`
 */
export function identify_functional_groups(mol: MolHandle): string;

/**
 * Generate InChI string from SMILES.
 *
 * Returns `"error:<msg>"` on parse failure.
 */
export function inchi_from_smiles(smiles: string): string;

/**
 * Generate InChIKey from SMILES (27-character identifier).
 *
 * Returns `"error:<msg>"` on parse failure.
 */
export function inchikey_from_smiles(smiles: string): string;

/**
 * Invert the stereochemistry of a tetrahedral stereocenter (U/D wedge bonds).
 *
 * If the atom has no wedge/dash bonds, returns an unchanged copy.
 * Returns error if atom_idx is invalid.
 */
export function invert_stereocenter_at(mol: MolHandle, atom_idx: number): MolHandle;

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
 * MCS with ring-awareness constraints.
 *
 * `smiles_json` — JSON array of at least 2 SMILES strings.
 * `ring_matches_ring_only` — ring atoms may only match ring atoms.
 * `complete_rings_only` — partial ring inclusion is removed from the result.
 * Returns the MCS SMILES, or `"null"` when no common substructure was found.
 */
export function mcs_smiles_json_with_ring_config(smiles_json: string, ring_matches_ring_only: boolean, complete_rings_only: boolean): string;

/**
 * MinHash fingerprint (128 hashes) as JSON.
 *
 * Returns `{"num_hashes":128,"hashes":[u64,...]}`.
 * Use `tanimoto_mhfp_smiles` for direct SMILES-to-SMILES similarity.
 */
export function mhfp_hashes_json(mol: MolHandle): string;

/**
 * Optimize molecular geometry using DREIDING force field.
 *
 * Performs geometry minimization with DREIDING force field parameters.
 * Returns minimized coordinate PDB.
 *
 * # Arguments
 * * `mol` - Molecule to optimize
 *
 * # Returns
 * PDB format string with optimized coordinates
 */
export function minimize_dreiding_json(mol: MolHandle): string;

/**
 * Minimize geometry using MMFF94 steepest descent (Halgren 1996 full parameters).
 * Generates 3D coords internally if needed.
 * Returns JSON: {"energy":E,"rmsd":R,"converged":true,"iterations":N} or {"error":"..."}.
 */
export function minimize_mmff94_json(mol: MolHandle, max_iter: number): string;

/**
 * Minimize geometry using MMFF94 L-BFGS (faster convergence than steepest descent).
 * Returns JSON: {"energy":E,"rmsd":R,"converged":true,"iterations":N} or {"error":"..."}.
 */
export function minimize_mmff94_lbfgs_json(mol: MolHandle, max_iter: number): string;

/**
 * MMFF94 partial charges (BCI table, ±0.1e accuracy) as a JSON array of f64.
 *
 * Uses Bond Charge Increment (BCI) model (Halgren 1996) for 25 common bond types.
 * Returns `[q0, q1, ..., qN]` — one value per heavy atom.
 * Total charge equals the sum of formal charges (charge conserved).
 */
export function mmff94_charges_json(mol: MolHandle): string;

/**
 * Compute MMFF94-style atom-typed partial charges (improved over element-pair BCI).
 * Returns JSON: {"charges":[f64,...]} or {"error":"..."}.
 * Uses atom-type classification (Csp3/Ccarbonyl/Ohydroxyl/Oester/Nar/NarH etc.)
 * for better accuracy (~±0.02e) vs element-pair BCI (~±0.05e).
 */
export function mmff94_charges_typed_json(mol: MolHandle): string;

/**
 * Compute MMFF94 energy breakdown for current rule-based 3D geometry.
 * Returns JSON: {"bond":B,"angle":A,"torsion":T,"vdw":V,"elec":E,"total":X} or {"error":"..."}.
 */
export function mmff94_energy_breakdown_json(mol: MolHandle): string;

/**
 * Compute MMFF94 partial charges using numeric atom types (Halgren 1996 eq. 15).
 * Returns JSON: {"charges":[-0.28,0.15,...]} or {"error":"..."}.
 */
export function mmff94_partial_charges_json(mol: MolHandle): string;

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
 * Generate a complete molecular report (JSON string) from a SMILES.
 * Returns the JSON representation of a `MoleculeReport` struct.
 *
 * # Example (JS)
 * ```javascript
 * const json = module.molecule_report_json("CC(=O)Oc1ccccc1C(=O)O");
 * const report = JSON.parse(json);
 * console.log(report.canonical_smiles, report.descriptors.tpsa);
 * ```
 */
export function molecule_report_json(smiles: string): string;

/**
 * MQN descriptor (42 integer values: Molecular Quantum Numbers).
 */
export function mqn_json(mol: MolHandle): string;

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
 * Parse and re-serialize CXSMILES, preserving supported CX metadata.
 * Returns error if atom count exceeds 10,000.
 */
export function normalize_cxsmiles(s: string): string;

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
 * Parse CXSMARTS and return preserved metadata as JSON.
 * Returns error if atom count exceeds 10,000.
 */
export function parse_cxsmarts_json(s: string): string;

/**
 * Parse CXSMILES and return preserved metadata as JSON.
 *
 * Supported CX fields: atom labels (`$...$`), `atomProp`, atom radicals (`^n:`),
 * and zero-order bonds (`Z:`). The `cxsmiles` field is a re-serialized
 * round-trip form using the supported fields.
 * Returns error if atom count exceeds 10,000.
 */
export function parse_cxsmiles_json(s: string): string;

/**
 * Parse a SMILES string into a `MolHandle`.
 *
 * Returns a JS error string on parse failure or if atom count exceeds 10,000.
 */
export function parse_smiles(s: string): MolHandle;

/**
 * PEOE_VSA descriptors (14 bins) as a JSON array.
 */
export function peoe_vsa_json(mol: MolHandle): string;

/**
 * Detect pharmacophore features for virtual screening and lead optimization.
 * Returns JSON array of features: [{type, atom_idx, neighbor_count}, ...]
 */
export function pharmacophore_features_json(mol: MolHandle): string;

/**
 * Compute 2D pharmacophore fingerprint (2048 bits) as a JSON feature count summary.
 * Returns simplified JSON with feature type counts: {Donor, Acceptor, Aromatic, Hydrophobic, Positive, Negative}
 */
export function pharmacophore_fp_2d_summary(mol: MolHandle): string;

/**
 * Compute 3D pharmacophore fingerprint from generated 3D coordinates.
 * Returns simplified JSON with feature type counts (3D-aware version).
 */
export function pharmacophore_fp_3d_summary(mol: MolHandle): string;

/**
 * Predict pKa for all ionizable sites in a molecule.
 *
 * Returns a JSON array: `[{"atom_idx":8,"pka":4.0,"type":"acid","group":"carboxylic_acid"},...]`
 *
 * Returns `[]` if no ionizable sites are found, or `{"error":"..."}` on parse failure.
 */
export function predict_pka_json(smiles: string): string;

/**
 * Generate `count` random SMILES from a SMILES string using the given seed.
 * Atoms are permuted based on xorshift64 RNG. Each variant should parse back
 * to the same molecule. Returns a JSON array of SMILES strings.
 *
 * # Arguments
 * - `smiles`: input SMILES string
 * - `count`: number of variants to generate (capped at 100)
 * - `seed`: xorshift64 seed
 *
 * # Example
 * ```javascript
 * const variants = random_smiles_json("CC(C)O", 5, 42);
 * // variants: ["CC(C)O", "C(C)(O)C", ...]
 * ```
 */
export function random_smiles_json(smiles: string, count: number, seed: bigint): string;

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
 * Ring family classification and detection as JSON.
 * Returns an array of ring families with their atoms, ring indices, and topology kind.
 */
export function ring_families_json(mol: MolHandle): string;

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
 * Screen a batch of SMILES strings (JSON string output).
 * Returns per-record results including pass/fail with error details.
 * Includes MaxMin diversity picking and Butina clustering by default.
 *
 * # Example (JS)
 * ```javascript
 * const smilesList = "c1ccccc1\nCC\nCCC";
 * const json = module.screen_smiles_json(smilesList, "\n");
 * const report = JSON.parse(json);
 * console.log(report.records); // Array of ScreeningRecord
 * console.log(report.maxmin_picks); // Diversity-selected indices
 * console.log(report.butina_clusters); // Clustering result
 * ```
 */
export function screen_smiles_json(smiles_batch: string, delimiter: string): string;

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
 * Set dihedral angle A—B—C—D and return PDB block with modified coordinates.
 * Rotates the D-side subtree around the B—C bond.
 * Returns a JS error if parsing fails or atom indices are invalid.
 *
 * # Arguments
 * - `smiles`: SMILES string
 * - `a`, `b`, `c`, `d`: atom indices
 * - `angle_deg`: target dihedral angle in degrees
 *
 * # Example
 * ```javascript
 * const pdbBlock = set_dihedral_json("CCCC", 0, 1, 2, 3, 120.0);
 * ```
 */
export function set_dihedral_json(smiles: string, a: number, b: number, c: number, d: number, angle_deg: number): string;

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
 * Like `smarts_match_atoms` but with explicit chirality matching control.
 *
 * When `use_chirality=true`, SMARTS chirality primitives `[@]` and `[@@]` are
 * matched against the target molecule's stereochemistry. When `false`, chirality
 * is ignored (RDKit default).
 */
export function smarts_match_atoms_with_chirality(smarts: string, mol: MolHandle, use_chirality: boolean): string;

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

/**
 * Standardize a SMILES string and return result SMILES plus an audit report as JSON.
 *
 * Boolean flags map directly to `StandardizeOptions`.
 * Returns `"error:<msg>"` on parse or serialization failure.
 */
export function standardize_smiles_report_json(smiles: string, largest_fragment_only: boolean, neutralize_charges: boolean, remove_explicit_h: boolean, canonical_tautomer: boolean): string;

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
 * Tanimoto-like similarity between two SMILES via MHFP (MinHash Jaccard approximation).
 */
export function tanimoto_mhfp_smiles(smi1: string, smi2: string): number;

/**
 * Compute ECFP4 Tanimoto similarity from one query SMILES to all db SMILES (dense output).
 *
 * `db_smiles_json`: JSON array of SMILES strings (max 1024 via WASM_MAX_BATCH_ITEMS).
 *
 * Returns a flat JSON array of f32 scores, one per db entry, e.g. `[0.12,0.0,0.85]`.
 * No zero-filtering: the length always equals the number of db entries.
 * Returns `"error:<msg>"` on parse failure or oversized input.
 */
export function tanimoto_row_json(query_smi: string, db_smiles_json: string): string;

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
 * Scan a torsion dihedral i-j-k-l from 0° to 360° in `steps` increments.
 * Returns JSON array: [{"angle":0.0,"energy":E},...] or {"error":"..."}.
 */
export function torsion_scan_json(mol: MolHandle, i: number, j: number, k: number, l: number, steps: number): string;

/**
 * Virtual screen a query SMILES against a database of SMILES using ECFP4 Tanimoto.
 *
 * `db_smiles_json`: JSON array of SMILES strings (max 1024 via WASM_MAX_BATCH_ITEMS).
 * `k`: number of top hits to return; clamped to db size if larger.
 *
 * Returns JSON: `{"results":[{"rank":1,"score":0.85,"smiles":"CCO","idx":42},...]}`.
 * Returns `"error:<msg>"` on any parse failure or oversized input.
 */
export function virtual_screen_ecfp4_json(query_smi: string, db_smiles_json: string, k: number): string;

/**
 * Compute WHIM descriptors (Weighted Holistic Invariant Molecular) from 3D coordinates.
 * Returns JSON array of 10 values: [L1, L2, L3, P1, P2, P3, ALPHA, BETA, GAMMA, DELTA]
 * where L* = inertia tensor eigenvalues, P* = principal moments, ALPHA = sum of moments,
 * BETA = average pairwise interaction, GAMMA = geometric mean, DELTA = anisotropy.
 */
export function whim_descriptors_json(mol: MolHandle): string;

/**
 * Compute combined WHIM + GETAWAY descriptors (19 values total) as JSON array.
 * Useful for ML pipelines requiring both shape and topologic features.
 */
export function whim_getaway_combined_json(mol: MolHandle): string;

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
    readonly __wbg_mhfplshhandle_free: (a: number, b: number) => void;
    readonly __wbg_molhandle_free: (a: number, b: number) => void;
    readonly add_hydrogens: (a: number) => number;
    readonly admet_profile_json: (a: number, b: number) => [number, number];
    readonly atom_pair_bitvec: (a: number) => [number, number];
    readonly autocorr_2d_json: (a: number) => [number, number];
    readonly autocorr_3d_json: (a: number) => [number, number];
    readonly balance_check_json: (a: number, b: number) => [number, number];
    readonly brics_fragment_count: (a: number) => number;
    readonly brics_fragments_json: (a: number) => [number, number];
    readonly butina_cluster_ecfp4_json: (a: number, b: number, c: number) => [number, number, number, number];
    readonly canonical_tautomer: (a: number) => number;
    readonly canonical_tautomer_with_blocked_atoms_json: (a: number, b: number, c: number) => [number, number];
    readonly cdxml_to_smiles_json: (a: number, b: number) => [number, number, number, number];
    readonly cip_assignments_json: (a: number) => [number, number];
    readonly conformerhandle_add_generated_conformer: (a: number) => number;
    readonly conformerhandle_add_minimized_conformer: (a: number) => number;
    readonly conformerhandle_cluster_conformers_json: (a: number, b: number) => [number, number];
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
    readonly determine_bonds_from_xyz_json: (a: number, b: number) => [number, number];
    readonly dice_ecfp4: (a: number, b: number) => number;
    readonly dice_ecfp6: (a: number, b: number) => number;
    readonly dice_maccs: (a: number, b: number) => number;
    readonly ecfp4_bitvec: (a: number) => [number, number];
    readonly ecfp4_bitvec_with_chirality: (a: number, b: number) => [number, number];
    readonly ecfp6_bitvec: (a: number) => [number, number];
    readonly ecfp6_bitvec_with_chirality: (a: number, b: number) => [number, number];
    readonly ecfp_bitvec_custom: (a: number, b: number, c: number, d: number) => [number, number];
    readonly enumerate_library_2way: (a: number, b: number, c: number, d: number, e: number, f: number) => [number, number, number, number];
    readonly enumerate_stereo_isomers_json: (a: number) => [number, number, number, number];
    readonly enumerate_tautomers_json: (a: number) => [number, number];
    readonly erg_vec_json: (a: number) => [number, number];
    readonly estate_indices_json: (a: number) => [number, number];
    readonly fcfp4_bitvec: (a: number) => [number, number];
    readonly fcfp6_bitvec: (a: number) => [number, number];
    readonly find_reaction_center_json: (a: number, b: number) => [number, number];
    readonly gasteiger_charges_json: (a: number) => [number, number];
    readonly generate_3d_etkdg_minimized_pdb: (a: number) => [number, number];
    readonly generate_3d_etkdg_pdb: (a: number) => [number, number];
    readonly generate_3d_minimized_pdb: (a: number) => [number, number];
    readonly generate_3d_pdb: (a: number) => [number, number];
    readonly generic_murcko_scaffold: (a: number) => number;
    readonly get_atom_info: (a: number, b: number) => [number, number];
    readonly get_bond_between: (a: number, b: number, c: number) => [number, number];
    readonly get_bond_info: (a: number, b: number) => [number, number];
    readonly get_bond_length_json: (a: number, b: number, c: number, d: number) => number;
    readonly get_descriptors_json: (a: number) => [number, number];
    readonly get_dihedral_json: (a: number, b: number, c: number, d: number, e: number, f: number) => any;
    readonly getaway_descriptors_json: (a: number) => [number, number];
    readonly identify_functional_groups: (a: number) => [number, number];
    readonly inchi_from_smiles: (a: number, b: number) => [number, number];
    readonly inchikey_from_smiles: (a: number, b: number) => [number, number];
    readonly invert_stereocenter_at: (a: number, b: number) => [number, number, number];
    readonly is_valid_smiles: (a: number, b: number) => number;
    readonly labute_asa_per_atom_json: (a: number) => [number, number];
    readonly largest_fragment: (a: number) => number;
    readonly logp_per_atom_json: (a: number) => [number, number];
    readonly maccs_bitvec: (a: number) => [number, number];
    readonly match_smarts_smiles: (a: number, b: number, c: number, d: number) => [number, number, number, number];
    readonly maxmin_picks_ecfp4_json: (a: number, b: number, c: number) => [number, number, number, number];
    readonly mcs_smiles_json: (a: number, b: number) => [number, number, number, number];
    readonly mcs_smiles_json_with_ring_config: (a: number, b: number, c: number, d: number) => [number, number, number, number];
    readonly mhfp_hashes_json: (a: number) => [number, number];
    readonly mhfplshhandle_add_smiles: (a: number, b: number, c: number) => [number, number, number];
    readonly mhfplshhandle_is_empty: (a: number) => number;
    readonly mhfplshhandle_len: (a: number) => number;
    readonly mhfplshhandle_new: (a: number) => number;
    readonly mhfplshhandle_query_json: (a: number, b: number, c: number, d: number) => [number, number, number, number];
    readonly minimize_dreiding_json: (a: number) => [number, number];
    readonly minimize_mmff94_json: (a: number, b: number) => [number, number];
    readonly minimize_mmff94_lbfgs_json: (a: number, b: number) => [number, number];
    readonly mmff94_charges_json: (a: number) => [number, number];
    readonly mmff94_charges_typed_json: (a: number) => [number, number];
    readonly mmff94_energy_breakdown_json: (a: number) => [number, number];
    readonly mmff94_partial_charges_json: (a: number) => [number, number];
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
    readonly molhandle_assign_cip_json: (a: number) => [number, number];
    readonly molhandle_bbb_passes: (a: number) => number;
    readonly molhandle_bbb_score: (a: number) => number;
    readonly molhandle_bertz_ct: (a: number) => number;
    readonly molhandle_bond_count: (a: number) => number;
    readonly molhandle_caco2_permeability: (a: number) => number;
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
    readonly molhandle_cyp3a4_inhibition_risk: (a: number) => number;
    readonly molhandle_depict_png: (a: number) => [number, number];
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
    readonly molhandle_herg_risk_score: (a: number) => number;
    readonly molhandle_isotope_distribution_json: (a: number, b: number) => [number, number];
    readonly molhandle_iupac_name: (a: number) => [number, number];
    readonly molhandle_kappa1: (a: number) => number;
    readonly molhandle_kappa2: (a: number) => number;
    readonly molhandle_kappa3: (a: number) => number;
    readonly molhandle_labute_asa: (a: number) => number;
    readonly molhandle_lipinski_passes: (a: number) => number;
    readonly molhandle_logd_at_ph: (a: number, b: number) => number;
    readonly molhandle_logd_profile_json: (a: number, b: number, c: number, d: number) => [number, number];
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
    readonly molhandle_pka_acid_value: (a: number) => number;
    readonly molhandle_pka_base_value: (a: number) => number;
    readonly molhandle_qed: (a: number) => number;
    readonly molhandle_randic_index: (a: number) => number;
    readonly molhandle_reos_passes: (a: number) => number;
    readonly molhandle_ring_count: (a: number) => number;
    readonly molhandle_rotatable_bond_count: (a: number) => number;
    readonly molhandle_sum_estate: (a: number) => number;
    readonly molhandle_to_inchi: (a: number) => [number, number];
    readonly molhandle_to_inchikey: (a: number) => [number, number];
    readonly molhandle_tpsa: (a: number) => number;
    readonly molhandle_veber_passes: (a: number) => number;
    readonly molhandle_wiener_index: (a: number) => number;
    readonly molhandle_zagreb_index_m1: (a: number) => number;
    readonly mqn_json: (a: number) => [number, number];
    readonly mr_per_atom_json: (a: number) => [number, number];
    readonly murcko_scaffold: (a: number) => number;
    readonly nearest_neighbors_json: (a: number, b: number, c: number, d: number, e: number) => [number, number];
    readonly neutralize_charges: (a: number) => number;
    readonly normalize_cxsmiles: (a: number, b: number) => [number, number, number, number];
    readonly normalize_reaction_smiles: (a: number, b: number) => [number, number, number, number];
    readonly pains_matches_json: (a: number) => [number, number];
    readonly parse_cxsmarts_json: (a: number, b: number) => [number, number, number, number];
    readonly parse_cxsmiles_json: (a: number, b: number) => [number, number, number, number];
    readonly parse_smiles: (a: number, b: number) => [number, number, number];
    readonly peoe_vsa_json: (a: number) => [number, number];
    readonly pharmacophore_features_json: (a: number) => [number, number];
    readonly pharmacophore_fp_2d_summary: (a: number) => [number, number];
    readonly pharmacophore_fp_3d_summary: (a: number) => [number, number];
    readonly predict_pka_json: (a: number, b: number) => [number, number];
    readonly random_smiles_json: (a: number, b: number, c: number, d: bigint) => [number, number, number, number];
    readonly remove_hydrogens: (a: number) => number;
    readonly rgroup_decompose_json: (a: number, b: number, c: number, d: number) => [number, number, number, number];
    readonly ring_families_json: (a: number) => [number, number, number, number];
    readonly run_md_json: (a: number, b: number, c: number) => [number, number];
    readonly run_reactants: (a: number, b: number, c: number, d: number) => [number, number, number, number];
    readonly sa_score: (a: number) => number;
    readonly sdf_from_records_json: (a: number, b: number, c: number, d: number, e: number, f: number) => [number, number, number, number];
    readonly sdf_to_records_json: (a: number, b: number) => [number, number];
    readonly sdf_to_smiles_json: (a: number, b: number) => [number, number];
    readonly set_dihedral_json: (a: number, b: number, c: number, d: number, e: number, f: number, g: number) => [number, number, number, number];
    readonly shape_descriptors_json: (a: number) => [number, number];
    readonly slogp_vsa_json: (a: number) => [number, number];
    readonly smarts_match_atoms: (a: number, b: number, c: number) => [number, number, number, number];
    readonly smarts_match_atoms_with_chirality: (a: number, b: number, c: number, d: number) => [number, number, number, number];
    readonly smiles_array_to_sdf: (a: number, b: number) => [number, number, number, number];
    readonly smiles_to_mol2: (a: number, b: number) => [number, number];
    readonly smiles_to_svg_highlighted: (a: number, b: number, c: number, d: number, e: number, f: number, g: number, h: number) => [number, number, number, number];
    readonly smr_vsa_json: (a: number) => [number, number];
    readonly sssr_rings_json: (a: number) => [number, number];
    readonly standardize_smiles: (a: number, b: number) => [number, number];
    readonly standardize_smiles_report_json: (a: number, b: number, c: number, d: number, e: number, f: number) => [number, number];
    readonly tanimoto_atom_pair: (a: number, b: number) => number;
    readonly tanimoto_ecfp4: (a: number, b: number) => number;
    readonly tanimoto_ecfp6: (a: number, b: number) => number;
    readonly tanimoto_fcfp4: (a: number, b: number) => number;
    readonly tanimoto_fcfp6: (a: number, b: number) => number;
    readonly tanimoto_maccs: (a: number, b: number) => number;
    readonly tanimoto_mhfp_smiles: (a: number, b: number, c: number, d: number) => [number, number, number];
    readonly tanimoto_row_json: (a: number, b: number, c: number, d: number) => [number, number];
    readonly tanimoto_smiles: (a: number, b: number, c: number, d: number) => [number, number, number];
    readonly tanimoto_topo_path: (a: number, b: number) => number;
    readonly tanimoto_torsion: (a: number, b: number) => number;
    readonly to_cml: (a: number) => [number, number];
    readonly to_mol_block: (a: number) => [number, number];
    readonly to_mol_v3000_block: (a: number) => [number, number];
    readonly to_xyz: (a: number) => [number, number];
    readonly torsion_bitvec: (a: number) => [number, number];
    readonly torsion_scan_json: (a: number, b: number, c: number, d: number, e: number, f: number) => [number, number];
    readonly virtual_screen_ecfp4_json: (a: number, b: number, c: number, d: number, e: number) => [number, number];
    readonly whim_descriptors_json: (a: number) => [number, number];
    readonly whim_getaway_combined_json: (a: number) => [number, number];
    readonly write_smiles: (a: number) => [number, number];
    readonly start: () => void;
    readonly molhandle_atom_count: (a: number) => number;
    readonly compare_molecules_batch_json: (a: number, b: number, c: number, d: number) => [number, number, number, number];
    readonly compare_molecules_json: (a: number, b: number, c: number, d: number) => [number, number, number, number];
    readonly generate_3d_from_smiles: (a: number, b: number) => [number, number, number, number];
    readonly generate_3d_optimized_pdb: (a: number, b: number) => [number, number, number, number];
    readonly molecule_report_json: (a: number, b: number) => [number, number, number, number];
    readonly screen_smiles_json: (a: number, b: number, c: number, d: number) => [number, number];
    readonly __wbindgen_malloc: (a: number, b: number) => number;
    readonly __wbindgen_realloc: (a: number, b: number, c: number, d: number) => number;
    readonly __wbindgen_free: (a: number, b: number, c: number) => void;
    readonly __wbindgen_exn_store: (a: number) => void;
    readonly __externref_table_alloc: () => number;
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
