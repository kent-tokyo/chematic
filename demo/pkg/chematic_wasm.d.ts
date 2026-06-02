/* tslint:disable */
/* eslint-disable */

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
     * Number of spiro atoms (sole shared atom between exactly 2 rings).
     */
    num_spiro_atoms(): number;
    /**
     * Number of assigned stereocenters (R/S).
     */
    num_stereocenters(): number;
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
 * Number of BRICS fragments produced by fragmenting the molecule.
 *
 * Returns 1 if no BRICS-breakable bonds exist (whole molecule is one fragment).
 */
export function brics_fragment_count(mol: MolHandle): number;

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
 * Detect named functional groups in `mol`.
 *
 * Returns a JSON array of `{"name":"hydroxyl","atoms":[3]}` objects.
 * Multiple matches of the same group (e.g. two hydroxyl groups) each appear
 * as a separate entry.  Overlapping groups (carboxylic acid → "carboxyl" +
 * "hydroxyl" + "carbonyl") are all returned.
 */
export function detect_functional_groups(mol: MolHandle): string;

/**
 * Compute the ECFP4 fingerprint as a bit-packed byte vector (256 bytes = 2048 bits).
 */
export function ecfp4_bitvec(mol: MolHandle): Uint8Array;

/**
 * Per-atom EState values as a JSON array of f64.
 *
 * Indices match `mol.atoms()` order.  Hydrogen atoms get 0.0.
 */
export function estate_indices_json(mol: MolHandle): string;

/**
 * Gasteiger-Marsili PEOE partial charges as a JSON array of f64.
 */
export function gasteiger_charges_json(mol: MolHandle): string;

/**
 * Generate 3D coordinates for the molecule and return a PDB string.
 *
 * Coordinates are generated using distance-geometry placement with ring templates.
 * Returns heavy-atom PDB (HETATM records, no explicit H).
 */
export function generate_3d_pdb(mol: MolHandle): string;

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
 * Identify functional groups. Returns a JSON array of objects:
 * `[{"atoms":[0,2,3],"type":"C,N,O"}, …]`
 */
export function identify_functional_groups(mol: MolHandle): string;

/**
 * Returns `true` if the SMILES string can be parsed without error.
 */
export function is_valid_smiles(s: string): boolean;

/**
 * Find all SMARTS matches in a molecule given only SMILES strings.
 *
 * Convenience wrapper around `smarts_match_atoms` that accepts raw SMILES
 * instead of a `MolHandle`.  Returns the same JSON format: `[[0,1],[3,4]]`.
 * Returns a JS error on SMILES or SMARTS parse failure.
 */
export function match_smarts_smiles(smiles: string, smarts: string): string;

/**
 * Serialize a SMILES string directly to a MOL V2000 block.
 *
 * Convenience wrapper; all atom coordinates are 0.0.
 * Returns a JS error on SMILES parse failure.
 */
export function mol_block_from_smiles(smiles: string): string;

/**
 * Parse a MOL V2000 block and return a `MolHandle`.
 *
 * Returns a JS error string on parse failure.
 */
export function mol_from_sdf_block(block: string): MolHandle;

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
 * Parse an SDF string and return a JSON array of canonical SMILES strings.
 *
 * Invalid records are represented as `null` in the array.
 */
export function sdf_to_smiles_json(sdf: string): string;

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
 * Tanimoto similarity between two molecules using FCFP4 fingerprints (pharmacophore-based).
 */
export function tanimoto_fcfp4(a: MolHandle, b: MolHandle): number;

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
 * Serialize a molecule to a MOL V2000 block.
 *
 * All atom coordinates are written as 0.0 (the `Molecule` type has no 2D
 * coordinate storage; real coordinates would require a separate layout pass).
 */
export function to_mol_block(mol: MolHandle): string;

export type InitInput = RequestInfo | URL | Response | BufferSource | WebAssembly.Module;

export interface InitOutput {
    readonly memory: WebAssembly.Memory;
    readonly __wbg_depictoptions_free: (a: number, b: number) => void;
    readonly __wbg_molhandle_free: (a: number, b: number) => void;
    readonly add_hydrogens: (a: number) => number;
    readonly brics_fragment_count: (a: number) => number;
    readonly depict_reaction_svg: (a: number, b: number) => [number, number, number, number];
    readonly depict_svg_grid: (a: number, b: number, c: number) => [number, number];
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
    readonly ecfp4_bitvec: (a: number) => [number, number];
    readonly estate_indices_json: (a: number) => [number, number];
    readonly gasteiger_charges_json: (a: number) => [number, number];
    readonly generate_3d_pdb: (a: number) => [number, number];
    readonly get_atom_info: (a: number, b: number) => [number, number];
    readonly get_bond_between: (a: number, b: number, c: number) => [number, number];
    readonly get_bond_info: (a: number, b: number) => [number, number];
    readonly identify_functional_groups: (a: number) => [number, number];
    readonly is_valid_smiles: (a: number, b: number) => number;
    readonly match_smarts_smiles: (a: number, b: number, c: number, d: number) => [number, number, number, number];
    readonly mol_block_from_smiles: (a: number, b: number) => [number, number, number, number];
    readonly mol_from_sdf_block: (a: number, b: number) => [number, number, number];
    readonly molhandle_aromatic_ring_count: (a: number) => number;
    readonly molhandle_atom_count: (a: number) => number;
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
    readonly molhandle_num_aromatic_heterocycles: (a: number) => number;
    readonly molhandle_num_bridgehead_atoms: (a: number) => number;
    readonly molhandle_num_heteroatoms: (a: number) => number;
    readonly molhandle_num_saturated_heterocycles: (a: number) => number;
    readonly molhandle_num_spiro_atoms: (a: number) => number;
    readonly molhandle_num_stereocenters: (a: number) => number;
    readonly molhandle_pains_passes: (a: number) => number;
    readonly molhandle_qed: (a: number) => number;
    readonly molhandle_reos_passes: (a: number) => number;
    readonly molhandle_ring_count: (a: number) => number;
    readonly molhandle_rotatable_bond_count: (a: number) => number;
    readonly molhandle_sum_estate: (a: number) => number;
    readonly molhandle_tpsa: (a: number) => number;
    readonly molhandle_veber_passes: (a: number) => number;
    readonly molhandle_wiener_index: (a: number) => number;
    readonly parse_smiles: (a: number, b: number) => [number, number, number];
    readonly peoe_vsa_json: (a: number) => [number, number];
    readonly remove_hydrogens: (a: number) => number;
    readonly run_reactants: (a: number, b: number, c: number, d: number) => [number, number, number, number];
    readonly sa_score: (a: number) => number;
    readonly sdf_to_smiles_json: (a: number, b: number) => [number, number];
    readonly slogp_vsa_json: (a: number) => [number, number];
    readonly smarts_match_atoms: (a: number, b: number, c: number) => [number, number, number, number];
    readonly smiles_to_svg_highlighted: (a: number, b: number, c: number, d: number, e: number, f: number, g: number, h: number) => [number, number, number, number];
    readonly smr_vsa_json: (a: number) => [number, number];
    readonly tanimoto_atom_pair: (a: number, b: number) => number;
    readonly tanimoto_ecfp4: (a: number, b: number) => number;
    readonly tanimoto_fcfp4: (a: number, b: number) => number;
    readonly tanimoto_smiles: (a: number, b: number, c: number, d: number) => [number, number, number];
    readonly tanimoto_topo_path: (a: number, b: number) => number;
    readonly tanimoto_torsion: (a: number, b: number) => number;
    readonly to_mol_block: (a: number) => [number, number];
    readonly start: () => void;
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
