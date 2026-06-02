/* @ts-self-types="./chematic_wasm.d.ts" */

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
    __destroy_into_raw() {
        const ptr = this.__wbg_ptr;
        this.__wbg_ptr = 0;
        DepictOptionsFinalization.unregister(this);
        return ptr;
    }
    free() {
        const ptr = this.__destroy_into_raw();
        wasm.__wbg_depictoptions_free(ptr, 0);
    }
    constructor() {
        const ret = wasm.depictoptions_new();
        this.__wbg_ptr = ret;
        DepictOptionsFinalization.register(this, this.__wbg_ptr, this);
        return this;
    }
    /**
     * Set a per-atom color override (CSS color string).  Calling multiple times
     * for the same `idx` uses the last value.  The atom is highlighted even if
     * not in `set_highlight_atoms`.
     * @param {number} idx
     * @param {string} color
     */
    set_atom_color(idx, color) {
        const ptr0 = passStringToWasm0(color, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        wasm.depictoptions_set_atom_color(this.__wbg_ptr, idx, ptr0, len0);
    }
    /**
     * @param {boolean} v
     */
    set_atom_ids(v) {
        wasm.depictoptions_set_atom_ids(this.__wbg_ptr, v);
    }
    /**
     * @param {string} bg
     */
    set_background(bg) {
        const ptr0 = passStringToWasm0(bg, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        wasm.depictoptions_set_background(this.__wbg_ptr, ptr0, len0);
    }
    /**
     * @param {boolean} dark
     */
    set_dark(dark) {
        wasm.depictoptions_set_dark(this.__wbg_ptr, dark);
    }
    /**
     * @param {number} h
     */
    set_height(h) {
        wasm.depictoptions_set_height(this.__wbg_ptr, h);
    }
    /**
     * @param {Uint32Array} atoms
     */
    set_highlight_atoms(atoms) {
        const ptr0 = passArray32ToWasm0(atoms, wasm.__wbindgen_malloc);
        const len0 = WASM_VECTOR_LEN;
        wasm.depictoptions_set_highlight_atoms(this.__wbg_ptr, ptr0, len0);
    }
    /**
     * @param {Uint32Array} bonds
     */
    set_highlight_bonds(bonds) {
        const ptr0 = passArray32ToWasm0(bonds, wasm.__wbindgen_malloc);
        const len0 = WASM_VECTOR_LEN;
        wasm.depictoptions_set_highlight_bonds(this.__wbg_ptr, ptr0, len0);
    }
    /**
     * @param {string} color
     */
    set_highlight_color(color) {
        const ptr0 = passStringToWasm0(color, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        wasm.depictoptions_set_highlight_color(this.__wbg_ptr, ptr0, len0);
    }
    /**
     * @param {boolean} v
     */
    set_kekulize(v) {
        wasm.depictoptions_set_kekulize(this.__wbg_ptr, v);
    }
    /**
     * @param {number} p
     */
    set_padding(p) {
        wasm.depictoptions_set_padding(this.__wbg_ptr, p);
    }
    /**
     * @param {boolean} v
     */
    set_show_atom_indices(v) {
        wasm.depictoptions_set_show_atom_indices(this.__wbg_ptr, v);
    }
    /**
     * @param {number} w
     */
    set_width(w) {
        wasm.depictoptions_set_width(this.__wbg_ptr, w);
    }
}
if (Symbol.dispose) DepictOptions.prototype[Symbol.dispose] = DepictOptions.prototype.free;

/**
 * A handle to a parsed molecule.  Owns the molecule behind an `Rc` so that
 * it can be cheaply cloned on the JS side without copying atom/bond data.
 */
export class MolHandle {
    static __wrap(ptr) {
        const obj = Object.create(MolHandle.prototype);
        obj.__wbg_ptr = ptr;
        MolHandleFinalization.register(obj, obj.__wbg_ptr, obj);
        return obj;
    }
    __destroy_into_raw() {
        const ptr = this.__wbg_ptr;
        this.__wbg_ptr = 0;
        MolHandleFinalization.unregister(this);
        return ptr;
    }
    free() {
        const ptr = this.__destroy_into_raw();
        wasm.__wbg_molhandle_free(ptr, 0);
    }
    /**
     * Number of aromatic rings (all ring atoms aromatic).
     * @returns {number}
     */
    aromatic_ring_count() {
        const ret = wasm.molhandle_aromatic_ring_count(this.__wbg_ptr);
        return ret >>> 0;
    }
    /**
     * Number of heavy atoms (explicit atoms in the graph; does not count implicit H).
     * @returns {number}
     */
    atom_count() {
        const ret = wasm.molhandle_atom_count(this.__wbg_ptr);
        return ret >>> 0;
    }
    /**
     * Bertz complexity index (BertzCT).
     * @returns {number}
     */
    bertz_ct() {
        const ret = wasm.molhandle_bertz_ct(this.__wbg_ptr);
        return ret;
    }
    /**
     * Number of bonds.
     * @returns {number}
     */
    bond_count() {
        const ret = wasm.molhandle_bond_count(this.__wbg_ptr);
        return ret >>> 0;
    }
    /**
     * Canonical SMILES string.
     * @returns {string}
     */
    canonical_smiles() {
        let deferred1_0;
        let deferred1_1;
        try {
            const ret = wasm.molhandle_canonical_smiles(this.__wbg_ptr);
            deferred1_0 = ret[0];
            deferred1_1 = ret[1];
            return getStringFromWasm0(ret[0], ret[1]);
        } finally {
            wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
        }
    }
    /**
     * Kier–Hall χ0 molecular connectivity index.
     * @returns {number}
     */
    chi0() {
        const ret = wasm.molhandle_chi0(this.__wbg_ptr);
        return ret;
    }
    /**
     * Kier–Hall χ0v valence-weighted connectivity index.
     * @returns {number}
     */
    chi0v() {
        const ret = wasm.molhandle_chi0v(this.__wbg_ptr);
        return ret;
    }
    /**
     * Kier–Hall χ1 molecular connectivity index.
     * @returns {number}
     */
    chi1() {
        const ret = wasm.molhandle_chi1(this.__wbg_ptr);
        return ret;
    }
    /**
     * Kier–Hall χ1v valence-weighted connectivity index.
     * @returns {number}
     */
    chi1v() {
        const ret = wasm.molhandle_chi1v(this.__wbg_ptr);
        return ret;
    }
    /**
     * Kier–Hall χ2 molecular connectivity index.
     * @returns {number}
     */
    chi2() {
        const ret = wasm.molhandle_chi2(this.__wbg_ptr);
        return ret;
    }
    /**
     * Kier–Hall χ2v valence-weighted connectivity index.
     * @returns {number}
     */
    chi2v() {
        const ret = wasm.molhandle_chi2v(this.__wbg_ptr);
        return ret;
    }
    /**
     * Kier–Hall χ3 molecular connectivity index.
     * @returns {number}
     */
    chi3() {
        const ret = wasm.molhandle_chi3(this.__wbg_ptr);
        return ret;
    }
    /**
     * Kier–Hall χ3v valence-weighted connectivity index.
     * @returns {number}
     */
    chi3v() {
        const ret = wasm.molhandle_chi3v(this.__wbg_ptr);
        return ret;
    }
    /**
     * Kier–Hall χ4 molecular connectivity index.
     * @returns {number}
     */
    chi4() {
        const ret = wasm.molhandle_chi4(this.__wbg_ptr);
        return ret;
    }
    /**
     * Kier–Hall χ4v valence-weighted connectivity index.
     * @returns {number}
     */
    chi4v() {
        const ret = wasm.molhandle_chi4v(this.__wbg_ptr);
        return ret;
    }
    /**
     * 2D SVG depiction of the molecule (CPK coloring).
     * @returns {string}
     */
    depict_svg() {
        let deferred1_0;
        let deferred1_1;
        try {
            const ret = wasm.molhandle_depict_svg(this.__wbg_ptr);
            deferred1_0 = ret[0];
            deferred1_1 = ret[1];
            return getStringFromWasm0(ret[0], ret[1]);
        } finally {
            wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
        }
    }
    /**
     * 2D SVG depiction with style options.
     * @param {DepictOptions} opts
     * @returns {string}
     */
    depict_svg_opts(opts) {
        let deferred1_0;
        let deferred1_1;
        try {
            _assertClass(opts, DepictOptions);
            const ret = wasm.molhandle_depict_svg_opts(this.__wbg_ptr, opts.__wbg_ptr);
            deferred1_0 = ret[0];
            deferred1_1 = ret[1];
            return getStringFromWasm0(ret[0], ret[1]);
        } finally {
            wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
        }
    }
    /**
     * Returns `true` if the molecule passes Egan's absorption criteria
     * (TPSA ≤ 131.6 Å² and LogP ≤ 5.88).
     * @returns {boolean}
     */
    egan_passes() {
        const ret = wasm.molhandle_egan_passes(this.__wbg_ptr);
        return ret !== 0;
    }
    /**
     * Monoisotopic (exact) mass.
     * @returns {number}
     */
    exact_mass() {
        const ret = wasm.molhandle_exact_mass(this.__wbg_ptr);
        return ret;
    }
    /**
     * Sum of formal charges.
     * @returns {number}
     */
    formal_charge_sum() {
        const ret = wasm.molhandle_formal_charge_sum(this.__wbg_ptr);
        return ret;
    }
    /**
     * Molecular formula string (Hill notation: C first, H second, then alphabetical).
     * @returns {string}
     */
    formula() {
        let deferred1_0;
        let deferred1_1;
        try {
            const ret = wasm.molhandle_formula(this.__wbg_ptr);
            deferred1_0 = ret[0];
            deferred1_1 = ret[1];
            return getStringFromWasm0(ret[0], ret[1]);
        } finally {
            wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
        }
    }
    /**
     * Fraction of sp3 carbons (Fsp3).
     * @returns {number}
     */
    fsp3() {
        const ret = wasm.molhandle_fsp3(this.__wbg_ptr);
        return ret;
    }
    /**
     * Returns `true` if the molecule passes Ghose's drug-likeness filter
     * (MW 160–480, LogP −0.4–5.6, HeavyAtoms 20–70, MR 40–130).
     * @returns {boolean}
     */
    ghose_passes() {
        const ret = wasm.molhandle_ghose_passes(this.__wbg_ptr);
        return ret !== 0;
    }
    /**
     * Number of hydrogen bond acceptors (Lipinski: all N and O atoms).
     * @returns {number}
     */
    hba_count() {
        const ret = wasm.molhandle_hba_count(this.__wbg_ptr);
        return ret >>> 0;
    }
    /**
     * Number of hydrogen bond donors (N-H or O-H groups).
     * @returns {number}
     */
    hbd_count() {
        const ret = wasm.molhandle_hbd_count(this.__wbg_ptr);
        return ret >>> 0;
    }
    /**
     * Number of non-hydrogen heavy atoms.
     * @returns {number}
     */
    heavy_atom_count() {
        const ret = wasm.molhandle_heavy_atom_count(this.__wbg_ptr);
        return ret >>> 0;
    }
    /**
     * Hall–Kier κ1 shape index.
     * @returns {number}
     */
    kappa1() {
        const ret = wasm.molhandle_kappa1(this.__wbg_ptr);
        return ret;
    }
    /**
     * Hall–Kier κ2 shape index.
     * @returns {number}
     */
    kappa2() {
        const ret = wasm.molhandle_kappa2(this.__wbg_ptr);
        return ret;
    }
    /**
     * Hall–Kier κ3 shape index.
     * @returns {number}
     */
    kappa3() {
        const ret = wasm.molhandle_kappa3(this.__wbg_ptr);
        return ret;
    }
    /**
     * Labute approximate surface area (Å²).
     * @returns {number}
     */
    labute_asa() {
        const ret = wasm.molhandle_labute_asa(this.__wbg_ptr);
        return ret;
    }
    /**
     * Returns `true` if the molecule satisfies Lipinski's Rule of Five.
     * @returns {boolean}
     */
    lipinski_passes() {
        const ret = wasm.molhandle_lipinski_passes(this.__wbg_ptr);
        return ret !== 0;
    }
    /**
     * Crippen–Wildman octanol/water partition coefficient (LogP).
     * @returns {number}
     */
    logp_crippen() {
        const ret = wasm.molhandle_logp_crippen(this.__wbg_ptr);
        return ret;
    }
    /**
     * Maximum EState index across all heavy atoms.
     * @returns {number}
     */
    max_estate() {
        const ret = wasm.molhandle_max_estate(this.__wbg_ptr);
        return ret;
    }
    /**
     * Minimum EState index across all heavy atoms.
     * @returns {number}
     */
    min_estate() {
        const ret = wasm.molhandle_min_estate(this.__wbg_ptr);
        return ret;
    }
    /**
     * Wildman–Crippen molar refractivity (MR).
     * @returns {number}
     */
    molar_refractivity() {
        const ret = wasm.molhandle_molar_refractivity(this.__wbg_ptr);
        return ret;
    }
    /**
     * Average molecular weight (Da).
     * @returns {number}
     */
    molecular_weight() {
        const ret = wasm.molhandle_molecular_weight(this.__wbg_ptr);
        return ret;
    }
    /**
     * Morgan count fingerprint as a JSON object string (`{"<hash>": count, …}`).
     *
     * `radius` controls the ECFP radius (2 = ECFP4-equivalent).
     * @param {number} radius
     * @returns {string}
     */
    morgan_fp_counts_json(radius) {
        let deferred1_0;
        let deferred1_1;
        try {
            const ret = wasm.molhandle_morgan_fp_counts_json(this.__wbg_ptr, radius);
            deferred1_0 = ret[0];
            deferred1_1 = ret[1];
            return getStringFromWasm0(ret[0], ret[1]);
        } finally {
            wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
        }
    }
    /**
     * Number of non-aromatic rings containing at least one heteroatom.
     * @returns {number}
     */
    num_aliphatic_heterocycles() {
        const ret = wasm.molhandle_num_aliphatic_heterocycles(this.__wbg_ptr);
        return ret >>> 0;
    }
    /**
     * Number of aromatic rings containing at least one heteroatom (N, O, S, …).
     * @returns {number}
     */
    num_aromatic_heterocycles() {
        const ret = wasm.molhandle_num_aromatic_heterocycles(this.__wbg_ptr);
        return ret >>> 0;
    }
    /**
     * Number of bridgehead atoms (shared by ≥2 rings with ≥3 ring bonds).
     * @returns {number}
     */
    num_bridgehead_atoms() {
        const ret = wasm.molhandle_num_bridgehead_atoms(this.__wbg_ptr);
        return ret >>> 0;
    }
    /**
     * Number of heteroatoms (non-C, non-H heavy atoms).
     * @returns {number}
     */
    num_heteroatoms() {
        const ret = wasm.molhandle_num_heteroatoms(this.__wbg_ptr);
        return ret >>> 0;
    }
    /**
     * Number of fully saturated rings containing at least one heteroatom.
     * @returns {number}
     */
    num_saturated_heterocycles() {
        const ret = wasm.molhandle_num_saturated_heterocycles(this.__wbg_ptr);
        return ret >>> 0;
    }
    /**
     * Number of spiro atoms (sole shared atom between exactly 2 rings).
     * @returns {number}
     */
    num_spiro_atoms() {
        const ret = wasm.molhandle_num_spiro_atoms(this.__wbg_ptr);
        return ret >>> 0;
    }
    /**
     * Number of assigned stereocenters (R/S).
     * @returns {number}
     */
    num_stereocenters() {
        const ret = wasm.molhandle_num_stereocenters(this.__wbg_ptr);
        return ret >>> 0;
    }
    /**
     * Returns `true` if the molecule has no PAINS structural alerts.
     * @returns {boolean}
     */
    pains_passes() {
        const ret = wasm.molhandle_pains_passes(this.__wbg_ptr);
        return ret !== 0;
    }
    /**
     * Quantitative Estimate of Drug-likeness (QED); range [0, 1].
     * @returns {number}
     */
    qed() {
        const ret = wasm.molhandle_qed(this.__wbg_ptr);
        return ret;
    }
    /**
     * Returns `true` if the molecule passes the REOS (Rapid Elimination Of Swill) filter.
     * @returns {boolean}
     */
    reos_passes() {
        const ret = wasm.molhandle_reos_passes(this.__wbg_ptr);
        return ret !== 0;
    }
    /**
     * Total number of rings (SSSR count).
     * @returns {number}
     */
    ring_count() {
        const ret = wasm.molhandle_ring_count(this.__wbg_ptr);
        return ret >>> 0;
    }
    /**
     * Number of rotatable bonds.
     * @returns {number}
     */
    rotatable_bond_count() {
        const ret = wasm.molhandle_rotatable_bond_count(this.__wbg_ptr);
        return ret >>> 0;
    }
    /**
     * Sum of EState indices over all heavy atoms.
     * @returns {number}
     */
    sum_estate() {
        const ret = wasm.molhandle_sum_estate(this.__wbg_ptr);
        return ret;
    }
    /**
     * Topological polar surface area (Å²).
     * @returns {number}
     */
    tpsa() {
        const ret = wasm.molhandle_tpsa(this.__wbg_ptr);
        return ret;
    }
    /**
     * Returns `true` if the molecule passes Veber's oral bioavailability criteria
     * (TPSA ≤ 140 Å² and rotatable bonds ≤ 10).
     * @returns {boolean}
     */
    veber_passes() {
        const ret = wasm.molhandle_veber_passes(this.__wbg_ptr);
        return ret !== 0;
    }
    /**
     * Wiener topological index (sum of all pairwise shortest-path distances).
     * @returns {number}
     */
    wiener_index() {
        const ret = wasm.molhandle_wiener_index(this.__wbg_ptr);
        return ret;
    }
}
if (Symbol.dispose) MolHandle.prototype[Symbol.dispose] = MolHandle.prototype.free;

/**
 * Return a copy of the molecule with all implicit hydrogens converted to explicit H atoms.
 * @param {MolHandle} mol
 * @returns {MolHandle}
 */
export function add_hydrogens(mol) {
    _assertClass(mol, MolHandle);
    const ret = wasm.add_hydrogens(mol.__wbg_ptr);
    return MolHandle.__wrap(ret);
}

/**
 * Number of BRICS fragments produced by fragmenting the molecule.
 *
 * Returns 1 if no BRICS-breakable bonds exist (whole molecule is one fragment).
 * @param {MolHandle} mol
 * @returns {number}
 */
export function brics_fragment_count(mol) {
    _assertClass(mol, MolHandle);
    const ret = wasm.brics_fragment_count(mol.__wbg_ptr);
    return ret >>> 0;
}

/**
 * Render a reaction SMILES string (e.g. `"CC(=O)O.CCO>>CC(=O)OCC.O"`) as a
 * single SVG showing reactants → products with `+` separators.
 *
 * Returns a self-contained SVG string.  Returns a JS error on invalid input.
 * @param {string} rxn_smiles
 * @returns {string}
 */
export function depict_reaction_svg(rxn_smiles) {
    let deferred3_0;
    let deferred3_1;
    try {
        const ptr0 = passStringToWasm0(rxn_smiles, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.depict_reaction_svg(ptr0, len0);
        var ptr2 = ret[0];
        var len2 = ret[1];
        if (ret[3]) {
            ptr2 = 0; len2 = 0;
            throw takeFromExternrefTable0(ret[2]);
        }
        deferred3_0 = ptr2;
        deferred3_1 = len2;
        return getStringFromWasm0(ptr2, len2);
    } finally {
        wasm.__wbindgen_free(deferred3_0, deferred3_1, 1);
    }
}

/**
 * Render a grid SVG from newline-separated SMILES (one per line).
 *
 * Lines that fail to parse are silently skipped.
 * `cols` controls the number of columns (each cell is 200×200 px).
 * @param {string} smiles_block
 * @param {number} cols
 * @returns {string}
 */
export function depict_svg_grid(smiles_block, cols) {
    let deferred2_0;
    let deferred2_1;
    try {
        const ptr0 = passStringToWasm0(smiles_block, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.depict_svg_grid(ptr0, len0, cols);
        deferred2_0 = ret[0];
        deferred2_1 = ret[1];
        return getStringFromWasm0(ret[0], ret[1]);
    } finally {
        wasm.__wbindgen_free(deferred2_0, deferred2_1, 1);
    }
}

/**
 * Detect named functional groups in `mol`.
 *
 * Returns a JSON array of `{"name":"hydroxyl","atoms":[3]}` objects.
 * Multiple matches of the same group (e.g. two hydroxyl groups) each appear
 * as a separate entry.  Overlapping groups (carboxylic acid → "carboxyl" +
 * "hydroxyl" + "carbonyl") are all returned.
 * @param {MolHandle} mol
 * @returns {string}
 */
export function detect_functional_groups(mol) {
    let deferred1_0;
    let deferred1_1;
    try {
        _assertClass(mol, MolHandle);
        const ret = wasm.detect_functional_groups(mol.__wbg_ptr);
        deferred1_0 = ret[0];
        deferred1_1 = ret[1];
        return getStringFromWasm0(ret[0], ret[1]);
    } finally {
        wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
    }
}

/**
 * Compute the ECFP4 fingerprint as a bit-packed byte vector (256 bytes = 2048 bits).
 * @param {MolHandle} mol
 * @returns {Uint8Array}
 */
export function ecfp4_bitvec(mol) {
    _assertClass(mol, MolHandle);
    const ret = wasm.ecfp4_bitvec(mol.__wbg_ptr);
    var v1 = getArrayU8FromWasm0(ret[0], ret[1]).slice();
    wasm.__wbindgen_free(ret[0], ret[1] * 1, 1);
    return v1;
}

/**
 * Per-atom EState values as a JSON array of f64.
 *
 * Indices match `mol.atoms()` order.  Hydrogen atoms get 0.0.
 * @param {MolHandle} mol
 * @returns {string}
 */
export function estate_indices_json(mol) {
    let deferred1_0;
    let deferred1_1;
    try {
        _assertClass(mol, MolHandle);
        const ret = wasm.estate_indices_json(mol.__wbg_ptr);
        deferred1_0 = ret[0];
        deferred1_1 = ret[1];
        return getStringFromWasm0(ret[0], ret[1]);
    } finally {
        wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
    }
}

/**
 * Gasteiger-Marsili PEOE partial charges as a JSON array of f64.
 * @param {MolHandle} mol
 * @returns {string}
 */
export function gasteiger_charges_json(mol) {
    let deferred1_0;
    let deferred1_1;
    try {
        _assertClass(mol, MolHandle);
        const ret = wasm.gasteiger_charges_json(mol.__wbg_ptr);
        deferred1_0 = ret[0];
        deferred1_1 = ret[1];
        return getStringFromWasm0(ret[0], ret[1]);
    } finally {
        wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
    }
}

/**
 * Generate 3D coordinates for the molecule and return a PDB string.
 *
 * Coordinates are generated using distance-geometry placement with ring templates.
 * Returns heavy-atom PDB (HETATM records, no explicit H).
 * @param {MolHandle} mol
 * @returns {string}
 */
export function generate_3d_pdb(mol) {
    let deferred1_0;
    let deferred1_1;
    try {
        _assertClass(mol, MolHandle);
        const ret = wasm.generate_3d_pdb(mol.__wbg_ptr);
        deferred1_0 = ret[0];
        deferred1_1 = ret[1];
        return getStringFromWasm0(ret[0], ret[1]);
    } finally {
        wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
    }
}

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
 * @param {MolHandle} mol
 * @param {number} idx
 * @returns {string}
 */
export function get_atom_info(mol, idx) {
    let deferred1_0;
    let deferred1_1;
    try {
        _assertClass(mol, MolHandle);
        const ret = wasm.get_atom_info(mol.__wbg_ptr, idx);
        deferred1_0 = ret[0];
        deferred1_1 = ret[1];
        return getStringFromWasm0(ret[0], ret[1]);
    } finally {
        wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
    }
}

/**
 * Return bond information as a JSON object, looked up by the two bonded atom indices.
 *
 * Useful when you know the atom indices from SMARTS matching or `data-atom-idx` SVG
 * attributes but not the bond index.  Returns `"null"` if no bond exists between them.
 *
 * Fields: same as `get_bond_info` plus `bondIdx` (u32).
 * @param {MolHandle} mol
 * @param {number} atom1
 * @param {number} atom2
 * @returns {string}
 */
export function get_bond_between(mol, atom1, atom2) {
    let deferred1_0;
    let deferred1_1;
    try {
        _assertClass(mol, MolHandle);
        const ret = wasm.get_bond_between(mol.__wbg_ptr, atom1, atom2);
        deferred1_0 = ret[0];
        deferred1_1 = ret[1];
        return getStringFromWasm0(ret[0], ret[1]);
    } finally {
        wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
    }
}

/**
 * Return bond information as a JSON object, looked up by bond index.
 *
 * `idx` is the 0-based bond index (order matches `mol.bonds()` iteration).
 * Returns `"null"` if `idx` is out of range.
 *
 * Fields: `bondOrder` (1.0/1.5/2.0/3.0), `isAromatic` (bool),
 * `isInRing` (bool), `atomFrom` (u32), `atomTo` (u32).
 * @param {MolHandle} mol
 * @param {number} idx
 * @returns {string}
 */
export function get_bond_info(mol, idx) {
    let deferred1_0;
    let deferred1_1;
    try {
        _assertClass(mol, MolHandle);
        const ret = wasm.get_bond_info(mol.__wbg_ptr, idx);
        deferred1_0 = ret[0];
        deferred1_1 = ret[1];
        return getStringFromWasm0(ret[0], ret[1]);
    } finally {
        wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
    }
}

/**
 * Identify functional groups. Returns a JSON array of objects:
 * `[{"atoms":[0,2,3],"type":"C,N,O"}, …]`
 * @param {MolHandle} mol
 * @returns {string}
 */
export function identify_functional_groups(mol) {
    let deferred1_0;
    let deferred1_1;
    try {
        _assertClass(mol, MolHandle);
        const ret = wasm.identify_functional_groups(mol.__wbg_ptr);
        deferred1_0 = ret[0];
        deferred1_1 = ret[1];
        return getStringFromWasm0(ret[0], ret[1]);
    } finally {
        wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
    }
}

/**
 * Returns `true` if the SMILES string can be parsed without error.
 * @param {string} s
 * @returns {boolean}
 */
export function is_valid_smiles(s) {
    const ptr0 = passStringToWasm0(s, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
    const len0 = WASM_VECTOR_LEN;
    const ret = wasm.is_valid_smiles(ptr0, len0);
    return ret !== 0;
}

/**
 * Find all SMARTS matches in a molecule given only SMILES strings.
 *
 * Convenience wrapper around `smarts_match_atoms` that accepts raw SMILES
 * instead of a `MolHandle`.  Returns the same JSON format: `[[0,1],[3,4]]`.
 * Returns a JS error on SMILES or SMARTS parse failure.
 * @param {string} smiles
 * @param {string} smarts
 * @returns {string}
 */
export function match_smarts_smiles(smiles, smarts) {
    let deferred4_0;
    let deferred4_1;
    try {
        const ptr0 = passStringToWasm0(smiles, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ptr1 = passStringToWasm0(smarts, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len1 = WASM_VECTOR_LEN;
        const ret = wasm.match_smarts_smiles(ptr0, len0, ptr1, len1);
        var ptr3 = ret[0];
        var len3 = ret[1];
        if (ret[3]) {
            ptr3 = 0; len3 = 0;
            throw takeFromExternrefTable0(ret[2]);
        }
        deferred4_0 = ptr3;
        deferred4_1 = len3;
        return getStringFromWasm0(ptr3, len3);
    } finally {
        wasm.__wbindgen_free(deferred4_0, deferred4_1, 1);
    }
}

/**
 * Serialize a SMILES string directly to a MOL V2000 block.
 *
 * Convenience wrapper; all atom coordinates are 0.0.
 * Returns a JS error on SMILES parse failure.
 * @param {string} smiles
 * @returns {string}
 */
export function mol_block_from_smiles(smiles) {
    let deferred3_0;
    let deferred3_1;
    try {
        const ptr0 = passStringToWasm0(smiles, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.mol_block_from_smiles(ptr0, len0);
        var ptr2 = ret[0];
        var len2 = ret[1];
        if (ret[3]) {
            ptr2 = 0; len2 = 0;
            throw takeFromExternrefTable0(ret[2]);
        }
        deferred3_0 = ptr2;
        deferred3_1 = len2;
        return getStringFromWasm0(ptr2, len2);
    } finally {
        wasm.__wbindgen_free(deferred3_0, deferred3_1, 1);
    }
}

/**
 * Parse a MOL V2000 block and return a `MolHandle`.
 *
 * Returns a JS error string on parse failure.
 * @param {string} block
 * @returns {MolHandle}
 */
export function mol_from_sdf_block(block) {
    const ptr0 = passStringToWasm0(block, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
    const len0 = WASM_VECTOR_LEN;
    const ret = wasm.mol_from_sdf_block(ptr0, len0);
    if (ret[2]) {
        throw takeFromExternrefTable0(ret[1]);
    }
    return MolHandle.__wrap(ret[0]);
}

/**
 * Parse a SMILES string into a `MolHandle`.
 *
 * Returns a JS error string on parse failure.
 * @param {string} s
 * @returns {MolHandle}
 */
export function parse_smiles(s) {
    const ptr0 = passStringToWasm0(s, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
    const len0 = WASM_VECTOR_LEN;
    const ret = wasm.parse_smiles(ptr0, len0);
    if (ret[2]) {
        throw takeFromExternrefTable0(ret[1]);
    }
    return MolHandle.__wrap(ret[0]);
}

/**
 * PEOE_VSA descriptors (14 bins) as a JSON array.
 * @param {MolHandle} mol
 * @returns {string}
 */
export function peoe_vsa_json(mol) {
    let deferred1_0;
    let deferred1_1;
    try {
        _assertClass(mol, MolHandle);
        const ret = wasm.peoe_vsa_json(mol.__wbg_ptr);
        deferred1_0 = ret[0];
        deferred1_1 = ret[1];
        return getStringFromWasm0(ret[0], ret[1]);
    } finally {
        wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
    }
}

/**
 * Return a copy of the molecule with all explicit hydrogen atoms removed.
 * @param {MolHandle} mol
 * @returns {MolHandle}
 */
export function remove_hydrogens(mol) {
    _assertClass(mol, MolHandle);
    const ret = wasm.remove_hydrogens(mol.__wbg_ptr);
    return MolHandle.__wrap(ret);
}

/**
 * Apply a SMIRKS reaction template and return product SMILES as a JSON string.
 *
 * `reactants_smiles`: pipe-separated SMILES, one per reactant slot in the SMIRKS.
 * Returns a JSON array of arrays: `[["product_smi", …], …]`.
 * Returns a JS error on parse failure or arity mismatch.
 * @param {string} smirks
 * @param {string} reactants_smiles
 * @returns {string}
 */
export function run_reactants(smirks, reactants_smiles) {
    let deferred4_0;
    let deferred4_1;
    try {
        const ptr0 = passStringToWasm0(smirks, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ptr1 = passStringToWasm0(reactants_smiles, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len1 = WASM_VECTOR_LEN;
        const ret = wasm.run_reactants(ptr0, len0, ptr1, len1);
        var ptr3 = ret[0];
        var len3 = ret[1];
        if (ret[3]) {
            ptr3 = 0; len3 = 0;
            throw takeFromExternrefTable0(ret[2]);
        }
        deferred4_0 = ptr3;
        deferred4_1 = len3;
        return getStringFromWasm0(ptr3, len3);
    } finally {
        wasm.__wbindgen_free(deferred4_0, deferred4_1, 1);
    }
}

/**
 * Synthetic Accessibility Score (1 = easy, 10 = hard).
 * @param {MolHandle} mol
 * @returns {number}
 */
export function sa_score(mol) {
    _assertClass(mol, MolHandle);
    const ret = wasm.sa_score(mol.__wbg_ptr);
    return ret;
}

/**
 * Parse an SDF string and return a JSON array of canonical SMILES strings.
 *
 * Invalid records are represented as `null` in the array.
 * @param {string} sdf
 * @returns {string}
 */
export function sdf_to_smiles_json(sdf) {
    let deferred2_0;
    let deferred2_1;
    try {
        const ptr0 = passStringToWasm0(sdf, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.sdf_to_smiles_json(ptr0, len0);
        deferred2_0 = ret[0];
        deferred2_1 = ret[1];
        return getStringFromWasm0(ret[0], ret[1]);
    } finally {
        wasm.__wbindgen_free(deferred2_0, deferred2_1, 1);
    }
}

/**
 * SlogP_VSA descriptors (12 bins) as a JSON array.
 * @param {MolHandle} mol
 * @returns {string}
 */
export function slogp_vsa_json(mol) {
    let deferred1_0;
    let deferred1_1;
    try {
        _assertClass(mol, MolHandle);
        const ret = wasm.slogp_vsa_json(mol.__wbg_ptr);
        deferred1_0 = ret[0];
        deferred1_1 = ret[1];
        return getStringFromWasm0(ret[0], ret[1]);
    } finally {
        wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
    }
}

/**
 * Find all substructure matches of a SMARTS pattern in `mol`.
 *
 * Returns JSON array of arrays of atom indices (sorted, 0-based).
 * Example: `[[0,1,2],[3,4,5]]` — two matches.
 * Returns `"[]"` if no match. Returns a JS error on invalid SMARTS.
 * @param {string} smarts
 * @param {MolHandle} mol
 * @returns {string}
 */
export function smarts_match_atoms(smarts, mol) {
    let deferred3_0;
    let deferred3_1;
    try {
        const ptr0 = passStringToWasm0(smarts, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        _assertClass(mol, MolHandle);
        const ret = wasm.smarts_match_atoms(ptr0, len0, mol.__wbg_ptr);
        var ptr2 = ret[0];
        var len2 = ret[1];
        if (ret[3]) {
            ptr2 = 0; len2 = 0;
            throw takeFromExternrefTable0(ret[2]);
        }
        deferred3_0 = ptr2;
        deferred3_1 = len2;
        return getStringFromWasm0(ptr2, len2);
    } finally {
        wasm.__wbindgen_free(deferred3_0, deferred3_1, 1);
    }
}

/**
 * Render a highlighted SVG from a SMILES string in one call.
 *
 * `atoms` — 0-based atom indices to highlight (Uint32Array in JS).
 * `bonds` — 0-based bond indices to highlight (Uint32Array in JS).
 * `color` — CSS color for highlights (e.g. `"#ef4444"`); empty string uses default yellow.
 *
 * Returns a JS error on SMILES parse failure.
 * @param {string} smiles
 * @param {Uint32Array} atoms
 * @param {Uint32Array} bonds
 * @param {string} color
 * @returns {string}
 */
export function smiles_to_svg_highlighted(smiles, atoms, bonds, color) {
    let deferred6_0;
    let deferred6_1;
    try {
        const ptr0 = passStringToWasm0(smiles, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ptr1 = passArray32ToWasm0(atoms, wasm.__wbindgen_malloc);
        const len1 = WASM_VECTOR_LEN;
        const ptr2 = passArray32ToWasm0(bonds, wasm.__wbindgen_malloc);
        const len2 = WASM_VECTOR_LEN;
        const ptr3 = passStringToWasm0(color, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len3 = WASM_VECTOR_LEN;
        const ret = wasm.smiles_to_svg_highlighted(ptr0, len0, ptr1, len1, ptr2, len2, ptr3, len3);
        var ptr5 = ret[0];
        var len5 = ret[1];
        if (ret[3]) {
            ptr5 = 0; len5 = 0;
            throw takeFromExternrefTable0(ret[2]);
        }
        deferred6_0 = ptr5;
        deferred6_1 = len5;
        return getStringFromWasm0(ptr5, len5);
    } finally {
        wasm.__wbindgen_free(deferred6_0, deferred6_1, 1);
    }
}

/**
 * SMR_VSA descriptors (10 bins) as a JSON array.
 * @param {MolHandle} mol
 * @returns {string}
 */
export function smr_vsa_json(mol) {
    let deferred1_0;
    let deferred1_1;
    try {
        _assertClass(mol, MolHandle);
        const ret = wasm.smr_vsa_json(mol.__wbg_ptr);
        deferred1_0 = ret[0];
        deferred1_1 = ret[1];
        return getStringFromWasm0(ret[0], ret[1]);
    } finally {
        wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
    }
}

export function start() {
    wasm.start();
}

/**
 * Tanimoto similarity between two molecules using AtomPair fingerprints.
 * @param {MolHandle} a
 * @param {MolHandle} b
 * @returns {number}
 */
export function tanimoto_atom_pair(a, b) {
    _assertClass(a, MolHandle);
    _assertClass(b, MolHandle);
    const ret = wasm.tanimoto_atom_pair(a.__wbg_ptr, b.__wbg_ptr);
    return ret;
}

/**
 * Tanimoto similarity between two molecules using ECFP4 fingerprints.
 * @param {MolHandle} a
 * @param {MolHandle} b
 * @returns {number}
 */
export function tanimoto_ecfp4(a, b) {
    _assertClass(a, MolHandle);
    _assertClass(b, MolHandle);
    const ret = wasm.tanimoto_ecfp4(a.__wbg_ptr, b.__wbg_ptr);
    return ret;
}

/**
 * Tanimoto similarity between two molecules using FCFP4 fingerprints (pharmacophore-based).
 * @param {MolHandle} a
 * @param {MolHandle} b
 * @returns {number}
 */
export function tanimoto_fcfp4(a, b) {
    _assertClass(a, MolHandle);
    _assertClass(b, MolHandle);
    const ret = wasm.tanimoto_fcfp4(a.__wbg_ptr, b.__wbg_ptr);
    return ret;
}

/**
 * Tanimoto similarity between two molecules given only SMILES strings (ECFP4).
 *
 * Returns a JS error on parse failure.
 * @param {string} smiles1
 * @param {string} smiles2
 * @returns {number}
 */
export function tanimoto_smiles(smiles1, smiles2) {
    const ptr0 = passStringToWasm0(smiles1, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
    const len0 = WASM_VECTOR_LEN;
    const ptr1 = passStringToWasm0(smiles2, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
    const len1 = WASM_VECTOR_LEN;
    const ret = wasm.tanimoto_smiles(ptr0, len0, ptr1, len1);
    if (ret[2]) {
        throw takeFromExternrefTable0(ret[1]);
    }
    return ret[0];
}

/**
 * Tanimoto similarity between two molecules using topological path fingerprints.
 * @param {MolHandle} a
 * @param {MolHandle} b
 * @returns {number}
 */
export function tanimoto_topo_path(a, b) {
    _assertClass(a, MolHandle);
    _assertClass(b, MolHandle);
    const ret = wasm.tanimoto_topo_path(a.__wbg_ptr, b.__wbg_ptr);
    return ret;
}

/**
 * Tanimoto similarity between two molecules using Topological Torsion fingerprints.
 * @param {MolHandle} a
 * @param {MolHandle} b
 * @returns {number}
 */
export function tanimoto_torsion(a, b) {
    _assertClass(a, MolHandle);
    _assertClass(b, MolHandle);
    const ret = wasm.tanimoto_torsion(a.__wbg_ptr, b.__wbg_ptr);
    return ret;
}

/**
 * Serialize a molecule to a MOL V2000 block.
 *
 * All atom coordinates are written as 0.0 (the `Molecule` type has no 2D
 * coordinate storage; real coordinates would require a separate layout pass).
 * @param {MolHandle} mol
 * @returns {string}
 */
export function to_mol_block(mol) {
    let deferred1_0;
    let deferred1_1;
    try {
        _assertClass(mol, MolHandle);
        const ret = wasm.to_mol_block(mol.__wbg_ptr);
        deferred1_0 = ret[0];
        deferred1_1 = ret[1];
        return getStringFromWasm0(ret[0], ret[1]);
    } finally {
        wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
    }
}
function __wbg_get_imports() {
    const import0 = {
        __proto__: null,
        __wbg___wbindgen_throw_1506f2235d1bdba0: function(arg0, arg1) {
            throw new Error(getStringFromWasm0(arg0, arg1));
        },
        __wbg_error_a6fa202b58aa1cd3: function(arg0, arg1) {
            let deferred0_0;
            let deferred0_1;
            try {
                deferred0_0 = arg0;
                deferred0_1 = arg1;
                console.error(getStringFromWasm0(arg0, arg1));
            } finally {
                wasm.__wbindgen_free(deferred0_0, deferred0_1, 1);
            }
        },
        __wbg_new_227d7c05414eb861: function() {
            const ret = new Error();
            return ret;
        },
        __wbg_stack_3b0d974bbf31e44f: function(arg0, arg1) {
            const ret = arg1.stack;
            const ptr1 = passStringToWasm0(ret, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
            const len1 = WASM_VECTOR_LEN;
            getDataViewMemory0().setInt32(arg0 + 4 * 1, len1, true);
            getDataViewMemory0().setInt32(arg0 + 4 * 0, ptr1, true);
        },
        __wbindgen_cast_0000000000000001: function(arg0, arg1) {
            // Cast intrinsic for `Ref(String) -> Externref`.
            const ret = getStringFromWasm0(arg0, arg1);
            return ret;
        },
        __wbindgen_init_externref_table: function() {
            const table = wasm.__wbindgen_externrefs;
            const offset = table.grow(4);
            table.set(0, undefined);
            table.set(offset + 0, undefined);
            table.set(offset + 1, null);
            table.set(offset + 2, true);
            table.set(offset + 3, false);
        },
    };
    return {
        __proto__: null,
        "./chematic_wasm_bg.js": import0,
    };
}

const DepictOptionsFinalization = (typeof FinalizationRegistry === 'undefined')
    ? { register: () => {}, unregister: () => {} }
    : new FinalizationRegistry(ptr => wasm.__wbg_depictoptions_free(ptr, 1));
const MolHandleFinalization = (typeof FinalizationRegistry === 'undefined')
    ? { register: () => {}, unregister: () => {} }
    : new FinalizationRegistry(ptr => wasm.__wbg_molhandle_free(ptr, 1));

function _assertClass(instance, klass) {
    if (!(instance instanceof klass)) {
        throw new Error(`expected instance of ${klass.name}`);
    }
}

function getArrayU8FromWasm0(ptr, len) {
    ptr = ptr >>> 0;
    return getUint8ArrayMemory0().subarray(ptr / 1, ptr / 1 + len);
}

let cachedDataViewMemory0 = null;
function getDataViewMemory0() {
    if (cachedDataViewMemory0 === null || cachedDataViewMemory0.buffer.detached === true || (cachedDataViewMemory0.buffer.detached === undefined && cachedDataViewMemory0.buffer !== wasm.memory.buffer)) {
        cachedDataViewMemory0 = new DataView(wasm.memory.buffer);
    }
    return cachedDataViewMemory0;
}

function getStringFromWasm0(ptr, len) {
    return decodeText(ptr >>> 0, len);
}

let cachedUint32ArrayMemory0 = null;
function getUint32ArrayMemory0() {
    if (cachedUint32ArrayMemory0 === null || cachedUint32ArrayMemory0.byteLength === 0) {
        cachedUint32ArrayMemory0 = new Uint32Array(wasm.memory.buffer);
    }
    return cachedUint32ArrayMemory0;
}

let cachedUint8ArrayMemory0 = null;
function getUint8ArrayMemory0() {
    if (cachedUint8ArrayMemory0 === null || cachedUint8ArrayMemory0.byteLength === 0) {
        cachedUint8ArrayMemory0 = new Uint8Array(wasm.memory.buffer);
    }
    return cachedUint8ArrayMemory0;
}

function passArray32ToWasm0(arg, malloc) {
    const ptr = malloc(arg.length * 4, 4) >>> 0;
    getUint32ArrayMemory0().set(arg, ptr / 4);
    WASM_VECTOR_LEN = arg.length;
    return ptr;
}

function passStringToWasm0(arg, malloc, realloc) {
    if (realloc === undefined) {
        const buf = cachedTextEncoder.encode(arg);
        const ptr = malloc(buf.length, 1) >>> 0;
        getUint8ArrayMemory0().subarray(ptr, ptr + buf.length).set(buf);
        WASM_VECTOR_LEN = buf.length;
        return ptr;
    }

    let len = arg.length;
    let ptr = malloc(len, 1) >>> 0;

    const mem = getUint8ArrayMemory0();

    let offset = 0;

    for (; offset < len; offset++) {
        const code = arg.charCodeAt(offset);
        if (code > 0x7F) break;
        mem[ptr + offset] = code;
    }
    if (offset !== len) {
        if (offset !== 0) {
            arg = arg.slice(offset);
        }
        ptr = realloc(ptr, len, len = offset + arg.length * 3, 1) >>> 0;
        const view = getUint8ArrayMemory0().subarray(ptr + offset, ptr + len);
        const ret = cachedTextEncoder.encodeInto(arg, view);

        offset += ret.written;
        ptr = realloc(ptr, len, offset, 1) >>> 0;
    }

    WASM_VECTOR_LEN = offset;
    return ptr;
}

function takeFromExternrefTable0(idx) {
    const value = wasm.__wbindgen_externrefs.get(idx);
    wasm.__externref_table_dealloc(idx);
    return value;
}

let cachedTextDecoder = new TextDecoder('utf-8', { ignoreBOM: true, fatal: true });
cachedTextDecoder.decode();
const MAX_SAFARI_DECODE_BYTES = 2146435072;
let numBytesDecoded = 0;
function decodeText(ptr, len) {
    numBytesDecoded += len;
    if (numBytesDecoded >= MAX_SAFARI_DECODE_BYTES) {
        cachedTextDecoder = new TextDecoder('utf-8', { ignoreBOM: true, fatal: true });
        cachedTextDecoder.decode();
        numBytesDecoded = len;
    }
    return cachedTextDecoder.decode(getUint8ArrayMemory0().subarray(ptr, ptr + len));
}

const cachedTextEncoder = new TextEncoder();

if (!('encodeInto' in cachedTextEncoder)) {
    cachedTextEncoder.encodeInto = function (arg, view) {
        const buf = cachedTextEncoder.encode(arg);
        view.set(buf);
        return {
            read: arg.length,
            written: buf.length
        };
    };
}

let WASM_VECTOR_LEN = 0;

let wasmModule, wasmInstance, wasm;
function __wbg_finalize_init(instance, module) {
    wasmInstance = instance;
    wasm = instance.exports;
    wasmModule = module;
    cachedDataViewMemory0 = null;
    cachedUint32ArrayMemory0 = null;
    cachedUint8ArrayMemory0 = null;
    wasm.__wbindgen_start();
    return wasm;
}

async function __wbg_load(module, imports) {
    if (typeof Response === 'function' && module instanceof Response) {
        if (typeof WebAssembly.instantiateStreaming === 'function') {
            try {
                return await WebAssembly.instantiateStreaming(module, imports);
            } catch (e) {
                const validResponse = module.ok && expectedResponseType(module.type);

                if (validResponse && module.headers.get('Content-Type') !== 'application/wasm') {
                    console.warn("`WebAssembly.instantiateStreaming` failed because your server does not serve Wasm with `application/wasm` MIME type. Falling back to `WebAssembly.instantiate` which is slower. Original error:\n", e);

                } else { throw e; }
            }
        }

        const bytes = await module.arrayBuffer();
        return await WebAssembly.instantiate(bytes, imports);
    } else {
        const instance = await WebAssembly.instantiate(module, imports);

        if (instance instanceof WebAssembly.Instance) {
            return { instance, module };
        } else {
            return instance;
        }
    }

    function expectedResponseType(type) {
        switch (type) {
            case 'basic': case 'cors': case 'default': return true;
        }
        return false;
    }
}

function initSync(module) {
    if (wasm !== undefined) return wasm;


    if (module !== undefined) {
        if (Object.getPrototypeOf(module) === Object.prototype) {
            ({module} = module)
        } else {
            console.warn('using deprecated parameters for `initSync()`; pass a single object instead')
        }
    }

    const imports = __wbg_get_imports();
    if (!(module instanceof WebAssembly.Module)) {
        module = new WebAssembly.Module(module);
    }
    const instance = new WebAssembly.Instance(module, imports);
    return __wbg_finalize_init(instance, module);
}

async function __wbg_init(module_or_path) {
    if (wasm !== undefined) return wasm;


    if (module_or_path !== undefined) {
        if (Object.getPrototypeOf(module_or_path) === Object.prototype) {
            ({module_or_path} = module_or_path)
        } else {
            console.warn('using deprecated parameters for the initialization function; pass a single object instead')
        }
    }

    if (module_or_path === undefined) {
        module_or_path = new URL('chematic_wasm_bg.wasm', import.meta.url);
    }
    const imports = __wbg_get_imports();

    if (typeof module_or_path === 'string' || (typeof Request === 'function' && module_or_path instanceof Request) || (typeof URL === 'function' && module_or_path instanceof URL)) {
        module_or_path = fetch(module_or_path);
    }

    const { instance, module } = await __wbg_load(await module_or_path, imports);

    return __wbg_finalize_init(instance, module);
}

export { initSync, __wbg_init as default };
