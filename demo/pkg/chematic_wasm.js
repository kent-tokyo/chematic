/* @ts-self-types="./chematic_wasm.d.ts" */

/**
 * A conformer ensemble: one molecule geometry with multiple 3D coordinate sets.
 *
 * Create with `new(smiles)`, then add conformers with `add_generated_conformer`
 * or `add_minimized_conformer`.  Retrieve coordinates as PDB strings via
 * `get_conformer_pdb(idx)`.  Compare conformers with `conformer_rmsd`.
 */
export class ConformerHandle {
    __destroy_into_raw() {
        const ptr = this.__wbg_ptr;
        this.__wbg_ptr = 0;
        ConformerHandleFinalization.unregister(this);
        return ptr;
    }
    free() {
        const ptr = this.__destroy_into_raw();
        wasm.__wbg_conformerhandle_free(ptr, 0);
    }
    /**
     * Generate a new 3D conformer using distance-geometry and add it to the ensemble.
     *
     * Returns the index of the newly added conformer.
     * @returns {number}
     */
    add_generated_conformer() {
        const ret = wasm.conformerhandle_add_generated_conformer(this.__wbg_ptr);
        return ret >>> 0;
    }
    /**
     * Generate a new 3D conformer, run force-field minimization, and add it.
     *
     * Returns the index of the newly added conformer.
     * @returns {number}
     */
    add_minimized_conformer() {
        const ret = wasm.conformerhandle_add_minimized_conformer(this.__wbg_ptr);
        return ret >>> 0;
    }
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
     * @param {number} rms_threshold
     * @returns {string}
     */
    cluster_conformers_json(rms_threshold) {
        let deferred1_0;
        let deferred1_1;
        try {
            const ret = wasm.conformerhandle_cluster_conformers_json(this.__wbg_ptr, rms_threshold);
            deferred1_0 = ret[0];
            deferred1_1 = ret[1];
            return getStringFromWasm0(ret[0], ret[1]);
        } finally {
            wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
        }
    }
    /**
     * Number of conformers currently stored.
     * @returns {number}
     */
    conformer_count() {
        const ret = wasm.conformerhandle_conformer_count(this.__wbg_ptr);
        return ret >>> 0;
    }
    /**
     * Kabsch-aligned RMSD (Å) between conformers `a` and `b`.
     *
     * Returns `NaN` if either index is out of range.
     * @param {number} a
     * @param {number} b
     * @returns {number}
     */
    conformer_rmsd(a, b) {
        const ret = wasm.conformerhandle_conformer_rmsd(this.__wbg_ptr, a, b);
        return ret;
    }
    /**
     * Un-aligned (translation + rotation NOT removed) RMSD (Å) between conformers `a` and `b`.
     *
     * Returns `NaN` if either index is out of range.
     * @param {number} a
     * @param {number} b
     * @returns {number}
     */
    conformer_rmsd_no_align(a, b) {
        const ret = wasm.conformerhandle_conformer_rmsd_no_align(this.__wbg_ptr, a, b);
        return ret;
    }
    /**
     * Return conformer `idx` as a PDB string, or `null` if `idx` is out of range.
     * @param {number} idx
     * @returns {string | undefined}
     */
    get_conformer_pdb(idx) {
        const ret = wasm.conformerhandle_get_conformer_pdb(this.__wbg_ptr, idx);
        let v1;
        if (ret[0] !== 0) {
            v1 = getStringFromWasm0(ret[0], ret[1]).slice();
            wasm.__wbindgen_free(ret[0], ret[1] * 1, 1);
        }
        return v1;
    }
    /**
     * The ensemble's molecule as a `MolHandle`.
     * @returns {MolHandle}
     */
    mol() {
        const ret = wasm.conformerhandle_mol(this.__wbg_ptr);
        return MolHandle.__wrap(ret);
    }
    /**
     * Create a new empty ensemble for the molecule given by `smiles`.
     *
     * Returns a JS error on SMILES parse failure.
     * @param {string} smiles
     */
    constructor(smiles) {
        const ptr0 = passStringToWasm0(smiles, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.conformerhandle_new(ptr0, len0);
        if (ret[2]) {
            throw takeFromExternrefTable0(ret[1]);
        }
        this.__wbg_ptr = ret[0];
        ConformerHandleFinalization.register(this, this.__wbg_ptr, this);
        return this;
    }
    /**
     * Remove conformer `idx` and return `true`, or `false` if `idx` is out of range.
     * @param {number} idx
     * @returns {boolean}
     */
    remove_conformer(idx) {
        const ret = wasm.conformerhandle_remove_conformer(this.__wbg_ptr, idx);
        return ret !== 0;
    }
}
if (Symbol.dispose) ConformerHandle.prototype[Symbol.dispose] = ConformerHandle.prototype.free;

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

export class MhfpLshHandle {
    __destroy_into_raw() {
        const ptr = this.__wbg_ptr;
        this.__wbg_ptr = 0;
        MhfpLshHandleFinalization.unregister(this);
        return ptr;
    }
    free() {
        const ptr = this.__destroy_into_raw();
        wasm.__wbg_mhfplshhandle_free(ptr, 0);
    }
    /**
     * Add a molecule by SMILES; returns its 0-based index in the index.
     * @param {string} smiles
     * @returns {number}
     */
    add_smiles(smiles) {
        const ptr0 = passStringToWasm0(smiles, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.mhfplshhandle_add_smiles(this.__wbg_ptr, ptr0, len0);
        if (ret[2]) {
            throw takeFromExternrefTable0(ret[1]);
        }
        return ret[0] >>> 0;
    }
    /**
     * True if the index contains no molecules.
     * @returns {boolean}
     */
    is_empty() {
        const ret = wasm.mhfplshhandle_is_empty(this.__wbg_ptr);
        return ret !== 0;
    }
    /**
     * Number of molecules in the index.
     * @returns {number}
     */
    len() {
        const ret = wasm.mhfplshhandle_len(this.__wbg_ptr);
        return ret >>> 0;
    }
    /**
     * Create a new LSH index for MHFP fingerprints with `num_hashes` hash lanes.
     * Default band decomposition: 16 bands × (num_hashes / 16) rows.
     * `num_hashes` must be a multiple of 16 (e.g. 128).
     * @param {number} num_hashes
     */
    constructor(num_hashes) {
        const ret = wasm.mhfplshhandle_new(num_hashes);
        this.__wbg_ptr = ret;
        MhfpLshHandleFinalization.register(this, this.__wbg_ptr, this);
        return this;
    }
    /**
     * Query by SMILES for all entries with similarity ≥ threshold.
     *
     * Returns a JSON array `[{"index":N,"similarity":0.xxx},...]` sorted by
     * descending similarity.  Empty array `[]` when nothing qualifies.
     * @param {string} query_smiles
     * @param {number} threshold
     * @returns {string}
     */
    query_json(query_smiles, threshold) {
        let deferred3_0;
        let deferred3_1;
        try {
            const ptr0 = passStringToWasm0(query_smiles, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
            const len0 = WASM_VECTOR_LEN;
            const ret = wasm.mhfplshhandle_query_json(this.__wbg_ptr, ptr0, len0, threshold);
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
}
if (Symbol.dispose) MhfpLshHandle.prototype[Symbol.dispose] = MhfpLshHandle.prototype.free;

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
     * Assign CIP (R/S/E/Z) stereocenters and return JSON.
     *
     * Format: `{"centers":[{"atom":0,"code":"R"},{"atom":3,"code":"E"}]}`
     * @returns {string}
     */
    assign_cip_json() {
        let deferred1_0;
        let deferred1_1;
        try {
            const ret = wasm.molhandle_assign_cip_json(this.__wbg_ptr);
            deferred1_0 = ret[0];
            deferred1_1 = ret[1];
            return getStringFromWasm0(ret[0], ret[1]);
        } finally {
            wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
        }
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
     * Returns true when TPSA < 90 Å², MW < 400, HBD ≤ 3.
     * @returns {boolean}
     */
    bbb_passes() {
        const ret = wasm.molhandle_bbb_passes(this.__wbg_ptr);
        return ret !== 0;
    }
    /**
     * Clark (2000) blood-brain barrier logBB score.
     * logBB > −1.0 = likely CNS penetrant.
     * @returns {number}
     */
    bbb_score() {
        const ret = wasm.molhandle_bbb_score(this.__wbg_ptr);
        return ret;
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
     * Palm (1997) Caco-2 intestinal permeability (logPCaco2).
     * > −5.5 = high permeability.
     * @returns {number}
     */
    caco2_permeability() {
        const ret = wasm.molhandle_caco2_permeability(this.__wbg_ptr);
        return ret;
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
     * CYP3A4 metabolic inhibition risk score (0.0–1.0).
     * @returns {number}
     */
    cyp3a4_inhibition_risk() {
        const ret = wasm.molhandle_cyp3a4_inhibition_risk(this.__wbg_ptr);
        return ret;
    }
    /**
     * 2D PNG depiction — not available in the WASM build (PNG stack disabled to reduce bundle size).
     * Use `depict_svg()` in browser contexts; rasterize client-side if needed.
     * @returns {Uint8Array}
     */
    depict_png() {
        const ret = wasm.molhandle_depict_png(this.__wbg_ptr);
        var v1 = getArrayU8FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 1, 1);
        return v1;
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
     * hERG cardiac toxicity risk score (0.0–1.0).
     * @returns {number}
     */
    herg_risk_score() {
        const ret = wasm.molhandle_herg_risk_score(this.__wbg_ptr);
        return ret;
    }
    /**
     * Isotope distribution as JSON.
     *
     * Returns `[{"mass":100.0,"abundance":0.9},...]` sorted by mass.
     * `resolution`: m/z bin width in Da (e.g. `0.1` for nominal, `0.01` for high-res).
     * @param {number} resolution
     * @returns {string}
     */
    isotope_distribution_json(resolution) {
        let deferred1_0;
        let deferred1_1;
        try {
            const ret = wasm.molhandle_isotope_distribution_json(this.__wbg_ptr, resolution);
            deferred1_0 = ret[0];
            deferred1_1 = ret[1];
            return getStringFromWasm0(ret[0], ret[1]);
        } finally {
            wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
        }
    }
    /**
     * Generate IUPAC systematic name for the molecule.
     *
     * Returns the name string on success, or an empty string when the
     * structure is outside the supported naming scope (complex polycyclics,
     * multi-functional groups, etc.).
     * @returns {string}
     */
    iupac_name() {
        let deferred1_0;
        let deferred1_1;
        try {
            const ret = wasm.molhandle_iupac_name(this.__wbg_ptr);
            deferred1_0 = ret[0];
            deferred1_1 = ret[1];
            return getStringFromWasm0(ret[0], ret[1]);
        } finally {
            wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
        }
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
     * LogD (distribution coefficient) at a specific pH.
     *
     * Accounts for ionization state: neutral molecules return LogP unchanged,
     * ionizable molecules are adjusted by log(neutral_fraction).
     * @param {number} ph
     * @returns {number}
     */
    logd_at_ph(ph) {
        const ret = wasm.molhandle_logd_at_ph(this.__wbg_ptr, ph);
        return ret;
    }
    /**
     * LogD profile across a pH range as JSON.
     *
     * Returns `[{"ph":0.0,"logd":2.5}, ...]` with `steps` evenly-spaced pH points.
     * @param {number} ph_start
     * @param {number} ph_end
     * @param {number} steps
     * @returns {string}
     */
    logd_profile_json(ph_start, ph_end, steps) {
        let deferred1_0;
        let deferred1_1;
        try {
            const ret = wasm.molhandle_logd_profile_json(this.__wbg_ptr, ph_start, ph_end, steps);
            deferred1_0 = ret[0];
            deferred1_1 = ret[1];
            return getStringFromWasm0(ret[0], ret[1]);
        } finally {
            wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
        }
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
     * Count of aliphatic (non-aromatic) rings in the SSSR.
     * @returns {number}
     */
    num_aliphatic_rings() {
        const ret = wasm.molhandle_num_aliphatic_rings(this.__wbg_ptr);
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
     * Count of fully saturated rings in the SSSR.
     * @returns {number}
     */
    num_saturated_rings() {
        const ret = wasm.molhandle_num_saturated_rings(this.__wbg_ptr);
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
     * Count of tetrahedral stereocenters with unspecified configuration.
     * @returns {number}
     */
    num_unspecified_stereocenters() {
        const ret = wasm.molhandle_num_unspecified_stereocenters(this.__wbg_ptr);
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
     * Most acidic pKa in the molecule, or NaN if no acidic site.
     * @returns {number}
     */
    pka_acid_value() {
        const ret = wasm.molhandle_pka_acid_value(this.__wbg_ptr);
        return ret;
    }
    /**
     * Most basic pKa in the molecule, or NaN if no basic site.
     * @returns {number}
     */
    pka_base_value() {
        const ret = wasm.molhandle_pka_base_value(this.__wbg_ptr);
        return ret;
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
     * Randić connectivity index (χ₀).
     *
     * χ₀ = Σ 1/√(d_i × d_j) over all bonds, where d is heavy-atom degree.
     * @returns {number}
     */
    randic_index() {
        const ret = wasm.molhandle_randic_index(this.__wbg_ptr);
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
     * InChI string representation of the molecule.
     * @returns {string}
     */
    to_inchi() {
        let deferred1_0;
        let deferred1_1;
        try {
            const ret = wasm.molhandle_to_inchi(this.__wbg_ptr);
            deferred1_0 = ret[0];
            deferred1_1 = ret[1];
            return getStringFromWasm0(ret[0], ret[1]);
        } finally {
            wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
        }
    }
    /**
     * InChIKey (27-character identifier) for the molecule.
     * @returns {string}
     */
    to_inchikey() {
        let deferred1_0;
        let deferred1_1;
        try {
            const ret = wasm.molhandle_to_inchikey(this.__wbg_ptr);
            deferred1_0 = ret[0];
            deferred1_1 = ret[1];
            return getStringFromWasm0(ret[0], ret[1]);
        } finally {
            wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
        }
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
    /**
     * Zagreb index M1: Σ d_i² over all heavy atoms.
     * @returns {number}
     */
    zagreb_index_m1() {
        const ret = wasm.molhandle_zagreb_index_m1(this.__wbg_ptr);
        return ret >>> 0;
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
 * Compute a full ADMET property profile for a molecule.
 *
 * Returns a JSON object with fields:
 * `bbb_score`, `bbb_passes`, `caco2`, `herg_risk`, `cyp3a4_risk`,
 * `pka_acid` (null if absent), `pka_base` (null if absent),
 * `esol`, `logd74`, `mw`, `logp`, `tpsa`, `hbd`, `hba`, `rotatable_bonds`
 *
 * Returns `{"error":"..."}` on parse failure.
 * @param {string} smiles
 * @returns {string}
 */
export function admet_profile_json(smiles) {
    let deferred2_0;
    let deferred2_1;
    try {
        const ptr0 = passStringToWasm0(smiles, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.admet_profile_json(ptr0, len0);
        deferred2_0 = ret[0];
        deferred2_1 = ret[1];
        return getStringFromWasm0(ret[0], ret[1]);
    } finally {
        wasm.__wbindgen_free(deferred2_0, deferred2_1, 1);
    }
}

/**
 * AtomPair fingerprint as a bit-packed byte vector (256 bytes = 2048 bits).
 * @param {MolHandle} mol
 * @returns {Uint8Array}
 */
export function atom_pair_bitvec(mol) {
    _assertClass(mol, MolHandle);
    const ret = wasm.atom_pair_bitvec(mol.__wbg_ptr);
    var v1 = getArrayU8FromWasm0(ret[0], ret[1]).slice();
    wasm.__wbindgen_free(ret[0], ret[1] * 1, 1);
    return v1;
}

/**
 * AutoCorr2D descriptor (7 values: topological distance lags 1-7).
 * @param {MolHandle} mol
 * @returns {string}
 */
export function autocorr_2d_json(mol) {
    let deferred1_0;
    let deferred1_1;
    try {
        _assertClass(mol, MolHandle);
        const ret = wasm.autocorr_2d_json(mol.__wbg_ptr);
        deferred1_0 = ret[0];
        deferred1_1 = ret[1];
        return getStringFromWasm0(ret[0], ret[1]);
    } finally {
        wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
    }
}

/**
 * AutoCorr3D descriptor (8 values: Euclidean distance bins 1-8 Å).
 * Requires 3D coordinates (generated automatically).
 * @param {MolHandle} mol
 * @returns {string}
 */
export function autocorr_3d_json(mol) {
    let deferred1_0;
    let deferred1_1;
    try {
        _assertClass(mol, MolHandle);
        const ret = wasm.autocorr_3d_json(mol.__wbg_ptr);
        deferred1_0 = ret[0];
        deferred1_1 = ret[1];
        return getStringFromWasm0(ret[0], ret[1]);
    } finally {
        wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
    }
}

/**
 * Check whether a reaction SMILES is atom-balanced.
 *
 * Returns JSON: `{ "balanced": true|false, "diff": ["C: 1 reactant vs 2 product", ...] }`
 * Returns `"error:<msg>"` on parse failure.
 * @param {string} reaction_smiles
 * @returns {string}
 */
export function balance_check_json(reaction_smiles) {
    let deferred2_0;
    let deferred2_1;
    try {
        const ptr0 = passStringToWasm0(reaction_smiles, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.balance_check_json(ptr0, len0);
        deferred2_0 = ret[0];
        deferred2_1 = ret[1];
        return getStringFromWasm0(ret[0], ret[1]);
    } finally {
        wasm.__wbindgen_free(deferred2_0, deferred2_1, 1);
    }
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
 * Generate a self-contained HTML report for a newline-separated list of SMILES.
 *
 * Empty lines and invalid SMILES are silently skipped.
 * Returns the same card-grid HTML as Python's `chematic.report()`.
 * The input is limited to 1 MiB, 1,024 non-empty records, and 10,000 atoms
 * per successfully parsed molecule. Limit violations return an HTML error.
 *
 * ```js
 * const html = mod.batch_report_html("CCO\nc1ccccc1\nCC(=O)O");
 * const blob = new Blob([html], {type:'text/html'});
 * const url  = URL.createObjectURL(blob);
 * ```
 * @param {string} smiles_lines
 * @returns {string}
 */
export function batch_report_html(smiles_lines) {
    let deferred2_0;
    let deferred2_1;
    try {
        const ptr0 = passStringToWasm0(smiles_lines, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.batch_report_html(ptr0, len0);
        deferred2_0 = ret[0];
        deferred2_1 = ret[1];
        return getStringFromWasm0(ret[0], ret[1]);
    } finally {
        wasm.__wbindgen_free(deferred2_0, deferred2_1, 1);
    }
}

/**
 * Predict GI absorption and BBB penetration using the BOILED-Egg method
 * (Daina & Zoete 2016).
 *
 * Returns JSON: `{"gi_absorbed":bool,"bbb_penetrant":bool,"logp":f64,"tpsa":f64}`
 * @param {string} smiles
 * @returns {string}
 */
export function boiled_egg_json(smiles) {
    let deferred2_0;
    let deferred2_1;
    try {
        const ptr0 = passStringToWasm0(smiles, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.boiled_egg_json(ptr0, len0);
        deferred2_0 = ret[0];
        deferred2_1 = ret[1];
        return getStringFromWasm0(ret[0], ret[1]);
    } finally {
        wasm.__wbindgen_free(deferred2_0, deferred2_1, 1);
    }
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
 * BRICS fragment SMILES as a JSON array.
 *
 * Applies the BRICS fragmentation rules and returns the canonical SMILES of
 * every resulting fragment.  Returns `[]` for molecules with no BRICS-breakable
 * bonds (e.g. benzene).
 *
 * The count of fragments equals `brics_fragment_count`.
 * @param {MolHandle} mol
 * @returns {string}
 */
export function brics_fragments_json(mol) {
    let deferred1_0;
    let deferred1_1;
    try {
        _assertClass(mol, MolHandle);
        const ret = wasm.brics_fragments_json(mol.__wbg_ptr);
        deferred1_0 = ret[0];
        deferred1_1 = ret[1];
        return getStringFromWasm0(ret[0], ret[1]);
    } finally {
        wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
    }
}

/**
 * Cluster molecules by structural similarity (Butina algorithm, ECFP4 Tanimoto).
 *
 * `smiles_json` — a JSON array of SMILES strings.
 * `cutoff` — Tanimoto similarity threshold (0.0–1.0); molecules within this
 *   distance of a cluster centre are assigned to that cluster.
 * Returns a JSON array of clusters, each cluster being an array of 0-based input indices.
 * Returns a JS error if any SMILES fails to parse.
 * @param {string} smiles_json
 * @param {number} cutoff
 * @returns {string}
 */
export function butina_cluster_ecfp4_json(smiles_json, cutoff) {
    let deferred3_0;
    let deferred3_1;
    try {
        const ptr0 = passStringToWasm0(smiles_json, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.butina_cluster_ecfp4_json(ptr0, len0, cutoff);
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
 * Canonical tautomer of `mol`.
 *
 * Applies a rule-based tautomer normalisation and returns the canonical form
 * as a new `MolHandle`.
 * @param {MolHandle} mol
 * @returns {MolHandle}
 */
export function canonical_tautomer(mol) {
    _assertClass(mol, MolHandle);
    const ret = wasm.canonical_tautomer(mol.__wbg_ptr);
    return MolHandle.__wrap(ret);
}

/**
 * Compute the canonical tautomer with specific atoms blocked from H-transfer.
 *
 * `blocked_atom_indices_json`: JSON array of 0-based atom indices, e.g. `[0, 3]`.
 * Any tautomer move whose donor, bridge, or acceptor is in the blocked set is suppressed.
 *
 * Returns canonical SMILES of the result, or `{"error":"..."}` on failure.
 * Out-of-range indices are silently ignored (no effect).
 * @param {MolHandle} mol
 * @param {string} blocked_atom_indices_json
 * @returns {string}
 */
export function canonical_tautomer_with_blocked_atoms_json(mol, blocked_atom_indices_json) {
    let deferred2_0;
    let deferred2_1;
    try {
        _assertClass(mol, MolHandle);
        const ptr0 = passStringToWasm0(blocked_atom_indices_json, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.canonical_tautomer_with_blocked_atoms_json(mol.__wbg_ptr, ptr0, len0);
        deferred2_0 = ret[0];
        deferred2_1 = ret[1];
        return getStringFromWasm0(ret[0], ret[1]);
    } finally {
        wasm.__wbindgen_free(deferred2_0, deferred2_1, 1);
    }
}

/**
 * Parse all molecular fragments from a CDXML string.
 *
 * Returns a JSON array of SMILES strings, one per fragment:
 * `["CC","c1ccccc1"]`
 *
 * Stereochemistry (wedge/dash bonds) is read from the `Display` attribute
 * of bond elements.
 * @param {string} cdxml
 * @returns {string}
 */
export function cdxml_to_smiles_json(cdxml) {
    let deferred3_0;
    let deferred3_1;
    try {
        const ptr0 = passStringToWasm0(cdxml, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.cdxml_to_smiles_json(ptr0, len0);
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
 * Compute the charge Parent and return a status-shaped JSON result.
 * @param {MolHandle} mol
 * @returns {string}
 */
export function charge_parent_json(mol) {
    let deferred1_0;
    let deferred1_1;
    try {
        _assertClass(mol, MolHandle);
        const ret = wasm.charge_parent_json(mol.__wbg_ptr);
        deferred1_0 = ret[0];
        deferred1_1 = ret[1];
        return getStringFromWasm0(ret[0], ret[1]);
    } finally {
        wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
    }
}

/**
 * The `chematic-wasm` crate version (matches the workspace release version).
 *
 * Lets callers (e.g. the browser playground demo) display the running
 * version without hardcoding it — `demo/index.html` previously had a
 * static version string that silently went stale across releases.
 * @returns {string}
 */
export function chematic_version() {
    let deferred1_0;
    let deferred1_1;
    try {
        const ret = wasm.chematic_version();
        deferred1_0 = ret[0];
        deferred1_1 = ret[1];
        return getStringFromWasm0(ret[0], ret[1]);
    } finally {
        wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
    }
}

/**
 * CIP stereo assignments via the accurate hierarchical-digraph engine, as a JSON
 * array of `{atomIdx, cipCode}` objects -- same shape as [`cip_assignments_json`],
 * but merges the accurate engine's tetrahedral R/S (~99.6% oracle-stable agreement,
 * see `docs/rfcs/cip_accurate_rfc.md`) with legacy's E/Z and allene answers (the accurate
 * engine computes neither). Atoms it can't resolve are omitted here -- see
 * [`cip_unresolved_json`] -- never a silently-guessed label. Returns `"null"` on an
 * internal engine error (budget-independent computations should not normally hit this).
 * @param {MolHandle} mol
 * @returns {string}
 */
export function cip_assignments_accurate_json(mol) {
    let deferred1_0;
    let deferred1_1;
    try {
        _assertClass(mol, MolHandle);
        const ret = wasm.cip_assignments_accurate_json(mol.__wbg_ptr);
        deferred1_0 = ret[0];
        deferred1_1 = ret[1];
        return getStringFromWasm0(ret[0], ret[1]);
    } finally {
        wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
    }
}

/**
 * CIP stereo assignments as a JSON array of `{atomIdx, cipCode}` objects.
 *
 * `cipCode` is one of `"R"`, `"S"`, `"E"`, or `"Z"`.
 * Returns `[]` for molecules with no specified stereocenters.
 * @param {MolHandle} mol
 * @returns {string}
 */
export function cip_assignments_json(mol) {
    let deferred1_0;
    let deferred1_1;
    try {
        _assertClass(mol, MolHandle);
        const ret = wasm.cip_assignments_json(mol.__wbg_ptr);
        deferred1_0 = ret[0];
        deferred1_1 = ret[1];
        return getStringFromWasm0(ret[0], ret[1]);
    } finally {
        wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
    }
}

/**
 * Atoms the accurate CIP engine could not resolve a tetrahedral R/S for, as a JSON
 * array of `{atomIdx, reason}` objects. `reason` is `"tied"` (a genuine CIP-rule tie,
 * not a missing rule) or `"budgetExceeded"`. Always `[]` for the legacy engine (see
 * [`cip_assignments_json`]) -- it never reports "I don't know". Returns `"null"` on
 * an internal engine error.
 * @param {MolHandle} mol
 * @returns {string}
 */
export function cip_unresolved_json(mol) {
    let deferred1_0;
    let deferred1_1;
    try {
        _assertClass(mol, MolHandle);
        const ret = wasm.cip_unresolved_json(mol.__wbg_ptr);
        deferred1_0 = ret[0];
        deferred1_1 = ret[1];
        return getStringFromWasm0(ret[0], ret[1]);
    } finally {
        wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
    }
}

/**
 * Compare multiple SMILES strings (up to 256 by default).
 * Accepts a delimiter-separated list (e.g., newline or comma).
 * The input is limited to 1 MiB, 1,024 records, and 10,000 atoms per parsed
 * molecule.
 *
 * # Example (JS)
 * ```javascript
 * const smilesList = "c1ccccc1\nCc1ccccc1\nCCc1ccccc1";
 * const json = module.compare_molecules_batch_json(smilesList, "\n");
 * const comparison = JSON.parse(json);
 * ```
 * @param {string} smiles_batch
 * @param {string} delimiter
 * @returns {string}
 */
export function compare_molecules_batch_json(smiles_batch, delimiter) {
    let deferred4_0;
    let deferred4_1;
    try {
        const ptr0 = passStringToWasm0(smiles_batch, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ptr1 = passStringToWasm0(delimiter, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len1 = WASM_VECTOR_LEN;
        const ret = wasm.compare_molecules_batch_json(ptr0, len0, ptr1, len1);
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
 * Compare two or more SMILES strings (JSON string output).
 * Returns the JSON representation of a `MoleculeComparison` struct.
 *
 * # Example (JS)
 * ```javascript
 * const json = module.compare_molecules_json("c1ccccc1", "Cc1ccccc1");
 * const comparison = JSON.parse(json);
 * console.log(comparison.pairwise[0].similarities.ecfp4_tanimoto);
 * ```
 * @param {string} smiles1
 * @param {string} smiles2
 * @returns {string}
 */
export function compare_molecules_json(smiles1, smiles2) {
    let deferred4_0;
    let deferred4_1;
    try {
        const ptr0 = passStringToWasm0(smiles1, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ptr1 = passStringToWasm0(smiles2, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len1 = WASM_VECTOR_LEN;
        const ret = wasm.compare_molecules_json(ptr0, len0, ptr1, len1);
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
 * Generate multiple conformers with RMSD-based pruning.
 * Returns JSON: `{"conformers": [[[x,y,z],...], ...], "count": int}`.
 * @param {MolHandle} mol
 * @param {number} n
 * @param {number} rmsd_threshold
 * @returns {string}
 */
export function conformer_ensemble_json(mol, n, rmsd_threshold) {
    let deferred1_0;
    let deferred1_1;
    try {
        _assertClass(mol, MolHandle);
        const ret = wasm.conformer_ensemble_json(mol.__wbg_ptr, n, rmsd_threshold);
        deferred1_0 = ret[0];
        deferred1_1 = ret[1];
        return getStringFromWasm0(ret[0], ret[1]);
    } finally {
        wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
    }
}

/**
 * Convert between common topology-bearing molecular formats.
 *
 * This is the WASM counterpart of Python's `convert_format`. It supports
 * SMILES, MOL/SDF, MOL V3000, MOL2, CML, ChemicalJSON, MolJSON, and CDXML.
 * Coordinates and format-specific metadata are intentionally not carried
 * across this topology-only bridge; use the format-specific coordinate APIs
 * when those fields must be preserved.
 * @param {string} text
 * @param {string} input_format
 * @param {string} output_format
 * @returns {string}
 */
export function convert_common_format(text, input_format, output_format) {
    let deferred5_0;
    let deferred5_1;
    try {
        const ptr0 = passStringToWasm0(text, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ptr1 = passStringToWasm0(input_format, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len1 = WASM_VECTOR_LEN;
        const ptr2 = passStringToWasm0(output_format, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len2 = WASM_VECTOR_LEN;
        const ret = wasm.convert_common_format(ptr0, len0, ptr1, len1, ptr2, len2);
        var ptr4 = ret[0];
        var len4 = ret[1];
        if (ret[3]) {
            ptr4 = 0; len4 = 0;
            throw takeFromExternrefTable0(ret[2]);
        }
        deferred5_0 = ptr4;
        deferred5_1 = len4;
        return getStringFromWasm0(ptr4, len4);
    } finally {
        wasm.__wbindgen_free(deferred5_0, deferred5_1, 1);
    }
}

/**
 * Return the CPK color (CSS hex string) for the given element symbol.
 *
 * Returns `"#000000"` (black) for carbon and unknown elements.
 * @param {string} element_symbol
 * @returns {string}
 */
export function cpk_color(element_symbol) {
    let deferred2_0;
    let deferred2_1;
    try {
        const ptr0 = passStringToWasm0(element_symbol, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.cpk_color(ptr0, len0);
        deferred2_0 = ret[0];
        deferred2_1 = ret[1];
        return getStringFromWasm0(ret[0], ret[1]);
    } finally {
        wasm.__wbindgen_free(deferred2_0, deferred2_1, 1);
    }
}

/**
 * Parse a Gaussian Cube file and return its full [`chematic_mol::VolumetricGrid`]
 * as JSON: `{"origin":[x,y,z],"axes":[[..],[..],[..]],"shape":[nx,ny,nz],
 * "values":[...flat, row-major third-axis-fastest...],
 * "atoms":[{"element":"C","charge":6.0,"position":[x,y,z]}],
 * "units":"bohr"|"angstrom"}`. See module docs for the perf tradeoff of a
 * full `values` JSON round trip on a large grid.
 * @param {string} text
 * @returns {string}
 */
export function cube_grid_json(text) {
    let deferred3_0;
    let deferred3_1;
    try {
        const ptr0 = passStringToWasm0(text, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.cube_grid_json(ptr0, len0);
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
 * `[nx, ny, nz]` for a Gaussian Cube file's grid, as a `Uint32Array`.
 * @param {string} text
 * @returns {Uint32Array}
 */
export function cube_shape_u32(text) {
    const ptr0 = passStringToWasm0(text, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
    const len0 = WASM_VECTOR_LEN;
    const ret = wasm.cube_shape_u32(ptr0, len0);
    if (ret[2]) {
        throw takeFromExternrefTable0(ret[1]);
    }
    return takeFromExternrefTable0(ret[0]);
}

/**
 * Flat `values` from a Gaussian Cube file's grid, as a `Float64Array` --
 * same data [`cube_grid_json`]'s `"values"` field carries (row-major,
 * third-axis-fastest order -- see `chematic_mol::volumetric`'s module
 * docs for the exact index formula), as a real typed array instead of a
 * JSON number array.
 * @param {string} text
 * @returns {Float64Array}
 */
export function cube_values_f64(text) {
    const ptr0 = passStringToWasm0(text, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
    const len0 = WASM_VECTOR_LEN;
    const ret = wasm.cube_values_f64(ptr0, len0);
    if (ret[2]) {
        throw takeFromExternrefTable0(ret[1]);
    }
    return takeFromExternrefTable0(ret[0]);
}

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
 * @param {MolHandle} mol
 * @returns {string}
 */
export function depict_data_json(mol) {
    let deferred1_0;
    let deferred1_1;
    try {
        _assertClass(mol, MolHandle);
        const ret = wasm.depict_data_json(mol.__wbg_ptr);
        deferred1_0 = ret[0];
        deferred1_1 = ret[1];
        return getStringFromWasm0(ret[0], ret[1]);
    } finally {
        wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
    }
}

/**
 * Compute structured depiction data using caller-supplied 2D coordinates.
 *
 * `coords_json` — JSON array of `[x, y]` pairs, one per atom in order.
 *
 * Returns the same JSON format as `depict_data_json`.
 * @param {MolHandle} mol
 * @param {string} coords_json
 * @returns {string}
 */
export function depict_data_with_coords_json(mol, coords_json) {
    let deferred2_0;
    let deferred2_1;
    try {
        _assertClass(mol, MolHandle);
        const ptr0 = passStringToWasm0(coords_json, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.depict_data_with_coords_json(mol.__wbg_ptr, ptr0, len0);
        deferred2_0 = ret[0];
        deferred2_1 = ret[1];
        return getStringFromWasm0(ret[0], ret[1]);
    } finally {
        wasm.__wbindgen_free(deferred2_0, deferred2_1, 1);
    }
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
 * The input is limited to 1 MiB, 1,024 non-empty records, and 10,000 atoms
 * per successfully parsed molecule. Limit violations return an empty SVG
 * containing an error title.
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
 * Render a molecule grid with SMARTS-based atom highlighting.
 *
 * `smiles_block` — newline-separated SMILES strings (same format as `depict_svg_grid`).
 * `cols` — number of grid columns.
 * `match_smarts` — SMARTS pattern; matched atoms in each molecule are highlighted.
 *   Pass an empty string `""` to render without any highlighting.
 *
 * Invalid SMILES are rendered as empty cells; SMARTS parse failure returns an
 * unhighlighted grid (the SMARTS is silently ignored).
 * The input is limited to 1 MiB, 1,024 non-empty records, and 10,000 atoms
 * per successfully parsed molecule. Limit violations return an empty SVG
 * containing an error title.
 * @param {string} smiles_block
 * @param {number} cols
 * @param {string} match_smarts
 * @returns {string}
 */
export function depict_svg_grid_highlighted(smiles_block, cols, match_smarts) {
    let deferred3_0;
    let deferred3_1;
    try {
        const ptr0 = passStringToWasm0(smiles_block, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ptr1 = passStringToWasm0(match_smarts, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len1 = WASM_VECTOR_LEN;
        const ret = wasm.depict_svg_grid_highlighted(ptr0, len0, cols, ptr1, len1);
        deferred3_0 = ret[0];
        deferred3_1 = ret[1];
        return getStringFromWasm0(ret[0], ret[1]);
    } finally {
        wasm.__wbindgen_free(deferred3_0, deferred3_1, 1);
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
 * Dice similarity between `a` and `b` using ECFP4 fingerprints.
 * @param {MolHandle} a
 * @param {MolHandle} b
 * @returns {number}
 */
export function dice_ecfp4(a, b) {
    _assertClass(a, MolHandle);
    _assertClass(b, MolHandle);
    const ret = wasm.dice_ecfp4(a.__wbg_ptr, b.__wbg_ptr);
    return ret;
}

/**
 * Dice similarity between `a` and `b` using ECFP6 fingerprints.
 * @param {MolHandle} a
 * @param {MolHandle} b
 * @returns {number}
 */
export function dice_ecfp6(a, b) {
    _assertClass(a, MolHandle);
    _assertClass(b, MolHandle);
    const ret = wasm.dice_ecfp6(a.__wbg_ptr, b.__wbg_ptr);
    return ret;
}

/**
 * Dice similarity between `a` and `b` using MACCS 166-bit fingerprints.
 * @param {MolHandle} a
 * @param {MolHandle} b
 * @returns {number}
 */
export function dice_maccs(a, b) {
    _assertClass(a, MolHandle);
    _assertClass(b, MolHandle);
    const ret = wasm.dice_maccs(a.__wbg_ptr, b.__wbg_ptr);
    return ret;
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
 * Like `ecfp4_bitvec` but with explicit chirality control.
 *
 * When `use_chirality=true`, tetrahedral stereochemistry is included in the
 * initial atom hash, making enantiomers have different fingerprints.
 * When `false` (default), chirality is ignored.
 * @param {MolHandle} mol
 * @param {boolean} use_chirality
 * @returns {Uint8Array}
 */
export function ecfp4_bitvec_with_chirality(mol, use_chirality) {
    _assertClass(mol, MolHandle);
    const ret = wasm.ecfp4_bitvec_with_chirality(mol.__wbg_ptr, use_chirality);
    var v1 = getArrayU8FromWasm0(ret[0], ret[1]).slice();
    wasm.__wbindgen_free(ret[0], ret[1] * 1, 1);
    return v1;
}

/**
 * ECFP6 (radius-3) fingerprint as a bit-packed byte vector (256 bytes = 2048 bits).
 * @param {MolHandle} mol
 * @returns {Uint8Array}
 */
export function ecfp6_bitvec(mol) {
    _assertClass(mol, MolHandle);
    const ret = wasm.ecfp6_bitvec(mol.__wbg_ptr);
    var v1 = getArrayU8FromWasm0(ret[0], ret[1]).slice();
    wasm.__wbindgen_free(ret[0], ret[1] * 1, 1);
    return v1;
}

/**
 * Like `ecfp6_bitvec` but with explicit chirality control.
 *
 * When `use_chirality=true`, tetrahedral stereochemistry is included in the
 * initial atom hash, making enantiomers have different fingerprints.
 * When `false` (default), chirality is ignored.
 * @param {MolHandle} mol
 * @param {boolean} use_chirality
 * @returns {Uint8Array}
 */
export function ecfp6_bitvec_with_chirality(mol, use_chirality) {
    _assertClass(mol, MolHandle);
    const ret = wasm.ecfp6_bitvec_with_chirality(mol.__wbg_ptr, use_chirality);
    var v1 = getArrayU8FromWasm0(ret[0], ret[1]).slice();
    wasm.__wbindgen_free(ret[0], ret[1] * 1, 1);
    return v1;
}

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
 * @param {MolHandle} mol
 * @param {number} radius
 * @param {number} nbits
 * @param {boolean} use_chirality
 * @returns {Uint8Array}
 */
export function ecfp_bitvec_custom(mol, radius, nbits, use_chirality) {
    _assertClass(mol, MolHandle);
    const ret = wasm.ecfp_bitvec_custom(mol.__wbg_ptr, radius, nbits, use_chirality);
    var v1 = getArrayU8FromWasm0(ret[0], ret[1]).slice();
    wasm.__wbindgen_free(ret[0], ret[1] * 1, 1);
    return v1;
}

/**
 * Run `embed_ensemble_v2` on `mol`'s own atom order (never canonicalizes/
 * reparses, same convention as `embed_pipeline_v2_json`). See the module doc
 * for the config JSON shape and the success/failure envelope shape.
 *
 * Never throws. Always returns a JSON string tagged `schemaVersion: 1` and
 * `ok: true`/`false`. `ok: false` covers both a WASM-level rejection
 * (oversized input, too many atoms, malformed/incomplete config JSON) and a
 * config `embed_ensemble_v2` itself rejects (currently only an invalid
 * `rmsdThreshold`) -- distinguished by `error.stage` /
 * `error.cause.kind`, same as `embed_pipeline_v2_json`.
 * @param {MolHandle} mol
 * @param {string} config_json
 * @returns {string}
 */
export function embed_ensemble_v2_json(mol, config_json) {
    let deferred2_0;
    let deferred2_1;
    try {
        _assertClass(mol, MolHandle);
        const ptr0 = passStringToWasm0(config_json, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.embed_ensemble_v2_json(mol.__wbg_ptr, ptr0, len0);
        deferred2_0 = ret[0];
        deferred2_1 = ret[1];
        return getStringFromWasm0(ret[0], ret[1]);
    } finally {
        wasm.__wbindgen_free(deferred2_0, deferred2_1, 1);
    }
}

/**
 * Run the opt-in v2 embedding pipeline, applied directly to `mol`'s own atom
 * order (never canonicalizes/reparses -- see the module doc).
 *
 * `config_json` must be an object with the 15 required fields
 * `PipelineV2Config` requires (camelCase keys: `embedSeed`, `maxAttempts`,
 * `embedTimeoutMs`, `useExpTorsions`, `useSmallRingTorsions`,
 * `useMacrocycleTorsions`, `useMacrocycle14Bounds`,
 * `includeLegacyTorsionHeuristic`, `stereoPolicy`, `failOnUnevaluableStereo`,
 * `forceFieldPolicy`, `forceFieldMaxIterations`, `gateMmff94TorsionOop`,
 * `ringTorsionPolicy`, `totalTimeoutMs`), plus one optional field added in
 * Priority 2 (issue #227): `gateMmff94StretchBend` (`#[serde(default)]` ->
 * `false` if omitted, so pre-Priority-2 caller configs keep working
 * unmodified, matching `false`'s meaning of "existing/unchanged behavior"
 * everywhere else in this codebase). An unknown field, a missing *required*
 * field, an unknown `stereoPolicy`/`ringTorsionPolicy`/`forceFieldPolicy`
 * string, or a wrong-typed/out-of-range integer all fail closed rather than
 * silently defaulting.
 *
 * Never throws. Always returns a JSON string tagged with `schemaVersion: 1` and
 * `ok: true`/`false` -- see the module doc for both shapes.
 * @param {MolHandle} mol
 * @param {string} config_json
 * @returns {string}
 */
export function embed_pipeline_v2_json(mol, config_json) {
    let deferred2_0;
    let deferred2_1;
    try {
        _assertClass(mol, MolHandle);
        const ptr0 = passStringToWasm0(config_json, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.embed_pipeline_v2_json(mol.__wbg_ptr, ptr0, len0);
        deferred2_0 = ret[0];
        deferred2_1 = ret[1];
        return getStringFromWasm0(ret[0], ret[1]);
    } finally {
        wasm.__wbindgen_free(deferred2_0, deferred2_1, 1);
    }
}

/**
 * Enumerate a combinatorial library from a SMIRKS template and two fragment sets.
 *
 * Generates all products by combining every scaffold with every building block.
 * Input format: `scaffolds_smiles` and `building_blocks_smiles` are pipe-delimited
 * SMILES strings (e.g., `"c1ccccc1|Cc1ccccc1"`).
 *
 * Returns JSON array of product SMILES strings.
 * Example: `enumerate_library_2way("[C:1][Cl].[C:2][NH2]>>[C:1]N[C:2]", "c1ccccc1|Cc1ccccc1", "NCc1ccccc1|NCC")`
 * @param {string} template
 * @param {string} scaffolds_smiles
 * @param {string} building_blocks_smiles
 * @returns {string}
 */
export function enumerate_library_2way(template, scaffolds_smiles, building_blocks_smiles) {
    let deferred5_0;
    let deferred5_1;
    try {
        const ptr0 = passStringToWasm0(template, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ptr1 = passStringToWasm0(scaffolds_smiles, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len1 = WASM_VECTOR_LEN;
        const ptr2 = passStringToWasm0(building_blocks_smiles, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len2 = WASM_VECTOR_LEN;
        const ret = wasm.enumerate_library_2way(ptr0, len0, ptr1, len1, ptr2, len2);
        var ptr4 = ret[0];
        var len4 = ret[1];
        if (ret[3]) {
            ptr4 = 0; len4 = 0;
            throw takeFromExternrefTable0(ret[2]);
        }
        deferred5_0 = ptr4;
        deferred5_1 = len4;
        return getStringFromWasm0(ptr4, len4);
    } finally {
        wasm.__wbindgen_free(deferred5_0, deferred5_1, 1);
    }
}

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
 * @param {MolHandle} mol
 * @returns {string}
 */
export function enumerate_stereo_isomers_json(mol) {
    let deferred2_0;
    let deferred2_1;
    try {
        _assertClass(mol, MolHandle);
        const ret = wasm.enumerate_stereo_isomers_json(mol.__wbg_ptr);
        var ptr1 = ret[0];
        var len1 = ret[1];
        if (ret[3]) {
            ptr1 = 0; len1 = 0;
            throw takeFromExternrefTable0(ret[2]);
        }
        deferred2_0 = ptr1;
        deferred2_1 = len1;
        return getStringFromWasm0(ptr1, len1);
    } finally {
        wasm.__wbindgen_free(deferred2_0, deferred2_1, 1);
    }
}

/**
 * All enumerated tautomers of `mol` as a JSON array of canonical SMILES strings.
 *
 * Example return value: `["Oc1cccc2ccccc12","O=C1C=CC=Cc2ccccc21"]`
 * @param {MolHandle} mol
 * @returns {string}
 */
export function enumerate_tautomers_json(mol) {
    let deferred1_0;
    let deferred1_1;
    try {
        _assertClass(mol, MolHandle);
        const ret = wasm.enumerate_tautomers_json(mol.__wbg_ptr);
        deferred1_0 = ret[0];
        deferred1_1 = ret[1];
        return getStringFromWasm0(ret[0], ret[1]);
    } finally {
        wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
    }
}

/**
 * Compute ERG-style 315-element float histogram fingerprint.
 * Returns JSON: {"len":315,"values":[f64,...]} or {"error":"..."}.
 * Format: 21 pharmacophore-feature-pair × 15 distance bins with Gaussian fuzzing.
 * See `chematic_fp::erg_vec` for details.
 * @param {MolHandle} mol
 * @returns {string}
 */
export function erg_vec_json(mol) {
    let deferred1_0;
    let deferred1_1;
    try {
        _assertClass(mol, MolHandle);
        const ret = wasm.erg_vec_json(mol.__wbg_ptr);
        deferred1_0 = ret[0];
        deferred1_1 = ret[1];
        return getStringFromWasm0(ret[0], ret[1]);
    } finally {
        wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
    }
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
 * Extract an extxyz frame's coordinates, cell, per-atom properties and
 * frame metadata as JSON, in the SAME atom order [`mol_from_extxyz`]
 * returns topology for (both read the identical frame via
 * `chematic_mol::parse_extxyz`, so atom-index correspondence between the
 * two calls is structural, not just conventional).
 *
 * Returns JSON `{"coords":[[x,y,z],...],
 * "lattice":[9 numbers]|null,
 * "properties":{"name":[[...atom values...], ...], ...},
 * "info":{"key":"value", ...}}`.
 *
 * Returns a JS error on parse failure or if the frame exceeds the WASM
 * atom-count limit.
 * @param {string} text
 * @returns {string}
 */
export function extxyz_frame_json(text) {
    let deferred3_0;
    let deferred3_1;
    try {
        const ptr0 = passStringToWasm0(text, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.extxyz_frame_json(ptr0, len0);
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
 * FCFP4 (pharmacophore, radius-2) fingerprint as a bit-packed byte vector (256 bytes).
 * @param {MolHandle} mol
 * @returns {Uint8Array}
 */
export function fcfp4_bitvec(mol) {
    _assertClass(mol, MolHandle);
    const ret = wasm.fcfp4_bitvec(mol.__wbg_ptr);
    var v1 = getArrayU8FromWasm0(ret[0], ret[1]).slice();
    wasm.__wbindgen_free(ret[0], ret[1] * 1, 1);
    return v1;
}

/**
 * FCFP6 (pharmacophore, radius-3) fingerprint as a bit-packed byte vector (256 bytes).
 * @param {MolHandle} mol
 * @returns {Uint8Array}
 */
export function fcfp6_bitvec(mol) {
    _assertClass(mol, MolHandle);
    const ret = wasm.fcfp6_bitvec(mol.__wbg_ptr);
    var v1 = getArrayU8FromWasm0(ret[0], ret[1]).slice();
    wasm.__wbindgen_free(ret[0], ret[1] * 1, 1);
    return v1;
}

/**
 * Analyze a reaction SMILES and return the reaction center as JSON.
 *
 * JSON schema: `{ broken: [[a1,a2],...], formed: [[a1,a2],...], changed: [a,...] }`
 * where atom indices are 0-based within the first reactant molecule.
 * Returns an error string prefixed with `"error:"` on failure.
 * @param {string} reaction_smiles
 * @returns {string}
 */
export function find_reaction_center_json(reaction_smiles) {
    let deferred2_0;
    let deferred2_1;
    try {
        const ptr0 = passStringToWasm0(reaction_smiles, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.find_reaction_center_json(ptr0, len0);
        deferred2_0 = ret[0];
        deferred2_1 = ret[1];
        return getStringFromWasm0(ret[0], ret[1]);
    } finally {
        wasm.__wbindgen_free(deferred2_0, deferred2_1, 1);
    }
}

/**
 * Compute the fragment Parent and return a status-shaped JSON result.
 * @param {MolHandle} mol
 * @returns {string}
 */
export function fragment_parent_json(mol) {
    let deferred1_0;
    let deferred1_1;
    try {
        _assertClass(mol, MolHandle);
        const ret = wasm.fragment_parent_json(mol.__wbg_ptr);
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
 * Generate 3D coordinates as raw JSON array [[x,y,z], ...].
 *
 * Unlike `generate_3d_pdb`, this returns coordinates that can be passed
 * to descriptor functions like `whim_descriptors_json` or `shape_descriptors_json`.
 * @param {MolHandle} mol
 * @returns {string}
 */
export function generate_3d_coords_json(mol) {
    let deferred1_0;
    let deferred1_1;
    try {
        _assertClass(mol, MolHandle);
        const ret = wasm.generate_3d_coords_json(mol.__wbg_ptr);
        deferred1_0 = ret[0];
        deferred1_1 = ret[1];
        return getStringFromWasm0(ret[0], ret[1]);
    } finally {
        wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
    }
}

/**
 * Generate 3D coordinates using ETKDG as raw JSON array [[x,y,z], ...].
 * @param {MolHandle} mol
 * @returns {string}
 */
export function generate_3d_etkdg_coords_json(mol) {
    let deferred1_0;
    let deferred1_1;
    try {
        _assertClass(mol, MolHandle);
        const ret = wasm.generate_3d_etkdg_coords_json(mol.__wbg_ptr);
        deferred1_0 = ret[0];
        deferred1_1 = ret[1];
        return getStringFromWasm0(ret[0], ret[1]);
    } finally {
        wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
    }
}

/**
 * Generate 3D coordinates using ETKDG and minimize with DREIDING force field.
 * @param {MolHandle} mol
 * @returns {string}
 */
export function generate_3d_etkdg_minimized_pdb(mol) {
    let deferred1_0;
    let deferred1_1;
    try {
        _assertClass(mol, MolHandle);
        const ret = wasm.generate_3d_etkdg_minimized_pdb(mol.__wbg_ptr);
        deferred1_0 = ret[0];
        deferred1_1 = ret[1];
        return getStringFromWasm0(ret[0], ret[1]);
    } finally {
        wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
    }
}

/**
 * Generate 3D coordinates using ETKDG (torsion angle preferences) and return PDB block.
 * ETKDG produces higher-quality conformations than rule-based DG by applying
 * experimental torsion angle preferences to common structural patterns.
 * @param {MolHandle} mol
 * @returns {string}
 */
export function generate_3d_etkdg_pdb(mol) {
    let deferred1_0;
    let deferred1_1;
    try {
        _assertClass(mol, MolHandle);
        const ret = wasm.generate_3d_etkdg_pdb(mol.__wbg_ptr);
        deferred1_0 = ret[0];
        deferred1_1 = ret[1];
        return getStringFromWasm0(ret[0], ret[1]);
    } finally {
        wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
    }
}

/**
 * Generate 3D coordinates from SMILES (raw distance geometry, no minimization).
 * Returns PDB format string with atoms positioned in 3D space.
 *
 * # Example (JS)
 * ```javascript
 * const pdbStr = module.generate_3d_from_smiles("c1ccccc1");
 * console.log(pdbStr);  // PDB file content
 * ```
 * @param {string} smiles
 * @returns {string}
 */
export function generate_3d_from_smiles(smiles) {
    let deferred3_0;
    let deferred3_1;
    try {
        const ptr0 = passStringToWasm0(smiles, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.generate_3d_from_smiles(ptr0, len0);
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
 * Generate energy-minimized 3D coordinates and return a PDB string.
 *
 * Runs distance-geometry placement followed by gradient-descent force-field
 * minimization.  Geometry quality is better than `generate_3d_pdb` for
 * flexible molecules; the force field is approximate (not MMFF94/UFF).
 * @param {MolHandle} mol
 * @returns {string}
 */
export function generate_3d_minimized_pdb(mol) {
    let deferred1_0;
    let deferred1_1;
    try {
        _assertClass(mol, MolHandle);
        const ret = wasm.generate_3d_minimized_pdb(mol.__wbg_ptr);
        deferred1_0 = ret[0];
        deferred1_1 = ret[1];
        return getStringFromWasm0(ret[0], ret[1]);
    } finally {
        wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
    }
}

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
 * @param {string} smiles
 * @returns {string}
 */
export function generate_3d_optimized_pdb(smiles) {
    let deferred3_0;
    let deferred3_1;
    try {
        const ptr0 = passStringToWasm0(smiles, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.generate_3d_optimized_pdb(ptr0, len0);
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
 * Generic (atom-type-erased) Murcko scaffold of `mol`.
 *
 * All atoms become carbon and all bonds become single bonds, giving the pure
 * graph topology of the scaffold.
 * @param {MolHandle} mol
 * @returns {MolHandle}
 */
export function generic_murcko_scaffold(mol) {
    _assertClass(mol, MolHandle);
    const ret = wasm.generic_murcko_scaffold(mol.__wbg_ptr);
    return MolHandle.__wrap(ret);
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
 * @param {string} smiles
 * @param {number} a
 * @param {number} b
 * @returns {number}
 */
export function get_bond_length_json(smiles, a, b) {
    const ptr0 = passStringToWasm0(smiles, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
    const len0 = WASM_VECTOR_LEN;
    const ret = wasm.get_bond_length_json(ptr0, len0, a, b);
    return ret;
}

/**
 * All scalar molecular descriptors as a single JSON object.
 *
 * Keys use camelCase and match the individual `MolHandle` method names.
 * Drug-likeness rule outcomes are included as boolean fields.
 * @param {MolHandle} mol
 * @returns {string}
 */
export function get_descriptors_json(mol) {
    let deferred1_0;
    let deferred1_1;
    try {
        _assertClass(mol, MolHandle);
        const ret = wasm.get_descriptors_json(mol.__wbg_ptr);
        deferred1_0 = ret[0];
        deferred1_1 = ret[1];
        return getStringFromWasm0(ret[0], ret[1]);
    } finally {
        wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
    }
}

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
 * @param {string} smiles
 * @param {number} a
 * @param {number} b
 * @param {number} c
 * @param {number} d
 * @returns {any}
 */
export function get_dihedral_json(smiles, a, b, c, d) {
    const ptr0 = passStringToWasm0(smiles, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
    const len0 = WASM_VECTOR_LEN;
    const ret = wasm.get_dihedral_json(ptr0, len0, a, b, c, d);
    return ret;
}

/**
 * Compute GETAWAY descriptors (GEometry, Topology and Atom-Weights AssemblY) from 3D coords.
 *
 * Returns a JSON array of **19** values:
 * - `[0..7]`  H[1..8]  — leverage autocorrelation at topological lags 1–8
 * - `[8..15]` R[1..8]  — H[k] normalised by pair count W_k
 * - `[16]` Hmax, `[17]` Hmean, `[18]` Htot — per-atom leverage statistics
 *
 * Note: requires 3D coordinates (non-planar); for flat/2D structures the hat matrix
 * is degenerate and descriptors reflect squared centroid distances, not true leverage.
 * @param {MolHandle} mol
 * @returns {string}
 */
export function getaway_descriptors_json(mol) {
    let deferred1_0;
    let deferred1_1;
    try {
        _assertClass(mol, MolHandle);
        const ret = wasm.getaway_descriptors_json(mol.__wbg_ptr);
        deferred1_0 = ret[0];
        deferred1_1 = ret[1];
        return getStringFromWasm0(ret[0], ret[1]);
    } finally {
        wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
    }
}

/**
 * Parse a ChemDraw XML (CDXML) string into a `MolHandle`.
 *
 * Compute an HDF fingerprint and return it as a JSON array of float32 values.
 *
 * Returns a unit-norm vector of length `dim` as a JSON number array.
 * Use cosine dot product for similarity: `a · b = sum(a[i]*b[i])`.
 *
 * ```js
 * const fp = JSON.parse(hdf_json(mol));        // float[] of length 1024
 * const sim = fp.reduce((s, v, i) => s + v * fp2[i], 0);  // cosine similarity
 * ```
 * @param {MolHandle} mol
 * @param {number} dim
 * @param {number} radius
 * @param {bigint} seed
 * @returns {string}
 */
export function hdf_json(mol, dim, radius, seed) {
    let deferred1_0;
    let deferred1_1;
    try {
        _assertClass(mol, MolHandle);
        const ret = wasm.hdf_json(mol.__wbg_ptr, dim, radius, seed);
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
 * Generate InChI string from SMILES.
 *
 * Returns `"error:<msg>"` on parse failure.
 * @param {string} smiles
 * @returns {string}
 */
export function inchi_from_smiles(smiles) {
    let deferred2_0;
    let deferred2_1;
    try {
        const ptr0 = passStringToWasm0(smiles, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.inchi_from_smiles(ptr0, len0);
        deferred2_0 = ret[0];
        deferred2_1 = ret[1];
        return getStringFromWasm0(ret[0], ret[1]);
    } finally {
        wasm.__wbindgen_free(deferred2_0, deferred2_1, 1);
    }
}

/**
 * Generate InChIKey from SMILES (27-character identifier).
 *
 * Returns `"error:<msg>"` on parse failure.
 * @param {string} smiles
 * @returns {string}
 */
export function inchikey_from_smiles(smiles) {
    let deferred2_0;
    let deferred2_1;
    try {
        const ptr0 = passStringToWasm0(smiles, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.inchikey_from_smiles(ptr0, len0);
        deferred2_0 = ret[0];
        deferred2_1 = ret[1];
        return getStringFromWasm0(ret[0], ret[1]);
    } finally {
        wasm.__wbindgen_free(deferred2_0, deferred2_1, 1);
    }
}

/**
 * Invert the stereochemistry of a tetrahedral stereocenter (U/D wedge bonds).
 *
 * If the atom has no wedge/dash bonds, returns an unchanged copy.
 * Returns error if atom_idx is invalid.
 * @param {MolHandle} mol
 * @param {number} atom_idx
 * @returns {MolHandle}
 */
export function invert_stereocenter_at(mol, atom_idx) {
    _assertClass(mol, MolHandle);
    const ret = wasm.invert_stereocenter_at(mol.__wbg_ptr, atom_idx);
    if (ret[2]) {
        throw takeFromExternrefTable0(ret[1]);
    }
    return MolHandle.__wrap(ret[0]);
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
 * Compute the isotope Parent and return a status-shaped JSON result.
 * @param {MolHandle} mol
 * @returns {string}
 */
export function isotope_parent_json(mol) {
    let deferred1_0;
    let deferred1_1;
    try {
        _assertClass(mol, MolHandle);
        const ret = wasm.isotope_parent_json(mol.__wbg_ptr);
        deferred1_0 = ret[0];
        deferred1_1 = ret[1];
        return getStringFromWasm0(ret[0], ret[1]);
    } finally {
        wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
    }
}

/**
 * Per-atom Labute approximate surface area contributions as a JSON array of f64.
 *
 * Non-finite values (single-atom molecules etc.) are emitted as JSON `null`.
 * @param {MolHandle} mol
 * @returns {string}
 */
export function labute_asa_per_atom_json(mol) {
    let deferred1_0;
    let deferred1_1;
    try {
        _assertClass(mol, MolHandle);
        const ret = wasm.labute_asa_per_atom_json(mol.__wbg_ptr);
        deferred1_0 = ret[0];
        deferred1_1 = ret[1];
        return getStringFromWasm0(ret[0], ret[1]);
    } finally {
        wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
    }
}

/**
 * Parse a LAMMPS data file (`read_data` format) and return every section
 * as JSON: `{"counts":[["atoms",120],["atom types",4],...],
 * "atom_style":"atomic"|"charge"|"molecular"|"full"|"<other>",
 * "simulation_box":{"lo":[x,y,z],"hi":[x,y,z],"tilt":[xy,xz,yz]|null},
 * "masses":[{"atom_type":N,"mass":N}],
 * "atoms":[{"id":N,"molecule_id":N|null,"atom_type":N,"charge":N|null,"x":N,"y":N,"z":N,"image":[ix,iy,iz]|null}],
 * "velocities":[{"atom_id":N,"vx":N,"vy":N,"vz":N}],
 * "bonds":[{"id":N,"bond_type":N,"atom1":N,"atom2":N}],
 * "unparsed_sections":[["Angles","<raw row text>"],...]}`. `atom_type`
 * must be exactly `"atomic"`/`"charge"`/`"molecular"`/`"full"` -- LAMMPS's
 * atom style is not recoverable from the file itself (see
 * [`chematic_mol::LammpsData`]'s module doc comment); any other value is
 * rejected with a JS error, matching
 * [`chematic_mol::LammpsDataError::UnsupportedAtomStyle`].
 *
 * This module has no bond-perception step of its own: `Angles`/
 * `Dihedrals`/`Impropers`/`*Coeffs`/any other section not listed above
 * are preserved verbatim (byte-for-byte, `#` comments included) in
 * `unparsed_sections`, not modeled field-by-field.
 * @param {string} text
 * @param {string} atom_style
 * @returns {string}
 */
export function lammps_data_to_json(text, atom_style) {
    let deferred4_0;
    let deferred4_1;
    try {
        const ptr0 = passStringToWasm0(text, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ptr1 = passStringToWasm0(atom_style, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len1 = WASM_VECTOR_LEN;
        const ret = wasm.lammps_data_to_json(ptr0, len0, ptr1, len1);
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
 * Like [`lammps_dump_cartesian_positions_json`], but returns a flat
 * `Float64Array` (`[x0,y0,z0,x1,y1,z1,...]`, 3 values per atom) instead
 * of a JSON `[[x,y,z],...]` array.
 *
 * **Behavioral difference from the JSON sibling**: when the frame has no
 * recognized coordinate columns, [`lammps_dump_cartesian_positions_json`]
 * returns JSON `null`; a `Float64Array` has no `null`, so this function
 * returns `Err` instead, with a message naming the columns it looked for.
 * @param {string} frame_json
 * @returns {Float64Array}
 */
export function lammps_dump_cartesian_positions_f64(frame_json) {
    const ptr0 = passStringToWasm0(frame_json, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
    const len0 = WASM_VECTOR_LEN;
    const ret = wasm.lammps_dump_cartesian_positions_f64(ptr0, len0);
    if (ret[2]) {
        throw takeFromExternrefTable0(ret[1]);
    }
    return takeFromExternrefTable0(ret[0]);
}

/**
 * Real Cartesian positions for a LAMMPS dump frame (in the JSON shape
 * [`lammps_dump_frame_to_json_str`] returns), resolved by delegating
 * directly to [`chematic_mol::LammpsDumpFrame::cartesian_positions`] --
 * this function does not reimplement any part of the box-bounds or
 * scaled-coordinate math itself; that method is the single place this
 * crate gets the (orthogonal or triclinic) transform right, and every
 * WASM caller must go through it rather than re-deriving the transform in
 * JS (the same reasoning behind this crate's OpenDX fail-closed unit
 * handling and QCSchema's single Bohr<->Ångström conversion point).
 *
 * - `x y z` columns: passed straight through.
 * - `xs ys zs` columns: transformed through `frame.box_bounds` (including
 *   the triclinic shear terms when a tilt is present).
 * - Neither present (including an `xu yu zu`-only frame -- "unwrapped" is
 *   a materially different physical quantity from a scaled coordinate,
 *   never resolved by this method): returns JSON `null`, not an error and
 *   not an empty array, matching
 *   [`chematic_mol::LammpsDumpFrame::cartesian_positions`]'s own
 *   `Option` semantics exactly.
 *
 * Returns JSON `[[x,y,z],...]` on success, in the same atom order as
 * `frame.rows`. See [`lammps_dump_cartesian_positions_f64`] for a flat
 * `Float64Array` sibling -- note its `null` case becomes an `Err` there
 * instead, a disclosed, real API-shape difference (a typed array has no
 * `null`).
 * @param {string} frame_json
 * @returns {string}
 */
export function lammps_dump_cartesian_positions_json(frame_json) {
    let deferred3_0;
    let deferred3_1;
    try {
        const ptr0 = passStringToWasm0(frame_json, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.lammps_dump_cartesian_positions_json(ptr0, len0);
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
 * Parse a single LAMMPS dump/trajectory frame and return it as JSON:
 * `{"timestep":N,"num_atoms":N,
 * "box_bounds":{"lo":[x,y,z],"hi":[x,y,z],"tilt":[xy,xz,yz]|null},
 * "boundary_flags":["pp","pp","pp"],"column_names":[...],
 * "rows":[[...values, one per column_names entry...],...]}`.
 * `box_bounds` is already the resolved TRUE simulation box (the parser
 * applies [`chematic_mol::box_bounds_to_true`] internally before
 * `LammpsDumpFrame` is ever built) -- not the file's raw
 * `xlo_bound`/`xhi_bound`/... values. `rows` is the raw per-atom column
 * data as declared by `column_names`, which may be `x y z`
 * (already-Cartesian), `xs ys zs` (box-scaled), `xu yu zu` (unwrapped), or
 * any other dump-command column -- use
 * [`lammps_dump_cartesian_positions_json`] to resolve real Cartesian
 * positions from whichever convention is present, rather than
 * reimplementing that resolution/transform in JS.
 * @param {string} text
 * @returns {string}
 */
export function lammps_dump_frame_to_json_str(text) {
    let deferred3_0;
    let deferred3_1;
    try {
        const ptr0 = passStringToWasm0(text, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.lammps_dump_frame_to_json_str(ptr0, len0);
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
 * Flattens a LAMMPS dump frame's `rows` (JSON shape
 * [`lammps_dump_frame_to_json_str`] returns) into a single flat
 * `Float64Array`, row-major (atom 0's `column_names.len()` values, then
 * atom 1's, ...). The caller already has `column_names` from
 * [`lammps_dump_frame_to_json_str`] and can compute the row length
 * itself (`column_names.length`); no separate row-length accessor is
 * provided here.
 * @param {string} frame_json
 * @returns {Float64Array}
 */
export function lammps_dump_rows_f64(frame_json) {
    const ptr0 = passStringToWasm0(frame_json, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
    const len0 = WASM_VECTOR_LEN;
    const ret = wasm.lammps_dump_rows_f64(ptr0, len0);
    if (ret[2]) {
        throw takeFromExternrefTable0(ret[1]);
    }
    return takeFromExternrefTable0(ret[0]);
}

/**
 * Parse every frame of a LAMMPS dump/trajectory file and return them as a
 * JSON array (same per-frame shape as [`lammps_dump_frame_to_json_str`]).
 *
 * This reads the whole input, parses it fully, and returns every frame at
 * once -- [`chematic_mol::LammpsDumpReader`]'s per-frame streaming
 * iteration (reading one frame at a time from a `BufRead` without holding
 * the whole trajectory in memory) has no natural equivalent across the
 * JS/WASM boundary in this first pass and is deliberately not exposed
 * here, not silently dropped: a JS caller with a truly large trajectory
 * that needs bounded memory should process it server-side instead.
 * @param {string} text
 * @returns {string}
 */
export function lammps_trajectory_to_json(text) {
    let deferred3_0;
    let deferred3_1;
    try {
        const ptr0 = passStringToWasm0(text, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.lammps_trajectory_to_json(ptr0, len0);
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
 * Return the largest fragment of `mol` (salt/solvent stripping).
 *
 * For single-component molecules returns a copy of the same molecule.
 * @param {MolHandle} mol
 * @returns {MolHandle}
 */
export function largest_fragment(mol) {
    _assertClass(mol, MolHandle);
    const ret = wasm.largest_fragment(mol.__wbg_ptr);
    return MolHandle.__wrap(ret);
}

/**
 * Per-atom Crippen LogP contributions as a JSON array of f64.
 *
 * Index `i` corresponds to atom `i` in `mol.atoms()` order.
 * @param {MolHandle} mol
 * @returns {string}
 */
export function logp_per_atom_json(mol) {
    let deferred1_0;
    let deferred1_1;
    try {
        _assertClass(mol, MolHandle);
        const ret = wasm.logp_per_atom_json(mol.__wbg_ptr);
        deferred1_0 = ret[0];
        deferred1_1 = ret[1];
        return getStringFromWasm0(ret[0], ret[1]);
    } finally {
        wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
    }
}

/**
 * MACCS 166-bit structural keys fingerprint as a byte array (21 bytes, LSB-first).
 *
 * Bit `i` (0-indexed) corresponds to MACCS key `i+1`.
 * @param {MolHandle} mol
 * @returns {Uint8Array}
 */
export function maccs_bitvec(mol) {
    _assertClass(mol, MolHandle);
    const ret = wasm.maccs_bitvec(mol.__wbg_ptr);
    var v1 = getArrayU8FromWasm0(ret[0], ret[1]).slice();
    wasm.__wbindgen_free(ret[0], ret[1] * 1, 1);
    return v1;
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
 * Select `n` maximally-diverse molecules (MaxMin algorithm, ECFP4 Tanimoto).
 *
 * `smiles_json` — a JSON array of SMILES strings, e.g. `["CC","c1ccccc1","CCO"]`.
 * Returns a JSON array of 0-based indices into the input array.
 * Returns a JS error if any SMILES fails to parse (indices would otherwise shift).
 * @param {string} smiles_json
 * @param {number} n
 * @returns {string}
 */
export function maxmin_picks_ecfp4_json(smiles_json, n) {
    let deferred3_0;
    let deferred3_1;
    try {
        const ptr0 = passStringToWasm0(smiles_json, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.maxmin_picks_ecfp4_json(ptr0, len0, n);
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
 * Maximum Common Substructure of a set of molecules, returned as a canonical SMILES string.
 *
 * `smiles_json` — a JSON array of at least 2 SMILES strings.
 * Returns the MCS SMILES, or `"null"` when no common substructure was found.
 * Returns a JS error on SMILES parse failure.
 * @param {string} smiles_json
 * @returns {string}
 */
export function mcs_smiles_json(smiles_json) {
    let deferred3_0;
    let deferred3_1;
    try {
        const ptr0 = passStringToWasm0(smiles_json, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.mcs_smiles_json(ptr0, len0);
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
 * Full `McsConfig` + `McsOutcome`-aware MCS search.
 *
 * `smiles_json` — JSON array of at least 2 SMILES strings.
 * `config_json` — a JSON object with any subset of `McsConfig`'s fields
 * (camelCase keys, all optional, defaulting to `McsConfig::default()`):
 * `matchBonds`, `minAtoms`, `timeoutMs`, `ringMatchesRingOnly`,
 * `completeRingsOnly`, `atomCompare` (`"elements"` | `"any_heavy_atom"` |
 * `"any"`), `bondCompare` (`"order_or_aromatic"` | `"any"`),
 * `matchChiralTag`, `matchCharge`, `matchIsotope`, `maximizeBonds`.
 *
 * Returns a JSON object `{"smiles": string|null, "wasTimedOut": bool}` --
 * `smiles` is `null` when there is no common substructure; `wasTimedOut` is
 * `true` if `timeoutMs` was reached before the search finished exhaustively
 * (the returned `smiles`, if any, is then the best result found so far, not
 * proven optimal).
 * @param {string} smiles_json
 * @param {string} config_json
 * @returns {string}
 */
export function mcs_smiles_json_with_config(smiles_json, config_json) {
    let deferred4_0;
    let deferred4_1;
    try {
        const ptr0 = passStringToWasm0(smiles_json, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ptr1 = passStringToWasm0(config_json, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len1 = WASM_VECTOR_LEN;
        const ret = wasm.mcs_smiles_json_with_config(ptr0, len0, ptr1, len1);
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
 * MCS with ring-awareness constraints.
 *
 * `smiles_json` — JSON array of at least 2 SMILES strings.
 * `ring_matches_ring_only` — ring atoms may only match ring atoms.
 * `complete_rings_only` — partial ring inclusion is removed from the result.
 * Returns the MCS SMILES, or `"null"` when no common substructure was found.
 * @param {string} smiles_json
 * @param {boolean} ring_matches_ring_only
 * @param {boolean} complete_rings_only
 * @returns {string}
 */
export function mcs_smiles_json_with_ring_config(smiles_json, ring_matches_ring_only, complete_rings_only) {
    let deferred3_0;
    let deferred3_1;
    try {
        const ptr0 = passStringToWasm0(smiles_json, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.mcs_smiles_json_with_ring_config(ptr0, len0, ring_matches_ring_only, complete_rings_only);
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
 * MinHash fingerprint (128 hashes) as JSON.
 *
 * Returns `{"num_hashes":128,"hashes":[u64,...]}`.
 * Use `tanimoto_mhfp_smiles` for direct SMILES-to-SMILES similarity.
 * @param {MolHandle} mol
 * @returns {string}
 */
export function mhfp_hashes_json(mol) {
    let deferred1_0;
    let deferred1_1;
    try {
        _assertClass(mol, MolHandle);
        const ret = wasm.mhfp_hashes_json(mol.__wbg_ptr);
        deferred1_0 = ret[0];
        deferred1_1 = ret[1];
        return getStringFromWasm0(ret[0], ret[1]);
    } finally {
        wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
    }
}

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
 * @param {MolHandle} mol
 * @returns {string}
 */
export function minimize_dreiding_json(mol) {
    let deferred1_0;
    let deferred1_1;
    try {
        _assertClass(mol, MolHandle);
        const ret = wasm.minimize_dreiding_json(mol.__wbg_ptr);
        deferred1_0 = ret[0];
        deferred1_1 = ret[1];
        return getStringFromWasm0(ret[0], ret[1]);
    } finally {
        wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
    }
}

/**
 * Minimize geometry using MMFF94 steepest descent (Halgren 1996 full parameters).
 * Generates 3D coords internally if needed.
 * Returns JSON: {"energy":E,"rmsd":R,"converged":true,"iterations":N} or {"error":"..."}.
 * @param {MolHandle} mol
 * @param {number} max_iter
 * @returns {string}
 */
export function minimize_mmff94_json(mol, max_iter) {
    let deferred1_0;
    let deferred1_1;
    try {
        _assertClass(mol, MolHandle);
        const ret = wasm.minimize_mmff94_json(mol.__wbg_ptr, max_iter);
        deferred1_0 = ret[0];
        deferred1_1 = ret[1];
        return getStringFromWasm0(ret[0], ret[1]);
    } finally {
        wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
    }
}

/**
 * Minimize geometry using MMFF94 L-BFGS (faster convergence than steepest descent).
 * Returns JSON: {"energy":E,"rmsd":R,"converged":true,"iterations":N} or {"error":"..."}.
 * @param {MolHandle} mol
 * @param {number} max_iter
 * @returns {string}
 */
export function minimize_mmff94_lbfgs_json(mol, max_iter) {
    let deferred1_0;
    let deferred1_1;
    try {
        _assertClass(mol, MolHandle);
        const ret = wasm.minimize_mmff94_lbfgs_json(mol.__wbg_ptr, max_iter);
        deferred1_0 = ret[0];
        deferred1_1 = ret[1];
        return getStringFromWasm0(ret[0], ret[1]);
    } finally {
        wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
    }
}

/**
 * Minimise a molecule's geometry using the Universal Force Field (UFF).
 *
 * `coords_json` — JSON array of `[x,y,z]` arrays (Å), one per atom.
 * `max_iter` — maximum iterations (0 = default 500).
 *
 * Returns JSON: `{"coords":[[x,y,z],...], "energy":float, "iterations":int, "converged":bool, "sound":bool}`
 * or `{"error":"<msg>"}` on failure. `sound` is all-finite coordinates and
 * no bond stretched past a sane covalent-bond length — independent of
 * `converged`, since steepest descent often reports `converged:false` on
 * geometries that are perfectly fine but simply haven't hit the tight
 * RMS-gradient threshold yet. Check `sound`, not just `converged`, before
 * trusting a result.
 * @param {string} smiles
 * @param {string} coords_json
 * @param {number} max_iter
 * @returns {string}
 */
export function minimize_uff_json(smiles, coords_json, max_iter) {
    let deferred3_0;
    let deferred3_1;
    try {
        const ptr0 = passStringToWasm0(smiles, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ptr1 = passStringToWasm0(coords_json, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len1 = WASM_VECTOR_LEN;
        const ret = wasm.minimize_uff_json(ptr0, len0, ptr1, len1, max_iter);
        deferred3_0 = ret[0];
        deferred3_1 = ret[1];
        return getStringFromWasm0(ret[0], ret[1]);
    } finally {
        wasm.__wbindgen_free(deferred3_0, deferred3_1, 1);
    }
}

/**
 * Cartesian coordinates from an mmCIF file, in the SAME atom order
 * [`mol_from_mmcif`] returns topology for. Returns JSON `[[x,y,z],...]`
 * (Å).
 * @param {string} text
 * @returns {string}
 */
export function mmcif_coords_json(text) {
    let deferred3_0;
    let deferred3_1;
    try {
        const ptr0 = passStringToWasm0(text, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.mmcif_coords_json(ptr0, len0);
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
 * Parse an mmCIF file and return every `_atom_site` field (occupancy,
 * B-factor, chain/residue bookkeeping, formal charge, model number, ...),
 * the unit cell, space group, and any loop column this reader saw but
 * does not model, as JSON: `{"atoms":[{...}],"cell":{...}|null,
 * "space_group":"..."|null,"unhandled_columns":[...]}`. See
 * [`chematic_mol::MmcifAtomRecord`]'s doc comment for each atom field's
 * exact source column and defaulting rule.
 * @param {string} text
 * @returns {string}
 */
export function mmcif_to_json(text) {
    let deferred3_0;
    let deferred3_1;
    try {
        const ptr0 = passStringToWasm0(text, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.mmcif_to_json(ptr0, len0);
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
 * MMFF94 partial charges (BCI table, ±0.1e accuracy) as a JSON array of f64.
 *
 * Uses Bond Charge Increment (BCI) model (Halgren 1996) for 25 common bond types.
 * Returns `[q0, q1, ..., qN]` — one value per heavy atom.
 * Total charge equals the sum of formal charges (charge conserved).
 * @param {MolHandle} mol
 * @returns {string}
 */
export function mmff94_charges_json(mol) {
    let deferred1_0;
    let deferred1_1;
    try {
        _assertClass(mol, MolHandle);
        const ret = wasm.mmff94_charges_json(mol.__wbg_ptr);
        deferred1_0 = ret[0];
        deferred1_1 = ret[1];
        return getStringFromWasm0(ret[0], ret[1]);
    } finally {
        wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
    }
}

/**
 * Compute MMFF94-style atom-typed partial charges (improved over element-pair BCI).
 * Returns JSON: {"charges":[f64,...]} or {"error":"..."}.
 * Uses atom-type classification (Csp3/Ccarbonyl/Ohydroxyl/Oester/Nar/NarH etc.)
 * for better accuracy (~±0.02e) vs element-pair BCI (~±0.05e).
 * @param {MolHandle} mol
 * @returns {string}
 */
export function mmff94_charges_typed_json(mol) {
    let deferred1_0;
    let deferred1_1;
    try {
        _assertClass(mol, MolHandle);
        const ret = wasm.mmff94_charges_typed_json(mol.__wbg_ptr);
        deferred1_0 = ret[0];
        deferred1_1 = ret[1];
        return getStringFromWasm0(ret[0], ret[1]);
    } finally {
        wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
    }
}

/**
 * Compute MMFF94 energy breakdown for EXPLICIT, caller-supplied 3D
 * coordinates -- unlike `mmff94_energy_breakdown_json`, this reads the
 * geometry the caller actually has (e.g. from `pdb_coords_json`,
 * `generate_3d_coords_json`, `generate_3d_etkdg_coords_json`, or an
 * externally computed conformer) instead of silently generating a fresh
 * rule-based one, matching the Python binding's
 * `mol.mmff94_energy_breakdown(coords)` contract (issue #90).
 *
 * `coords_json` -- JSON array of `[x,y,z]` arrays (Å), one per heavy atom,
 * in the same atom order as `mol`.
 *
 * Returns JSON `{"bond":B,"angle":A,"stretch_bend":S,"torsion":T,"oop":O,
 * "vdw":V,"electrostatic":E,"total":X}` at full `f64` round-trip precision
 * (not rounded -- this API exists specifically for oracle comparison
 * against the Python binding, where 4-decimal rounding would mask
 * sub-1e-4 discrepancies), or `{"error":"<msg>"}` on malformed JSON, a
 * non-finite coordinate, or a coordinate-count/atom-count mismatch. Never
 * falls back to a generated conformer.
 * @param {MolHandle} mol
 * @param {string} coords_json
 * @returns {string}
 */
export function mmff94_energy_breakdown_from_coords_json(mol, coords_json) {
    let deferred2_0;
    let deferred2_1;
    try {
        _assertClass(mol, MolHandle);
        const ptr0 = passStringToWasm0(coords_json, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.mmff94_energy_breakdown_from_coords_json(mol.__wbg_ptr, ptr0, len0);
        deferred2_0 = ret[0];
        deferred2_1 = ret[1];
        return getStringFromWasm0(ret[0], ret[1]);
    } finally {
        wasm.__wbindgen_free(deferred2_0, deferred2_1, 1);
    }
}

/**
 * Computes energy on an internally generated conformer.
 * `MolHandle` stores topology only; coordinates previously read from PDB/XYZ
 * are not used by this function.
 * Use [`mmff94_energy_breakdown_from_coords_json`] for explicit coordinates.
 *
 * Returns JSON: {"bond":B,"angle":A,"torsion":T,"vdw":V,"elec":E,"total":X} or {"error":"..."}.
 * @param {MolHandle} mol
 * @returns {string}
 */
export function mmff94_energy_breakdown_json(mol) {
    let deferred1_0;
    let deferred1_1;
    try {
        _assertClass(mol, MolHandle);
        const ret = wasm.mmff94_energy_breakdown_json(mol.__wbg_ptr);
        deferred1_0 = ret[0];
        deferred1_1 = ret[1];
        return getStringFromWasm0(ret[0], ret[1]);
    } finally {
        wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
    }
}

/**
 * Compute MMFF94 partial charges using numeric atom types (Halgren 1996 eq. 15).
 * Returns JSON: {"charges":[-0.28,0.15,...]} or {"error":"..."}.
 * @param {MolHandle} mol
 * @returns {string}
 */
export function mmff94_partial_charges_json(mol) {
    let deferred1_0;
    let deferred1_1;
    try {
        _assertClass(mol, MolHandle);
        const ret = wasm.mmff94_partial_charges_json(mol.__wbg_ptr);
        deferred1_0 = ret[0];
        deferred1_1 = ret[1];
        return getStringFromWasm0(ret[0], ret[1]);
    } finally {
        wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
    }
}

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
 * @param {string} smiles_json
 * @returns {string}
 */
export function mmp_pairs_json(smiles_json) {
    let deferred3_0;
    let deferred3_1;
    try {
        const ptr0 = passStringToWasm0(smiles_json, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.mmp_pairs_json(ptr0, len0);
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
 * Parse a Tripos MOL2 string and return SMILES.
 *
 * Returns `"error:<msg>"` on failure.
 * @param {string} mol2_str
 * @returns {string}
 */
export function mol2_to_smiles(mol2_str) {
    let deferred2_0;
    let deferred2_1;
    try {
        const ptr0 = passStringToWasm0(mol2_str, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.mol2_to_smiles(ptr0, len0);
        deferred2_0 = ret[0];
        deferred2_1 = ret[1];
        return getStringFromWasm0(ret[0], ret[1]);
    } finally {
        wasm.__wbindgen_free(deferred2_0, deferred2_1, 1);
    }
}

/**
 * Parse a MOL V2000 string and return 2D coordinates as a JSON array.
 *
 * Returns `[[x0,y0],[x1,y1],...]` in atom-insertion order.
 * Coordinates are in Ångström as stored in the MOL file.
 * @param {string} mol_block
 * @returns {string}
 */
export function mol_block_coords_json(mol_block) {
    let deferred3_0;
    let deferred3_1;
    try {
        const ptr0 = passStringToWasm0(mol_block, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.mol_block_coords_json(ptr0, len0);
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
 * Serialize a SMILES string directly to a MOL V2000 block with 2D coordinates.
 *
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
 * Rejected wedge/hash stereocenters for a MOL V2000 block, as JSON.
 *
 * Companion to [`mol_from_sdf_block`]/[`mol_block_coords_json`] -- returns
 * `[{"atom_idx":N,"reason":"..."}]`, empty unless a wedge/hash bond was
 * present at some center and got rejected. See
 * `crates/chematic-py/src/formats.rs`'s `from_mol_block_with_diagnostics`
 * for the reason vocabulary (kept identical across bindings).
 * @param {string} mol_block
 * @returns {string}
 */
export function mol_block_stereo_diagnostics_json(mol_block) {
    let deferred3_0;
    let deferred3_1;
    try {
        const ptr0 = passStringToWasm0(mol_block, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.mol_block_stereo_diagnostics_json(ptr0, len0);
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
 * Only the first molecular fragment in the document is returned.
 * Returns a JS error if the document cannot be parsed.
 * @param {string} cdxml
 * @returns {MolHandle}
 */
export function mol_from_cdxml(cdxml) {
    const ptr0 = passStringToWasm0(cdxml, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
    const len0 = WASM_VECTOR_LEN;
    const ret = wasm.mol_from_cdxml(ptr0, len0);
    if (ret[2]) {
        throw takeFromExternrefTable0(ret[1]);
    }
    return MolHandle.__wrap(ret[0]);
}

/**
 * Parse a CML string into a `MolHandle`.
 *
 * Returns a JS error if the CML is invalid (unknown element, bad bond, etc.).
 * @param {string} cml
 * @returns {MolHandle}
 */
export function mol_from_cml(cml) {
    const ptr0 = passStringToWasm0(cml, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
    const len0 = WASM_VECTOR_LEN;
    const ret = wasm.mol_from_cml(ptr0, len0);
    if (ret[2]) {
        throw takeFromExternrefTable0(ret[1]);
    }
    return MolHandle.__wrap(ret[0]);
}

/**
 * Parse a Gaussian Cube file and return a `MolHandle` (topology only --
 * element list, no bonds; Cube carries no bond table). Use
 * [`cube_grid_json`] to recover coordinates, the scalar field, and the
 * grid geometry.
 * @param {string} text
 * @returns {MolHandle}
 */
export function mol_from_cube(text) {
    const ptr0 = passStringToWasm0(text, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
    const len0 = WASM_VECTOR_LEN;
    const ret = wasm.mol_from_cube(ptr0, len0);
    if (ret[2]) {
        throw takeFromExternrefTable0(ret[1]);
    }
    return MolHandle.__wrap(ret[0]);
}

/**
 * Parse an Extended XYZ (extxyz) frame and return a `MolHandle` (topology +
 * element/position only; use [`extxyz_frame_json`] to recover coordinates,
 * cell, per-atom properties and frame metadata in the SAME atom order).
 *
 * A plain XYZ file (free-form comment, no `Lattice=`/`Properties=`) parses
 * too. Only the first frame of a multi-frame file is read.
 *
 * Returns a JS error on parse failure or if the frame exceeds the WASM
 * atom-count limit.
 * @param {string} text
 * @returns {MolHandle}
 */
export function mol_from_extxyz(text) {
    const ptr0 = passStringToWasm0(text, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
    const len0 = WASM_VECTOR_LEN;
    const ret = wasm.mol_from_extxyz(ptr0, len0);
    if (ret[2]) {
        throw takeFromExternrefTable0(ret[1]);
    }
    return MolHandle.__wrap(ret[0]);
}

/**
 * Parse an mmCIF file and return a `MolHandle` (topology only -- element
 * list, no bonds; mmCIF's `_atom_site` category carries no connectivity).
 * Includes every model's atoms if the file has more than one -- use
 * [`mmcif_to_json`] to get each atom's `model_num` for filtering. Use
 * [`mmcif_coords_json`] to recover coordinates in the same atom order.
 * @param {string} text
 * @returns {MolHandle}
 */
export function mol_from_mmcif(text) {
    const ptr0 = passStringToWasm0(text, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
    const len0 = WASM_VECTOR_LEN;
    const ret = wasm.mol_from_mmcif(ptr0, len0);
    if (ret[2]) {
        throw takeFromExternrefTable0(ret[1]);
    }
    return MolHandle.__wrap(ret[0]);
}

/**
 * Parse a MolJSON string into a `MolHandle`.
 *
 * MolJSON is a JSON-based molecular representation designed for LLM
 * (large language model) compatibility.  Returns a JS error on invalid input.
 * @param {string} json
 * @returns {MolHandle}
 */
export function mol_from_moljson(json) {
    const ptr0 = passStringToWasm0(json, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
    const len0 = WASM_VECTOR_LEN;
    const ret = wasm.mol_from_moljson(ptr0, len0);
    if (ret[2]) {
        throw takeFromExternrefTable0(ret[1]);
    }
    return MolHandle.__wrap(ret[0]);
}

/**
 * Parse an ORCA input file (`.inp`) and return a `MolHandle` (topology
 * only -- element list, no bonds; ORCA input carries no bond table).
 * Returns a JS error unless the file's coordinate block is an embedded
 * `* xyz ... *` block -- `xyzfile`/`gzmtfile`/`int` (Z-matrix) blocks
 * carry no atom list to convert, or none is present at all. Use
 * [`orca_input_coords_json`] to recover coordinates + charge +
 * multiplicity in the same atom order, or [`orca_input_to_json`] for the
 * full input (comments/keywords/blocks/any coordinate-block kind).
 * @param {string} text
 * @returns {MolHandle}
 */
export function mol_from_orca_input(text) {
    const ptr0 = passStringToWasm0(text, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
    const len0 = WASM_VECTOR_LEN;
    const ret = wasm.mol_from_orca_input(ptr0, len0);
    if (ret[2]) {
        throw takeFromExternrefTable0(ret[1]);
    }
    return MolHandle.__wrap(ret[0]);
}

/**
 * Parse a PDB file and return a `MolHandle` (topology only; coordinates are
 * discarded -- use [`pdb_coords_json`] to recover them in the SAME atom
 * order, and [`mmff94_energy_breakdown_from_coords_json`] to score them
 * without chematic regenerating a fresh conformer).
 *
 * Uses CONECT records for connectivity if present; otherwise infers bonds from
 * atom distances (the same heuristic as the internal `pdb_to_molecule` function).
 * @param {string} pdb
 * @returns {MolHandle}
 */
export function mol_from_pdb(pdb) {
    const ptr0 = passStringToWasm0(pdb, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
    const len0 = WASM_VECTOR_LEN;
    const ret = wasm.mol_from_pdb(ptr0, len0);
    return MolHandle.__wrap(ret);
}

/**
 * Parse a PQR file and return a `MolHandle` (topology only -- element
 * list inferred per-atom, no bonds; PQR carries no connectivity). Use
 * [`pqr_coords_json`] to recover coordinates in the same atom order.
 * @param {string} text
 * @returns {MolHandle}
 */
export function mol_from_pqr(text) {
    const ptr0 = passStringToWasm0(text, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
    const len0 = WASM_VECTOR_LEN;
    const ret = wasm.mol_from_pqr(ptr0, len0);
    if (ret[2]) {
        throw takeFromExternrefTable0(ret[1]);
    }
    return MolHandle.__wrap(ret[0]);
}

/**
 * Parse a QCSchema `qcschema_molecule` JSON document and return a
 * `MolHandle` (topology + `atomic_numbers`-derived isotopes -- no bonds
 * unless the document's optional `connectivity` list is present, in which
 * case those bond orders are mapped onto the nearest
 * [`chematic_core::BondOrder`]). Use [`qcschema_molecule_coords_json`] to
 * recover coordinates (converted Bohr -> Å) plus molecular
 * charge/multiplicity, in the same atom order.
 * @param {string} json
 * @returns {MolHandle}
 */
export function mol_from_qcschema_molecule(json) {
    const ptr0 = passStringToWasm0(json, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
    const len0 = WASM_VECTOR_LEN;
    const ret = wasm.mol_from_qcschema_molecule(ptr0, len0);
    if (ret[2]) {
        throw takeFromExternrefTable0(ret[1]);
    }
    return MolHandle.__wrap(ret[0]);
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
 * Parse a MOL V3000 block and return a `MolHandle`.
 *
 * Returns a JS error string on parse failure.
 * @param {string} block
 * @returns {MolHandle}
 */
export function mol_from_v3000_block(block) {
    const ptr0 = passStringToWasm0(block, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
    const len0 = WASM_VECTOR_LEN;
    const ret = wasm.mol_from_v3000_block(ptr0, len0);
    if (ret[2]) {
        throw takeFromExternrefTable0(ret[1]);
    }
    return MolHandle.__wrap(ret[0]);
}

/**
 * Parse an XYZ file and return a `MolHandle` (topology only; coordinates are discarded).
 *
 * Returns a JS error on parse failure.
 * @param {string} xyz
 * @returns {MolHandle}
 */
export function mol_from_xyz(xyz) {
    const ptr0 = passStringToWasm0(xyz, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
    const len0 = WASM_VECTOR_LEN;
    const ret = wasm.mol_from_xyz(ptr0, len0);
    if (ret[2]) {
        throw takeFromExternrefTable0(ret[1]);
    }
    return MolHandle.__wrap(ret[0]);
}

/**
 * Return the index that would be assigned to an atom appended to `mol`.
 * @param {MolHandle} mol
 * @returns {number}
 */
export function mol_next_atom_idx(mol) {
    _assertClass(mol, MolHandle);
    const ret = wasm.mol_next_atom_idx(mol.__wbg_ptr);
    return ret >>> 0;
}

/**
 * Rejected wedge/hash stereocenters for a MOL V3000 block, as JSON. Same
 * shape as [`mol_block_stereo_diagnostics_json`] but for V3000 input.
 * @param {string} block
 * @returns {string}
 */
export function mol_v3000_stereo_diagnostics_json(block) {
    let deferred3_0;
    let deferred3_1;
    try {
        const ptr0 = passStringToWasm0(block, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.mol_v3000_stereo_diagnostics_json(ptr0, len0);
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
 * Return a new `MolHandle` with one atom appended.
 *
 * The second return value is the new atom's index (as a JS number).
 * Use `with_atom_added_idx` to retrieve the index.
 * @param {MolHandle} mol
 * @param {string} element_symbol
 * @returns {MolHandle}
 */
export function mol_with_atom_added(mol, element_symbol) {
    _assertClass(mol, MolHandle);
    const ptr0 = passStringToWasm0(element_symbol, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
    const len0 = WASM_VECTOR_LEN;
    const ret = wasm.mol_with_atom_added(mol.__wbg_ptr, ptr0, len0);
    if (ret[2]) {
        throw takeFromExternrefTable0(ret[1]);
    }
    return MolHandle.__wrap(ret[0]);
}

/**
 * Return a new `MolHandle` with the formal charge of atom `idx` changed.
 *
 * Returns a JS error if `idx` is out of range.
 * @param {MolHandle} mol
 * @param {number} idx
 * @param {number} charge
 * @returns {MolHandle}
 */
export function mol_with_atom_charge(mol, idx, charge) {
    _assertClass(mol, MolHandle);
    const ret = wasm.mol_with_atom_charge(mol.__wbg_ptr, idx, charge);
    if (ret[2]) {
        throw takeFromExternrefTable0(ret[1]);
    }
    return MolHandle.__wrap(ret[0]);
}

/**
 * Return a new `MolHandle` with the element of atom `idx` changed.
 *
 * `element_symbol` — periodic-table symbol, e.g. `"N"`, `"O"`, `"Cl"`.
 * Returns a JS error if `idx` is out of range or the symbol is unknown.
 * @param {MolHandle} mol
 * @param {number} idx
 * @param {string} element_symbol
 * @returns {MolHandle}
 */
export function mol_with_atom_element(mol, idx, element_symbol) {
    _assertClass(mol, MolHandle);
    const ptr0 = passStringToWasm0(element_symbol, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
    const len0 = WASM_VECTOR_LEN;
    const ret = wasm.mol_with_atom_element(mol.__wbg_ptr, idx, ptr0, len0);
    if (ret[2]) {
        throw takeFromExternrefTable0(ret[1]);
    }
    return MolHandle.__wrap(ret[0]);
}

/**
 * Return a new `MolHandle` with atom `idx` and all its bonds removed.
 *
 * Atom indices above `idx` shift down by 1.  Returns a JS error if `idx`
 * is out of range.
 * @param {MolHandle} mol
 * @param {number} idx
 * @returns {MolHandle}
 */
export function mol_with_atom_removed(mol, idx) {
    _assertClass(mol, MolHandle);
    const ret = wasm.mol_with_atom_removed(mol.__wbg_ptr, idx);
    if (ret[2]) {
        throw takeFromExternrefTable0(ret[1]);
    }
    return MolHandle.__wrap(ret[0]);
}

/**
 * Return a new `MolHandle` with one bond added between `a` and `b`.
 *
 * `order` — 1 = single, 2 = double, 3 = triple.
 * Returns a JS error if the bond already exists or `a == b`.
 * @param {MolHandle} mol
 * @param {number} a
 * @param {number} b
 * @param {number} order
 * @returns {MolHandle}
 */
export function mol_with_bond_added(mol, a, b, order) {
    _assertClass(mol, MolHandle);
    const ret = wasm.mol_with_bond_added(mol.__wbg_ptr, a, b, order);
    if (ret[2]) {
        throw takeFromExternrefTable0(ret[1]);
    }
    return MolHandle.__wrap(ret[0]);
}

/**
 * Return a new `MolHandle` with bond `idx` removed.
 *
 * Atom indices are unchanged; bond indices above `idx` shift down.
 * Returns a JS error if `idx` is out of range.
 * @param {MolHandle} mol
 * @param {number} idx
 * @returns {MolHandle}
 */
export function mol_with_bond_removed(mol, idx) {
    _assertClass(mol, MolHandle);
    const ret = wasm.mol_with_bond_removed(mol.__wbg_ptr, idx);
    if (ret[2]) {
        throw takeFromExternrefTable0(ret[1]);
    }
    return MolHandle.__wrap(ret[0]);
}

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
 * @param {string} smiles
 * @returns {string}
 */
export function molecule_report_json(smiles) {
    let deferred3_0;
    let deferred3_1;
    try {
        const ptr0 = passStringToWasm0(smiles, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.molecule_report_json(ptr0, len0);
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
 * MQN descriptor (42 integer values: Molecular Quantum Numbers).
 * @param {MolHandle} mol
 * @returns {string}
 */
export function mqn_json(mol) {
    let deferred1_0;
    let deferred1_1;
    try {
        _assertClass(mol, MolHandle);
        const ret = wasm.mqn_json(mol.__wbg_ptr);
        deferred1_0 = ret[0];
        deferred1_1 = ret[1];
        return getStringFromWasm0(ret[0], ret[1]);
    } finally {
        wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
    }
}

/**
 * Per-atom molar refractivity contributions as a JSON array of f64.
 * @param {MolHandle} mol
 * @returns {string}
 */
export function mr_per_atom_json(mol) {
    let deferred1_0;
    let deferred1_1;
    try {
        _assertClass(mol, MolHandle);
        const ret = wasm.mr_per_atom_json(mol.__wbg_ptr);
        deferred1_0 = ret[0];
        deferred1_1 = ret[1];
        return getStringFromWasm0(ret[0], ret[1]);
    } finally {
        wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
    }
}

/**
 * Murcko scaffold of `mol` — the ring system plus linkers, side-chains removed.
 *
 * Returns a new `MolHandle`.  For acyclic molecules returns an empty molecule.
 * @param {MolHandle} mol
 * @returns {MolHandle}
 */
export function murcko_scaffold(mol) {
    _assertClass(mol, MolHandle);
    const ret = wasm.murcko_scaffold(mol.__wbg_ptr);
    return MolHandle.__wrap(ret);
}

/**
 * Find the k nearest neighbours of a query SMILES in a list of db SMILES.
 *
 * `db_smiles_json`: JSON array of SMILES strings, e.g. `["CC","c1ccccc1"]`.
 * Returns JSON: `[{"index":0,"tanimoto":0.95},...]` sorted by descending Tanimoto.
 * Returns `"error:<msg>"` on parse failure.
 * @param {string} query_smiles
 * @param {string} db_smiles_json
 * @param {number} k
 * @returns {string}
 */
export function nearest_neighbors_json(query_smiles, db_smiles_json, k) {
    let deferred3_0;
    let deferred3_1;
    try {
        const ptr0 = passStringToWasm0(query_smiles, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ptr1 = passStringToWasm0(db_smiles_json, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len1 = WASM_VECTOR_LEN;
        const ret = wasm.nearest_neighbors_json(ptr0, len0, ptr1, len1, k);
        deferred3_0 = ret[0];
        deferred3_1 = ret[1];
        return getStringFromWasm0(ret[0], ret[1]);
    } finally {
        wasm.__wbindgen_free(deferred3_0, deferred3_1, 1);
    }
}

/**
 * Neutralize formal charges on `mol` by proton addition/removal.
 *
 * Returns a new `MolHandle` with all formal charges set to zero where possible.
 * @param {MolHandle} mol
 * @returns {MolHandle}
 */
export function neutralize_charges(mol) {
    _assertClass(mol, MolHandle);
    const ret = wasm.neutralize_charges(mol.__wbg_ptr);
    return MolHandle.__wrap(ret);
}

/**
 * Parse and re-serialize CXSMILES, preserving supported CX metadata.
 * Returns error if atom count exceeds 10,000.
 * @param {string} s
 * @returns {string}
 */
export function normalize_cxsmiles(s) {
    let deferred3_0;
    let deferred3_1;
    try {
        const ptr0 = passStringToWasm0(s, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.normalize_cxsmiles(ptr0, len0);
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
 * Parse and re-serialise a reaction SMILES string, returning the normalised form.
 *
 * Useful for validating reaction SMILES and obtaining a canonical representation.
 * Returns a JS error on parse failure.
 * @param {string} rxn_smiles
 * @returns {string}
 */
export function normalize_reaction_smiles(rxn_smiles) {
    let deferred3_0;
    let deferred3_1;
    try {
        const ptr0 = passStringToWasm0(rxn_smiles, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.normalize_reaction_smiles(ptr0, len0);
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
 * Parse an OpenDX (APBS scalar-field subset) file and return its full
 * [`chematic_mol::VolumetricGrid`] as JSON (same shape as
 * [`cube_grid_json`]; `atoms` is always empty -- OpenDX has no atom
 * section). No `mol_from_opendx` is provided: an OpenDX grid never carries
 * atoms, so a `MolHandle` from one would always be empty and is not a
 * useful binding.
 * @param {string} text
 * @returns {string}
 */
export function opendx_grid_json(text) {
    let deferred3_0;
    let deferred3_1;
    try {
        const ptr0 = passStringToWasm0(text, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.opendx_grid_json(ptr0, len0);
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
 * `[nx, ny, nz]` for an OpenDX file's grid, as a `Uint32Array`.
 * @param {string} text
 * @returns {Uint32Array}
 */
export function opendx_shape_u32(text) {
    const ptr0 = passStringToWasm0(text, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
    const len0 = WASM_VECTOR_LEN;
    const ret = wasm.opendx_shape_u32(ptr0, len0);
    if (ret[2]) {
        throw takeFromExternrefTable0(ret[1]);
    }
    return takeFromExternrefTable0(ret[0]);
}

/**
 * Flat `values` from an OpenDX file's grid, as a `Float64Array` -- same
 * data [`opendx_grid_json`]'s `"values"` field carries.
 * @param {string} text
 * @returns {Float64Array}
 */
export function opendx_values_f64(text) {
    const ptr0 = passStringToWasm0(text, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
    const len0 = WASM_VECTOR_LEN;
    const ret = wasm.opendx_values_f64(ptr0, len0);
    if (ret[2]) {
        throw takeFromExternrefTable0(ret[1]);
    }
    return takeFromExternrefTable0(ret[0]);
}

/**
 * Coordinates + charge + multiplicity from an ORCA input file's embedded
 * `* xyz ... *` block, in the SAME atom order [`mol_from_orca_input`]
 * returns topology for. Returns JSON
 * `{"coords":[[x,y,z],...],"charge":0,"multiplicity":1}`, or a JS error
 * under the same conditions as [`mol_from_orca_input`].
 * @param {string} text
 * @returns {string}
 */
export function orca_input_coords_json(text) {
    let deferred3_0;
    let deferred3_1;
    try {
        const ptr0 = passStringToWasm0(text, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.orca_input_coords_json(ptr0, len0);
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
 * Parse an ORCA input file and return every field as JSON:
 * `{"comments":[...],"keywords":[...],
 * "blocks":[{"name":"scf","raw":"...","has_end":true},...],
 * "coords":{"type":"xyz"|"xyzfile"|"gzmtfile"|"internal",...}|null}`.
 * See [`chematic_mol::OrcaInput`]'s doc comment for each field's meaning.
 * @param {string} text
 * @returns {string}
 */
export function orca_input_to_json(text) {
    let deferred3_0;
    let deferred3_1;
    try {
        const ptr0 = passStringToWasm0(text, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.orca_input_to_json(ptr0, len0);
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
 * Parse an ORCA output file (`.out`/`.log`) and return every extracted
 * field as JSON: `{"charge":N|null,"multiplicity":N|null,
 * "final_energy_hartree":N|null,
 * "trajectory":[{"elements":[...],"coords":[[x,y,z],...]},...],
 * "frequencies_cm1":[...],
 * "termination":{"kind":"normal"|"error"|"incomplete","detail":"..."?},
 * "optimization_convergence":"not_requested"|"converged"|"not_converged"|"unknown"}`.
 * No writer is provided -- an ORCA output file is a job log, not a
 * document this crate constructs.
 * @param {string} text
 * @returns {string}
 */
export function orca_output_to_json(text) {
    let deferred3_0;
    let deferred3_1;
    try {
        const ptr0 = passStringToWasm0(text, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.orca_output_to_json(ptr0, len0);
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
 * PAINS structural alert names matched by `mol` as a JSON array.
 *
 * Returns `[]` when no alerts fire, or e.g. `["ene_six_het_A(483)"]`.
 * Use alongside `pains_passes()` to know *which* alerts triggered.
 * @param {MolHandle} mol
 * @returns {string}
 */
export function pains_matches_json(mol) {
    let deferred1_0;
    let deferred1_1;
    try {
        _assertClass(mol, MolHandle);
        const ret = wasm.pains_matches_json(mol.__wbg_ptr);
        deferred1_0 = ret[0];
        deferred1_1 = ret[1];
        return getStringFromWasm0(ret[0], ret[1]);
    } finally {
        wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
    }
}

/**
 * Parse CXSMARTS and return preserved metadata as JSON.
 * Returns error if atom count exceeds 10,000.
 * @param {string} s
 * @returns {string}
 */
export function parse_cxsmarts_json(s) {
    let deferred3_0;
    let deferred3_1;
    try {
        const ptr0 = passStringToWasm0(s, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.parse_cxsmarts_json(ptr0, len0);
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
 * Parse CXSMILES and return preserved metadata as JSON.
 *
 * Supported CX fields: atom labels (`$...$`), `atomProp`, atom radicals (`^n:`),
 * and zero-order bonds (`Z:`). The `cxsmiles` field is a re-serialized
 * round-trip form using the supported fields.
 * Returns error if atom count exceeds 10,000.
 * @param {string} s
 * @returns {string}
 */
export function parse_cxsmiles_json(s) {
    let deferred3_0;
    let deferred3_1;
    try {
        const ptr0 = passStringToWasm0(s, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.parse_cxsmiles_json(ptr0, len0);
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
 * Parse a SMILES string into a `MolHandle`.
 *
 * Returns a JS error string on parse failure or if atom count exceeds 10,000.
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
 * Extract the atomic coordinates from a PDB block, in the SAME atom order
 * `mol_from_pdb` returns topology for (both read the identical underlying
 * parse via `parse_pdb_molecule_and_coords`, so atom-index correspondence
 * between the two calls is structural, not just conventional -- issue #90).
 *
 * Returns JSON `[[x,y,z],...]` (full `f64` precision, not rounded -- for
 * oracle comparison against the Python binding) or `{"error":"<msg>"}` --
 * including when a coordinate field parsed to a non-finite value (see
 * `coords_all_finite`'s doc comment), rather than emitting invalid JSON.
 * @param {string} pdb
 * @returns {string}
 */
export function pdb_coords_json(pdb) {
    let deferred2_0;
    let deferred2_1;
    try {
        const ptr0 = passStringToWasm0(pdb, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.pdb_coords_json(ptr0, len0);
        deferred2_0 = ret[0];
        deferred2_1 = ret[1];
        return getStringFromWasm0(ret[0], ret[1]);
    } finally {
        wasm.__wbindgen_free(deferred2_0, deferred2_1, 1);
    }
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
 * Detect pharmacophore features for virtual screening and lead optimization.
 * Returns JSON array of features: [{type, atom_idx, neighbor_count}, ...]
 * @param {MolHandle} mol
 * @returns {string}
 */
export function pharmacophore_features_json(mol) {
    let deferred1_0;
    let deferred1_1;
    try {
        _assertClass(mol, MolHandle);
        const ret = wasm.pharmacophore_features_json(mol.__wbg_ptr);
        deferred1_0 = ret[0];
        deferred1_1 = ret[1];
        return getStringFromWasm0(ret[0], ret[1]);
    } finally {
        wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
    }
}

/**
 * Compute 2D pharmacophore fingerprint (2048 bits) as a JSON feature count summary.
 * Returns simplified JSON with feature type counts: {Donor, Acceptor, Aromatic, Hydrophobic, Positive, Negative}
 * @param {MolHandle} mol
 * @returns {string}
 */
export function pharmacophore_fp_2d_summary(mol) {
    let deferred1_0;
    let deferred1_1;
    try {
        _assertClass(mol, MolHandle);
        const ret = wasm.pharmacophore_fp_2d_summary(mol.__wbg_ptr);
        deferred1_0 = ret[0];
        deferred1_1 = ret[1];
        return getStringFromWasm0(ret[0], ret[1]);
    } finally {
        wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
    }
}

/**
 * Compute 3D pharmacophore fingerprint from generated 3D coordinates.
 * Returns simplified JSON with feature type counts (3D-aware version).
 * @param {MolHandle} mol
 * @returns {string}
 */
export function pharmacophore_fp_3d_summary(mol) {
    let deferred1_0;
    let deferred1_1;
    try {
        _assertClass(mol, MolHandle);
        const ret = wasm.pharmacophore_fp_3d_summary(mol.__wbg_ptr);
        deferred1_0 = ret[0];
        deferred1_1 = ret[1];
        return getStringFromWasm0(ret[0], ret[1]);
    } finally {
        wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
    }
}

/**
 * Returns a ready-to-use `embed_pipeline_v2_json` config JSON string for the
 * "stereo-safe" configuration (issue #291/#383): `stereoPolicy:
 * "repair_and_verify"`, `enforceChirality: true`, and
 * `expandImplicitHThroughPipeline: true` together -- the exact combination
 * measured to correctly handle ring-fused declared stereocenters (e.g.
 * testosterone, cholesterol) that `enforceChirality` alone cannot repair.
 * Mirrors `PipelineV2Config::stereo_safe`/the Python binding's
 * `PipelineV2Config.stereo_safe(...)` -- prefer this over setting those three
 * fields individually: they only work correctly as a set, and forgetting one
 * silently falls back to a configuration issue #291 measured as unsound for
 * that molecule class. `forceFieldPolicy`/`ringTorsionPolicy` are still
 * required, explicit arguments; everything else takes the same conservative
 * defaults `embed_pipeline_v2_json`'s own documented examples do. The caller
 * may parse and further override individual fields before passing the result
 * to `embed_pipeline_v2_json` (e.g. a different `embedSeed`).
 *
 * Never throws, matching `embed_pipeline_v2_json`'s own convention: an
 * unknown `force_field`/`ring_torsion_policy` string returns the same
 * `{"ok": false, "error": {...}}` shape `embed_pipeline_v2_json` would for an
 * invalid config, tagged `schemaVersion: 1`.
 * @param {string} force_field
 * @param {string} ring_torsion_policy
 * @returns {string}
 */
export function pipeline_v2_stereo_safe_config_json(force_field, ring_torsion_policy) {
    let deferred3_0;
    let deferred3_1;
    try {
        const ptr0 = passStringToWasm0(force_field, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ptr1 = passStringToWasm0(ring_torsion_policy, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len1 = WASM_VECTOR_LEN;
        const ret = wasm.pipeline_v2_stereo_safe_config_json(ptr0, len0, ptr1, len1);
        deferred3_0 = ret[0];
        deferred3_1 = ret[1];
        return getStringFromWasm0(ret[0], ret[1]);
    } finally {
        wasm.__wbindgen_free(deferred3_0, deferred3_1, 1);
    }
}

/**
 * Cartesian coordinates from a PQR file, in the SAME atom order
 * [`mol_from_pqr`] returns topology for. Returns JSON `[[x,y,z],...]` (Å).
 * @param {string} text
 * @returns {string}
 */
export function pqr_coords_json(text) {
    let deferred3_0;
    let deferred3_1;
    try {
        const ptr0 = passStringToWasm0(text, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.pqr_coords_json(ptr0, len0);
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
 * Infer an element from a PQR atom name (see
 * [`chematic_mol::infer_element`]'s doc comment for the heuristic).
 * Returns `undefined` (JS) / `None` if no element could be inferred.
 * @param {string} group_pdb
 * @param {string} res_name
 * @param {string} atom_name
 * @returns {string | undefined}
 */
export function pqr_infer_element(group_pdb, res_name, atom_name) {
    const ptr0 = passStringToWasm0(group_pdb, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
    const len0 = WASM_VECTOR_LEN;
    const ptr1 = passStringToWasm0(res_name, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
    const len1 = WASM_VECTOR_LEN;
    const ptr2 = passStringToWasm0(atom_name, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
    const len2 = WASM_VECTOR_LEN;
    const ret = wasm.pqr_infer_element(ptr0, len0, ptr1, len1, ptr2, len2);
    let v4;
    if (ret[0] !== 0) {
        v4 = getStringFromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 1, 1);
    }
    return v4;
}

/**
 * Parse a PQR file and return every field (charge, radius, chain,
 * residue, inferred element, ...) as JSON: `{"atoms":[{...}]}`. See
 * [`chematic_mol::PqrAtomRecord`]'s doc comment for each field's meaning.
 * @param {string} text
 * @returns {string}
 */
export function pqr_to_json(text) {
    let deferred3_0;
    let deferred3_1;
    try {
        const ptr0 = passStringToWasm0(text, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.pqr_to_json(ptr0, len0);
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
 * Predict pKa for all ionizable sites in a molecule.
 *
 * Returns a JSON array: `[{"atom_idx":8,"pka":4.0,"type":"acid","group":"carboxylic_acid"},...]`
 *
 * Returns `[]` if no ionizable sites are found, or `{"error":"..."}` on parse failure.
 * @param {string} smiles
 * @returns {string}
 */
export function predict_pka_json(smiles) {
    let deferred2_0;
    let deferred2_1;
    try {
        const ptr0 = passStringToWasm0(smiles, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.predict_pka_json(ptr0, len0);
        deferred2_0 = ret[0];
        deferred2_1 = ret[1];
        return getStringFromWasm0(ret[0], ret[1]);
    } finally {
        wasm.__wbindgen_free(deferred2_0, deferred2_1, 1);
    }
}

/**
 * Coordinates (Å) plus molecular charge/multiplicity from a QCSchema
 * `qcschema_molecule` document, in the SAME atom order
 * [`mol_from_qcschema_molecule`] returns topology for. Returns JSON
 * `{"coords":[[x,y,z],...],"molecular_charge":0.0,"molecular_multiplicity":1}`.
 * @param {string} json
 * @returns {string}
 */
export function qcschema_molecule_coords_json(json) {
    let deferred3_0;
    let deferred3_1;
    try {
        const ptr0 = passStringToWasm0(json, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.qcschema_molecule_coords_json(ptr0, len0);
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
 * Parse a QCSchema `qcschema_input`/`qc_schema_input` JSON document
 * (molecule + driver + model + keywords) and re-emit it, validating and
 * canonicalizing field defaults in the process (e.g. a missing
 * `schema_name`/`schema_version` is filled in). Job-level fields
 * (`driver`, `model`, `keywords`, `protocols`, `extras`) are round-tripped
 * opaquely -- this binding validates/reformats the document; it does not
 * expose a separate JS-facing accessor for each field (out of scope for
 * this first pass, see module docs' "None of these formats carry a bond
 * table" section for the analogous molecule-centric scope choice made
 * elsewhere in this file).
 * @param {string} json
 * @returns {string}
 */
export function qcschema_validate_atomic_input(json) {
    let deferred3_0;
    let deferred3_1;
    try {
        const ptr0 = passStringToWasm0(json, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.qcschema_validate_atomic_input(ptr0, len0);
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
 * Like [`qcschema_validate_atomic_input`], for a QCSchema
 * `qcschema_output`/`qc_schema_output` (`AtomicResult`) document.
 * @param {string} json
 * @returns {string}
 */
export function qcschema_validate_atomic_result(json) {
    let deferred3_0;
    let deferred3_1;
    try {
        const ptr0 = passStringToWasm0(json, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.qcschema_validate_atomic_result(ptr0, len0);
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
 * @param {string} smiles
 * @param {number} count
 * @param {bigint} seed
 * @returns {string}
 */
export function random_smiles_json(smiles, count, seed) {
    let deferred3_0;
    let deferred3_1;
    try {
        const ptr0 = passStringToWasm0(smiles, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.random_smiles_json(ptr0, len0, count, seed);
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
 * RDKit-bit-exact ECFP4 (radius=2, 2048 bits, `useChirality=false`, `useBondTypes=true`,
 * RDKit's default atom invariant) as a bit-packed byte vector (256 bytes = 2048 bits).
 *
 * Bit-for-bit identical to
 * `rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048).GetFingerprint(mol)`
 * for every input this preprocessing handles. **Not** the same bits as `ecfp4_bitvec`
 * (that path uses chematic's own FNV-1a hash and is not RDKit-bit-compatible by
 * design -- the two are never silently interchanged).
 *
 * Returns a JS error (its string carrying `RdkitMorganError`'s `Display` text, e.g.
 * `"rdkit-exact ecfp4: aromaticity: ..."`) if RDKit-parity aromaticity preprocessing
 * fails, never a silent fallback to `ecfp4_bitvec`'s Hückel-based engine -- the two
 * engines are not bit-compatible, so a silent substitution would look successful
 * while actually returning the wrong hash. See `docs/rfcs/ecfp4_bitexact_api_rfc.md`.
 * @param {MolHandle} mol
 * @returns {Uint8Array}
 */
export function rdkit_ecfp4_bitvec(mol) {
    _assertClass(mol, MolHandle);
    const ret = wasm.rdkit_ecfp4_bitvec(mol.__wbg_ptr);
    if (ret[3]) {
        throw takeFromExternrefTable0(ret[2]);
    }
    var v1 = getArrayU8FromWasm0(ret[0], ret[1]).slice();
    wasm.__wbindgen_free(ret[0], ret[1] * 1, 1);
    return v1;
}

/**
 * Same fingerprint as `rdkit_ecfp4_bitvec`, plus the raw (unfolded) data behind it, as
 * JSON: `{"fingerprint":[u8,...],"sparseCounts":{"rawId":count,...},
 * "rawBitInfo":{"rawId":[[atomIdx,radius],...],...},
 * "foldedBitInfo":{"bit":[[atomIdx,radius],...],...}}`.
 *
 * Returns a JS error on the same preprocessing failures as `rdkit_ecfp4_bitvec`.
 * @param {MolHandle} mol
 * @returns {string}
 */
export function rdkit_ecfp4_detail_json(mol) {
    let deferred2_0;
    let deferred2_1;
    try {
        _assertClass(mol, MolHandle);
        const ret = wasm.rdkit_ecfp4_detail_json(mol.__wbg_ptr);
        var ptr1 = ret[0];
        var len1 = ret[1];
        if (ret[3]) {
            ptr1 = 0; len1 = 0;
            throw takeFromExternrefTable0(ret[2]);
        }
        deferred2_0 = ptr1;
        deferred2_1 = len1;
        return getStringFromWasm0(ptr1, len1);
    } finally {
        wasm.__wbindgen_free(deferred2_0, deferred2_1, 1);
    }
}

/**
 * RDKit-bit-exact Morgan/ECFP fingerprint at a caller-chosen radius/bit-width, as a
 * bit-packed byte vector (`nbits / 8` bytes, LSB-first).
 *
 * `radius` must be 0, 1, 2 (`rdkit_ecfp4_bitvec`'s ECFP4), or 3. `nbits` must be one
 * of 128, 256, 512, 1024, or 2048. Each of these 20 combinations is independently
 * re-verified against a live RDKit oracle (not assumed to generalize from
 * radius=2/2048 bits alone) -- see `validation/ecfp4_rdkit_stable_api_fixtures.json`.
 * An unsupported value returns a JS error rather than being silently coerced to the
 * nearest supported one.
 * @param {MolHandle} mol
 * @param {number} radius
 * @param {number} nbits
 * @returns {Uint8Array}
 */
export function rdkit_ecfp_config_bitvec(mol, radius, nbits) {
    _assertClass(mol, MolHandle);
    const ret = wasm.rdkit_ecfp_config_bitvec(mol.__wbg_ptr, radius, nbits);
    if (ret[3]) {
        throw takeFromExternrefTable0(ret[2]);
    }
    var v1 = getArrayU8FromWasm0(ret[0], ret[1]).slice();
    wasm.__wbindgen_free(ret[0], ret[1] * 1, 1);
    return v1;
}

/**
 * Chirality-enabled variant of [`rdkit_ecfp_config_bitvec`]. E/Z bond stereo
 * is not included by this API yet.
 * @param {MolHandle} mol
 * @param {number} radius
 * @param {number} nbits
 * @returns {Uint8Array}
 */
export function rdkit_ecfp_config_chiral_bitvec(mol, radius, nbits) {
    _assertClass(mol, MolHandle);
    const ret = wasm.rdkit_ecfp_config_chiral_bitvec(mol.__wbg_ptr, radius, nbits);
    if (ret[3]) {
        throw takeFromExternrefTable0(ret[2]);
    }
    var v1 = getArrayU8FromWasm0(ret[0], ret[1]).slice();
    wasm.__wbindgen_free(ret[0], ret[1] * 1, 1);
    return v1;
}

/**
 * Chirality-enabled variant of [`rdkit_ecfp_config_detail_json`].
 * @param {MolHandle} mol
 * @param {number} radius
 * @param {number} nbits
 * @returns {string}
 */
export function rdkit_ecfp_config_chiral_detail_json(mol, radius, nbits) {
    let deferred2_0;
    let deferred2_1;
    try {
        _assertClass(mol, MolHandle);
        const ret = wasm.rdkit_ecfp_config_chiral_detail_json(mol.__wbg_ptr, radius, nbits);
        var ptr1 = ret[0];
        var len1 = ret[1];
        if (ret[3]) {
            ptr1 = 0; len1 = 0;
            throw takeFromExternrefTable0(ret[2]);
        }
        deferred2_0 = ptr1;
        deferred2_1 = len1;
        return getStringFromWasm0(ptr1, len1);
    } finally {
        wasm.__wbindgen_free(deferred2_0, deferred2_1, 1);
    }
}

/**
 * Same fingerprint as `rdkit_ecfp_config_bitvec`, plus the raw (unfolded) data -- see
 * `rdkit_ecfp4_detail_json` for the JSON shape (identical, generalized to this
 * function's `radius`/`nbits`).
 * @param {MolHandle} mol
 * @param {number} radius
 * @param {number} nbits
 * @returns {string}
 */
export function rdkit_ecfp_config_detail_json(mol, radius, nbits) {
    let deferred2_0;
    let deferred2_1;
    try {
        _assertClass(mol, MolHandle);
        const ret = wasm.rdkit_ecfp_config_detail_json(mol.__wbg_ptr, radius, nbits);
        var ptr1 = ret[0];
        var len1 = ret[1];
        if (ret[3]) {
            ptr1 = 0; len1 = 0;
            throw takeFromExternrefTable0(ret[2]);
        }
        deferred2_0 = ptr1;
        deferred2_1 = len1;
        return getStringFromWasm0(ptr1, len1);
    } finally {
        wasm.__wbindgen_free(deferred2_0, deferred2_1, 1);
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
 * Single-step retrosynthetic disconnection (issue #91).
 *
 * Thin wrapper around [`chematic_rxn::retro::retro_disconnect`] -- applies
 * the same built-in 60-template SMIRKS library and returns identical
 * disconnections (same templates, same precursor sets, same ordering:
 * fewest precursors first) as the Rust and Python (`Mol.retro_disconnect()`)
 * APIs. This function changes nothing about the underlying algorithm; it
 * only serializes the result to JSON.
 *
 * `max_results` -- cap on returned disconnections (0 = the WASM safety cap).
 *
 * `reaction_class` -- filter to a single reaction class, or `""` for all
 * classes. Valid values: `"AmideBond"`, `"Ester"`, `"Ether"`, `"CNBond"`,
 * `"CCBond"`, `"CSBond"`, `"Other"`. An unrecognized non-empty value is a
 * JS error (not silently ignored).
 *
 * JSON schema: array of
 * `{"template":str,"reaction_class":str,"precursors":[str,...],"sa_scores":[number,...],"max_sa_score":number}`
 * -- same field names as the Python binding's dict output. Returns `[]`
 * when no template matches the molecule (e.g. it has no disconnectable
 * bond the template library recognizes) -- a valid, non-error result,
 * distinct from the `reaction_class` validation error above.
 * @param {MolHandle} mol
 * @param {number} max_results
 * @param {string} reaction_class
 * @returns {string}
 */
export function retro_disconnect_json(mol, max_results, reaction_class) {
    let deferred3_0;
    let deferred3_1;
    try {
        _assertClass(mol, MolHandle);
        const ptr0 = passStringToWasm0(reaction_class, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.retro_disconnect_json(mol.__wbg_ptr, max_results, ptr0, len0);
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
 * @param {string} smiles_json
 * @param {string} core_smarts
 * @returns {string}
 */
export function rgroup_decompose_json(smiles_json, core_smarts) {
    let deferred4_0;
    let deferred4_1;
    try {
        const ptr0 = passStringToWasm0(smiles_json, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ptr1 = passStringToWasm0(core_smarts, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len1 = WASM_VECTOR_LEN;
        const ret = wasm.rgroup_decompose_json(ptr0, len0, ptr1, len1);
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
 * Ring family classification and detection as JSON.
 * Returns an array of ring families with their atoms, ring indices, and topology kind.
 * @param {MolHandle} mol
 * @returns {string}
 */
export function ring_families_json(mol) {
    let deferred2_0;
    let deferred2_1;
    try {
        _assertClass(mol, MolHandle);
        const ret = wasm.ring_families_json(mol.__wbg_ptr);
        var ptr1 = ret[0];
        var len1 = ret[1];
        if (ret[3]) {
            ptr1 = 0; len1 = 0;
            throw takeFromExternrefTable0(ret[2]);
        }
        deferred2_0 = ptr1;
        deferred2_1 = len1;
        return getStringFromWasm0(ptr1, len1);
    } finally {
        wasm.__wbindgen_free(deferred2_0, deferred2_1, 1);
    }
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
 * Screen a batch of SMILES strings (JSON string output).
 * Returns per-record results including pass/fail with error details.
 * Includes MaxMin diversity picking and Butina clustering by default.
 * The input is limited to 1 MiB, 1,024 records, and 10,000 atoms per parsed
 * molecule.
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
 * @param {string} smiles_batch
 * @param {string} delimiter
 * @returns {string}
 */
export function screen_smiles_json(smiles_batch, delimiter) {
    let deferred3_0;
    let deferred3_1;
    try {
        const ptr0 = passStringToWasm0(smiles_batch, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ptr1 = passStringToWasm0(delimiter, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len1 = WASM_VECTOR_LEN;
        const ret = wasm.screen_smiles_json(ptr0, len0, ptr1, len1);
        deferred3_0 = ret[0];
        deferred3_1 = ret[1];
        return getStringFromWasm0(ret[0], ret[1]);
    } finally {
        wasm.__wbindgen_free(deferred3_0, deferred3_1, 1);
    }
}

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
 * @param {string} smiles_json
 * @param {string} names_json
 * @param {string} props_json
 * @returns {string}
 */
export function sdf_from_records_json(smiles_json, names_json, props_json) {
    let deferred5_0;
    let deferred5_1;
    try {
        const ptr0 = passStringToWasm0(smiles_json, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ptr1 = passStringToWasm0(names_json, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len1 = WASM_VECTOR_LEN;
        const ptr2 = passStringToWasm0(props_json, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len2 = WASM_VECTOR_LEN;
        const ret = wasm.sdf_from_records_json(ptr0, len0, ptr1, len1, ptr2, len2);
        var ptr4 = ret[0];
        var len4 = ret[1];
        if (ret[3]) {
            ptr4 = 0; len4 = 0;
            throw takeFromExternrefTable0(ret[2]);
        }
        deferred5_0 = ptr4;
        deferred5_1 = len4;
        return getStringFromWasm0(ptr4, len4);
    } finally {
        wasm.__wbindgen_free(deferred5_0, deferred5_1, 1);
    }
}

/**
 * Parse an SDF string and return a JSON array of record objects.
 *
 * Each record has the shape:
 * ```json
 * {"smiles":"CC(=O)O","name":"aspirin","properties":{"MW":"180.2","Activity":"high"},"stereo_diagnostics":[]}
 * ```
 *
 * Invalid records are represented as `null`.  SD data fields are included in
 * `properties`; multi-line values are joined with `\n`. `stereo_diagnostics`
 * is a list of `{"atom_idx":N,"reason":"..."}` objects, one per rejected
 * wedge/hash center (see [`mol_block_stereo_diagnostics_json`] for the
 * reason vocabulary) -- empty unless a wedge/hash bond was present at some
 * center and got rejected.
 * @param {string} sdf
 * @returns {string}
 */
export function sdf_to_records_json(sdf) {
    let deferred2_0;
    let deferred2_1;
    try {
        const ptr0 = passStringToWasm0(sdf, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.sdf_to_records_json(ptr0, len0);
        deferred2_0 = ret[0];
        deferred2_1 = ret[1];
        return getStringFromWasm0(ret[0], ret[1]);
    } finally {
        wasm.__wbindgen_free(deferred2_0, deferred2_1, 1);
    }
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
 * @param {string} smiles
 * @param {number} a
 * @param {number} b
 * @param {number} c
 * @param {number} d
 * @param {number} angle_deg
 * @returns {string}
 */
export function set_dihedral_json(smiles, a, b, c, d, angle_deg) {
    let deferred3_0;
    let deferred3_1;
    try {
        const ptr0 = passStringToWasm0(smiles, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.set_dihedral_json(ptr0, len0, a, b, c, d, angle_deg);
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
 * 3D shape descriptors as a JSON object.
 *
 * Keys: `pmi1`, `pmi2`, `pmi3`, `npr1`, `npr2`, `asphericity`, `eccentricity`,
 * `radiusOfGyration`, `planeOfBestFit`.  Non-finite values (e.g. single-atom
 * molecules where pmi3 = 0) are serialised as JSON `null`.
 * @param {MolHandle} mol
 * @returns {string}
 */
export function shape_descriptors_json(mol) {
    let deferred1_0;
    let deferred1_1;
    try {
        _assertClass(mol, MolHandle);
        const ret = wasm.shape_descriptors_json(mol.__wbg_ptr);
        deferred1_0 = ret[0];
        deferred1_1 = ret[1];
        return getStringFromWasm0(ret[0], ret[1]);
    } finally {
        wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
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
 * Like `smarts_match_atoms` but with explicit chirality matching control.
 *
 * When `use_chirality=true`, SMARTS chirality primitives `[@]` and `[@@]` are
 * matched against the target molecule's stereochemistry. When `false`, chirality
 * is ignored (RDKit default).
 * @param {string} smarts
 * @param {MolHandle} mol
 * @param {boolean} use_chirality
 * @returns {string}
 */
export function smarts_match_atoms_with_chirality(smarts, mol, use_chirality) {
    let deferred3_0;
    let deferred3_1;
    try {
        const ptr0 = passStringToWasm0(smarts, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        _assertClass(mol, MolHandle);
        const ret = wasm.smarts_match_atoms_with_chirality(ptr0, len0, mol.__wbg_ptr, use_chirality);
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
 * Serialise a JSON array of SMILES to an SDF string.
 *
 * Generates 2D coordinates for each molecule.  Property data can be
 * included by using `sdf_from_records_json` instead.
 * @param {string} smiles_json
 * @returns {string}
 */
export function smiles_array_to_sdf(smiles_json) {
    let deferred3_0;
    let deferred3_1;
    try {
        const ptr0 = passStringToWasm0(smiles_json, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.smiles_array_to_sdf(ptr0, len0);
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
 * Convert a SMILES to a minimal Tripos MOL2 string (no 3D coordinates).
 *
 * Returns `"error:<msg>"` on parse failure.
 * @param {string} smiles
 * @returns {string}
 */
export function smiles_to_mol2(smiles) {
    let deferred2_0;
    let deferred2_1;
    try {
        const ptr0 = passStringToWasm0(smiles, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.smiles_to_mol2(ptr0, len0);
        deferred2_0 = ret[0];
        deferred2_1 = ret[1];
        return getStringFromWasm0(ret[0], ret[1]);
    } finally {
        wasm.__wbindgen_free(deferred2_0, deferred2_1, 1);
    }
}

/**
 * Write a molecule to AutoDock PDBQT format.
 *
 * `coords_json` — JSON array of `[x,y,z]` arrays (Å). Pass `"[]"` for zero coords.
 * `charges_json` — JSON array of partial charges. Pass `"[]"` to write zeros.
 * `name` — ligand name for the REMARK header.
 *
 * Returns the PDBQT string, or `"error:<msg>"` on failure.
 * @param {string} smiles
 * @param {string} coords_json
 * @param {string} charges_json
 * @param {string} name
 * @returns {string}
 */
export function smiles_to_pdbqt(smiles, coords_json, charges_json, name) {
    let deferred5_0;
    let deferred5_1;
    try {
        const ptr0 = passStringToWasm0(smiles, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ptr1 = passStringToWasm0(coords_json, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len1 = WASM_VECTOR_LEN;
        const ptr2 = passStringToWasm0(charges_json, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len2 = WASM_VECTOR_LEN;
        const ptr3 = passStringToWasm0(name, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len3 = WASM_VECTOR_LEN;
        const ret = wasm.smiles_to_pdbqt(ptr0, len0, ptr1, len1, ptr2, len2, ptr3, len3);
        deferred5_0 = ret[0];
        deferred5_1 = ret[1];
        return getStringFromWasm0(ret[0], ret[1]);
    } finally {
        wasm.__wbindgen_free(deferred5_0, deferred5_1, 1);
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

/**
 * Smallest Set of Smallest Rings (SSSR) as a JSON array of atom-index arrays.
 *
 * Example return value for naphthalene:
 * `[[0,1,2,3,4,5],[5,6,7,8,9,4]]`
 * @param {MolHandle} mol
 * @returns {string}
 */
export function sssr_rings_json(mol) {
    let deferred1_0;
    let deferred1_1;
    try {
        _assertClass(mol, MolHandle);
        const ret = wasm.sssr_rings_json(mol.__wbg_ptr);
        deferred1_0 = ret[0];
        deferred1_1 = ret[1];
        return getStringFromWasm0(ret[0], ret[1]);
    } finally {
        wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
    }
}

/**
 * Standardize a SMILES string and return the canonical SMILES of the result.
 *
 * Applies: largest fragment extraction → charge neutralization.
 * Returns `"error:<msg>"` on parse failure.
 * @param {string} smiles
 * @returns {string}
 */
export function standardize_smiles(smiles) {
    let deferred2_0;
    let deferred2_1;
    try {
        const ptr0 = passStringToWasm0(smiles, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.standardize_smiles(ptr0, len0);
        deferred2_0 = ret[0];
        deferred2_1 = ret[1];
        return getStringFromWasm0(ret[0], ret[1]);
    } finally {
        wasm.__wbindgen_free(deferred2_0, deferred2_1, 1);
    }
}

/**
 * Standardize a SMILES string and return result SMILES plus an audit report as JSON.
 *
 * Boolean flags map directly to `StandardizeOptions`.
 * Returns `"error:<msg>"` on parse or serialization failure.
 * @param {string} smiles
 * @param {boolean} largest_fragment_only
 * @param {boolean} neutralize_charges
 * @param {boolean} remove_explicit_h
 * @param {boolean} canonical_tautomer
 * @returns {string}
 */
export function standardize_smiles_report_json(smiles, largest_fragment_only, neutralize_charges, remove_explicit_h, canonical_tautomer) {
    let deferred2_0;
    let deferred2_1;
    try {
        const ptr0 = passStringToWasm0(smiles, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.standardize_smiles_report_json(ptr0, len0, largest_fragment_only, neutralize_charges, remove_explicit_h, canonical_tautomer);
        deferred2_0 = ret[0];
        deferred2_1 = ret[1];
        return getStringFromWasm0(ret[0], ret[1]);
    } finally {
        wasm.__wbindgen_free(deferred2_0, deferred2_1, 1);
    }
}

export function start() {
    wasm.start();
}

/**
 * Compute the stereo Parent and return a status-shaped JSON result.
 * @param {MolHandle} mol
 * @returns {string}
 */
export function stereo_parent_json(mol) {
    let deferred1_0;
    let deferred1_1;
    try {
        _assertClass(mol, MolHandle);
        const ret = wasm.stereo_parent_json(mol.__wbg_ptr);
        deferred1_0 = ret[0];
        deferred1_1 = ret[1];
        return getStringFromWasm0(ret[0], ret[1]);
    } finally {
        wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
    }
}

/**
 * Compute the composed Super Parent with explicit resource limits.
 * @param {MolHandle} mol
 * @param {number} max_transforms
 * @param {number} max_tautomers
 * @param {bigint | null} [timeout_ms]
 * @returns {string}
 */
export function super_parent_json(mol, max_transforms, max_tautomers, timeout_ms) {
    let deferred1_0;
    let deferred1_1;
    try {
        _assertClass(mol, MolHandle);
        const ret = wasm.super_parent_json(mol.__wbg_ptr, max_transforms, max_tautomers, !isLikeNone(timeout_ms), isLikeNone(timeout_ms) ? BigInt(0) : timeout_ms);
        deferred1_0 = ret[0];
        deferred1_1 = ret[1];
        return getStringFromWasm0(ret[0], ret[1]);
    } finally {
        wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
    }
}

/**
 * Compute the composed Super Parent and expose every ordered stage.
 * @param {MolHandle} mol
 * @param {number} max_transforms
 * @param {number} max_tautomers
 * @param {bigint | null} [timeout_ms]
 * @returns {string}
 */
export function super_parent_report_json(mol, max_transforms, max_tautomers, timeout_ms) {
    let deferred1_0;
    let deferred1_1;
    try {
        _assertClass(mol, MolHandle);
        const ret = wasm.super_parent_report_json(mol.__wbg_ptr, max_transforms, max_tautomers, !isLikeNone(timeout_ms), isLikeNone(timeout_ms) ? BigInt(0) : timeout_ms);
        deferred1_0 = ret[0];
        deferred1_1 = ret[1];
        return getStringFromWasm0(ret[0], ret[1]);
    } finally {
        wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
    }
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
 * Tanimoto similarity between `a` and `b` using ECFP6 fingerprints.
 * @param {MolHandle} a
 * @param {MolHandle} b
 * @returns {number}
 */
export function tanimoto_ecfp6(a, b) {
    _assertClass(a, MolHandle);
    _assertClass(b, MolHandle);
    const ret = wasm.tanimoto_ecfp6(a.__wbg_ptr, b.__wbg_ptr);
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
 * Tanimoto similarity between `a` and `b` using FCFP6 (radius-3 pharmacophore) fingerprints.
 * @param {MolHandle} a
 * @param {MolHandle} b
 * @returns {number}
 */
export function tanimoto_fcfp6(a, b) {
    _assertClass(a, MolHandle);
    _assertClass(b, MolHandle);
    const ret = wasm.tanimoto_fcfp6(a.__wbg_ptr, b.__wbg_ptr);
    return ret;
}

/**
 * Tanimoto similarity between `a` and `b` using MACCS 166-bit fingerprints.
 * @param {MolHandle} a
 * @param {MolHandle} b
 * @returns {number}
 */
export function tanimoto_maccs(a, b) {
    _assertClass(a, MolHandle);
    _assertClass(b, MolHandle);
    const ret = wasm.tanimoto_maccs(a.__wbg_ptr, b.__wbg_ptr);
    return ret;
}

/**
 * Tanimoto-like similarity between two SMILES via MHFP (MinHash Jaccard approximation).
 * @param {string} smi1
 * @param {string} smi2
 * @returns {number}
 */
export function tanimoto_mhfp_smiles(smi1, smi2) {
    const ptr0 = passStringToWasm0(smi1, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
    const len0 = WASM_VECTOR_LEN;
    const ptr1 = passStringToWasm0(smi2, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
    const len1 = WASM_VECTOR_LEN;
    const ret = wasm.tanimoto_mhfp_smiles(ptr0, len0, ptr1, len1);
    if (ret[2]) {
        throw takeFromExternrefTable0(ret[1]);
    }
    return ret[0];
}

/**
 * Compute ECFP4 Tanimoto similarity from one query SMILES to all db SMILES (dense output).
 *
 * `db_smiles_json`: JSON array of SMILES strings (max 1024 via WASM_MAX_BATCH_ITEMS).
 *
 * Returns a flat JSON array of f32 scores, one per db entry, e.g. `[0.12,0.0,0.85]`.
 * No zero-filtering: the length always equals the number of db entries.
 * Returns `"error:<msg>"` on parse failure or oversized input.
 * @param {string} query_smi
 * @param {string} db_smiles_json
 * @returns {string}
 */
export function tanimoto_row_json(query_smi, db_smiles_json) {
    let deferred3_0;
    let deferred3_1;
    try {
        const ptr0 = passStringToWasm0(query_smi, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ptr1 = passStringToWasm0(db_smiles_json, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len1 = WASM_VECTOR_LEN;
        const ret = wasm.tanimoto_row_json(ptr0, len0, ptr1, len1);
        deferred3_0 = ret[0];
        deferred3_1 = ret[1];
        return getStringFromWasm0(ret[0], ret[1]);
    } finally {
        wasm.__wbindgen_free(deferred3_0, deferred3_1, 1);
    }
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
 * Compute the tautomer parent with explicit resource limits.
 * Returns `{"smiles":"...","status":"completed"}` (or a structured
 * error) so callers can distinguish a definite result from a budget-limited
 * one.
 * @param {MolHandle} mol
 * @param {number} max_transforms
 * @param {number} max_tautomers
 * @param {bigint | null} [timeout_ms]
 * @returns {string}
 */
export function tautomer_parent_json(mol, max_transforms, max_tautomers, timeout_ms) {
    let deferred1_0;
    let deferred1_1;
    try {
        _assertClass(mol, MolHandle);
        const ret = wasm.tautomer_parent_json(mol.__wbg_ptr, max_transforms, max_tautomers, !isLikeNone(timeout_ms), isLikeNone(timeout_ms) ? BigInt(0) : timeout_ms);
        deferred1_0 = ret[0];
        deferred1_1 = ret[1];
        return getStringFromWasm0(ret[0], ret[1]);
    } finally {
        wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
    }
}

/**
 * Serialise a `MolHandle` to a CML string with 2D coordinates.
 *
 * Coordinates are generated using the same 2D layout engine as `to_mol_block`.
 * @param {MolHandle} mol
 * @returns {string}
 */
export function to_cml(mol) {
    let deferred1_0;
    let deferred1_1;
    try {
        _assertClass(mol, MolHandle);
        const ret = wasm.to_cml(mol.__wbg_ptr);
        deferred1_0 = ret[0];
        deferred1_1 = ret[1];
        return getStringFromWasm0(ret[0], ret[1]);
    } finally {
        wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
    }
}

/**
 * Write a molecule + coordinates as an Extended XYZ (extxyz) frame.
 *
 * `coords_json`: `[[x,y,z],...]` (Å), same order and length as `mol`'s
 * atoms.
 *
 * `options_json`: an optional JSON object,
 * `{"lattice":[9 numbers]|null,"properties":{"name":[[...]],...},"info":{"key":"value",...}}`
 * -- pass `"{}"` for a plain (non-extended) frame. Only real-valued
 * per-atom `properties` columns are supported from this binding (matches
 * the Python `to_extxyz` binding's scope); build a
 * `chematic_mol::XyzProperty` directly from Rust for integer/string/logical
 * columns.
 *
 * Returns a JS error if `coords_json`/`options_json` are malformed, if
 * `coords_json`'s length doesn't match `mol`'s atom count, or if a
 * `properties` column's row count doesn't match it.
 * @param {MolHandle} mol
 * @param {string} coords_json
 * @param {string} options_json
 * @returns {string}
 */
export function to_extxyz_json(mol, coords_json, options_json) {
    let deferred4_0;
    let deferred4_1;
    try {
        _assertClass(mol, MolHandle);
        const ptr0 = passStringToWasm0(coords_json, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ptr1 = passStringToWasm0(options_json, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len1 = WASM_VECTOR_LEN;
        const ret = wasm.to_extxyz_json(mol.__wbg_ptr, ptr0, len0, ptr1, len1);
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
 * Serialize a molecule to a MOL V2000 block with 2D coordinates.
 *
 * Atom positions are computed via the same layout engine used for SVG depiction
 * and converted to Ångström units (`1.5 Å` per bond).
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

/**
 * Serialise a `MolHandle` to MOL V3000 format with 2D coordinates.
 * @param {MolHandle} mol
 * @returns {string}
 */
export function to_mol_v3000_block(mol) {
    let deferred1_0;
    let deferred1_1;
    try {
        _assertClass(mol, MolHandle);
        const ret = wasm.to_mol_v3000_block(mol.__wbg_ptr);
        deferred1_0 = ret[0];
        deferred1_1 = ret[1];
        return getStringFromWasm0(ret[0], ret[1]);
    } finally {
        wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
    }
}

/**
 * Serialise a `MolHandle` to a MolJSON string (pretty-printed).
 *
 * Atom IDs are assigned as `"a1"`, `"a2"`, … in molecule atom order.
 * The `hydrogens` field reflects computed implicit H count.
 * @param {MolHandle} mol
 * @returns {string}
 */
export function to_moljson(mol) {
    let deferred1_0;
    let deferred1_1;
    try {
        _assertClass(mol, MolHandle);
        const ret = wasm.to_moljson(mol.__wbg_ptr);
        deferred1_0 = ret[0];
        deferred1_1 = ret[1];
        return getStringFromWasm0(ret[0], ret[1]);
    } finally {
        wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
    }
}

/**
 * Serialize a `MolHandle` + coordinates (Å) + molecular charge/multiplicity
 * as a QCSchema `qcschema_molecule` JSON document (coordinates converted
 * to Bohr).
 *
 * `coords_json`: `[[x,y,z],...]` (Å), same order and length as `mol`'s
 * atoms.
 * @param {MolHandle} mol
 * @param {string} coords_json
 * @param {number} charge
 * @param {bigint} multiplicity
 * @returns {string}
 */
export function to_qcschema_molecule_json(mol, coords_json, charge, multiplicity) {
    let deferred3_0;
    let deferred3_1;
    try {
        _assertClass(mol, MolHandle);
        const ptr0 = passStringToWasm0(coords_json, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.to_qcschema_molecule_json(mol.__wbg_ptr, ptr0, len0, charge, multiplicity);
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
 * Serialize a molecule to XYZ format.
 *
 * 3D coordinates are generated via distance-geometry placement.
 * @param {MolHandle} mol
 * @returns {string}
 */
export function to_xyz(mol) {
    let deferred1_0;
    let deferred1_1;
    try {
        _assertClass(mol, MolHandle);
        const ret = wasm.to_xyz(mol.__wbg_ptr);
        deferred1_0 = ret[0];
        deferred1_1 = ret[1];
        return getStringFromWasm0(ret[0], ret[1]);
    } finally {
        wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
    }
}

/**
 * Torsion fingerprint as a bit-packed byte vector (256 bytes = 2048 bits).
 * @param {MolHandle} mol
 * @returns {Uint8Array}
 */
export function torsion_bitvec(mol) {
    _assertClass(mol, MolHandle);
    const ret = wasm.torsion_bitvec(mol.__wbg_ptr);
    var v1 = getArrayU8FromWasm0(ret[0], ret[1]).slice();
    wasm.__wbindgen_free(ret[0], ret[1] * 1, 1);
    return v1;
}

/**
 * Virtual screen a query SMILES against a database of SMILES using ECFP4 Tanimoto.
 *
 * `db_smiles_json`: JSON array of SMILES strings (max 1024 via WASM_MAX_BATCH_ITEMS).
 * `k`: number of top hits to return; clamped to db size if larger.
 *
 * Returns JSON: `{"results":[{"rank":1,"score":0.85,"smiles":"CCO","idx":42},...]}`.
 * Returns `"error:<msg>"` on any parse failure or oversized input.
 * @param {string} query_smi
 * @param {string} db_smiles_json
 * @param {number} k
 * @returns {string}
 */
export function virtual_screen_ecfp4_json(query_smi, db_smiles_json, k) {
    let deferred3_0;
    let deferred3_1;
    try {
        const ptr0 = passStringToWasm0(query_smi, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ptr1 = passStringToWasm0(db_smiles_json, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len1 = WASM_VECTOR_LEN;
        const ret = wasm.virtual_screen_ecfp4_json(ptr0, len0, ptr1, len1, k);
        deferred3_0 = ret[0];
        deferred3_1 = ret[1];
        return getStringFromWasm0(ret[0], ret[1]);
    } finally {
        wasm.__wbindgen_free(deferred3_0, deferred3_1, 1);
    }
}

/**
 * Compute WHIM descriptors (Weighted Holistic Invariant Molecular) from 3D coordinates.
 * Returns JSON array of 22 values: 11 unit-weight descriptors followed by 11 mass-weight
 * descriptors. Each 11-element block is [λ₁, λ₂, λ₃, ν₁, ν₂, ν₃, T, A, V, K, D].
 * @param {MolHandle} mol
 * @returns {string}
 */
export function whim_descriptors_json(mol) {
    let deferred1_0;
    let deferred1_1;
    try {
        _assertClass(mol, MolHandle);
        const ret = wasm.whim_descriptors_json(mol.__wbg_ptr);
        deferred1_0 = ret[0];
        deferred1_1 = ret[1];
        return getStringFromWasm0(ret[0], ret[1]);
    } finally {
        wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
    }
}

/**
 * Compute combined WHIM + GETAWAY descriptors (**41** values total) as JSON array.
 *
 * Returns WHIM[0..21] (22 values) followed by GETAWAY[0..18] (19 values) = 41 total.
 * Useful for ML pipelines requiring both shape and topologic features.
 * @param {MolHandle} mol
 * @returns {string}
 */
export function whim_getaway_combined_json(mol) {
    let deferred1_0;
    let deferred1_1;
    try {
        _assertClass(mol, MolHandle);
        const ret = wasm.whim_getaway_combined_json(mol.__wbg_ptr);
        deferred1_0 = ret[0];
        deferred1_1 = ret[1];
        return getStringFromWasm0(ret[0], ret[1]);
    } finally {
        wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
    }
}

/**
 * Write a grid (in the JSON shape [`cube_grid_json`] returns) as a
 * Gaussian Cube file.
 * @param {string} grid_json
 * @returns {string}
 */
export function write_cube_json(grid_json) {
    let deferred3_0;
    let deferred3_1;
    try {
        const ptr0 = passStringToWasm0(grid_json, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.write_cube_json(ptr0, len0);
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
 * Write a LAMMPS data file from the JSON shape [`lammps_data_to_json`]
 * returns.
 * @param {string} json
 * @returns {string}
 */
export function write_lammps_data_json(json) {
    let deferred3_0;
    let deferred3_1;
    try {
        const ptr0 = passStringToWasm0(json, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.write_lammps_data_json(ptr0, len0);
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
 * Write a single LAMMPS dump frame from the JSON shape
 * [`lammps_dump_frame_to_json_str`] returns.
 * @param {string} json
 * @returns {string}
 */
export function write_lammps_dump_frame_json(json) {
    let deferred3_0;
    let deferred3_1;
    try {
        const ptr0 = passStringToWasm0(json, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.write_lammps_dump_frame_json(ptr0, len0);
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
 * Write a LAMMPS trajectory (N frames concatenated back to back, matching
 * [`chematic_mol::write_lammps_trajectory`]) from a JSON array of frames
 * in the shape [`lammps_dump_frame_to_json_str`] returns.
 * @param {string} json
 * @returns {string}
 */
export function write_lammps_trajectory_json(json) {
    let deferred3_0;
    let deferred3_1;
    try {
        const ptr0 = passStringToWasm0(json, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.write_lammps_trajectory_json(ptr0, len0);
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
 * Write an mmCIF file from atom records in the JSON shape
 * [`mmcif_to_json`]'s `"atoms"` array uses (a full record per atom, not
 * just element+coordinates -- mmCIF has no equivalent of "build from a
 * bare `MolHandle`", since occupancy/B-factor/chain/residue fields have no
 * source in a plain [`MolHandle`]).
 *
 * `cell_json`: `"null"` or `{"a":...,"b":...,"c":...,"alpha":...,"beta":...,"gamma":...}`.
 * `space_group`: pass `""` for none.
 * @param {string} records_json
 * @param {string} cell_json
 * @param {string} space_group
 * @param {string} data_block_name
 * @returns {string}
 */
export function write_mmcif_json(records_json, cell_json, space_group, data_block_name) {
    let deferred6_0;
    let deferred6_1;
    try {
        const ptr0 = passStringToWasm0(records_json, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ptr1 = passStringToWasm0(cell_json, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len1 = WASM_VECTOR_LEN;
        const ptr2 = passStringToWasm0(space_group, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len2 = WASM_VECTOR_LEN;
        const ptr3 = passStringToWasm0(data_block_name, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len3 = WASM_VECTOR_LEN;
        const ret = wasm.write_mmcif_json(ptr0, len0, ptr1, len1, ptr2, len2, ptr3, len3);
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
 * Write a grid as an OpenDX file. Fails closed for a
 * [`chematic_mol::GridUnits::Bohr`]-tagged grid (OpenDX has no unit tag of
 * its own and is universally read back as Ångström -- see
 * `chematic_mol::opendx`'s module docs) and for a grid carrying any atoms
 * (OpenDX has no atom section). Use [`write_opendx_lossy_json`] to opt
 * into an explicit Bohr->Ångström conversion instead of failing.
 * @param {string} grid_json
 * @returns {string}
 */
export function write_opendx_json(grid_json) {
    let deferred3_0;
    let deferred3_1;
    try {
        const ptr0 = passStringToWasm0(grid_json, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.write_opendx_json(ptr0, len0);
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
 * Like [`write_opendx_json`], but a [`chematic_mol::GridUnits::Bohr`] grid
 * has its `origin`/`axes` explicitly converted to Ångström rather than
 * rejected (`values` -- the scalar-field samples themselves -- are never
 * rescaled; see `write_opendx_lossy`'s doc comment). Still fails for a
 * grid carrying any atoms.
 * @param {string} grid_json
 * @returns {string}
 */
export function write_opendx_lossy_json(grid_json) {
    let deferred3_0;
    let deferred3_1;
    try {
        const ptr0 = passStringToWasm0(grid_json, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.write_opendx_lossy_json(ptr0, len0);
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
 * Write an ORCA input file from the JSON shape [`orca_input_to_json`]
 * returns.
 * @param {string} json
 * @returns {string}
 */
export function write_orca_input_json(json) {
    let deferred3_0;
    let deferred3_1;
    try {
        const ptr0 = passStringToWasm0(json, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.write_orca_input_json(ptr0, len0);
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
 * Write a PQR file from atom records in the JSON shape [`pqr_to_json`]'s
 * `"atoms"` array uses. Each atom's `chain_id` independently controls
 * whether that line is written with or without the (optional) chain
 * column.
 * @param {string} records_json
 * @returns {string}
 */
export function write_pqr_json(records_json) {
    let deferred3_0;
    let deferred3_1;
    try {
        const ptr0 = passStringToWasm0(records_json, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
        const len0 = WASM_VECTOR_LEN;
        const ret = wasm.write_pqr_json(ptr0, len0);
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
 * Non-canonical SMILES for `mol`.
 *
 * Unlike `canonical_smiles`, the output depends on the internal atom ordering
 * and is not normalised.  Useful when round-trip fidelity (preserving atom
 * order) matters more than a canonical form.
 * @param {MolHandle} mol
 * @returns {string}
 */
export function write_smiles(mol) {
    let deferred1_0;
    let deferred1_1;
    try {
        _assertClass(mol, MolHandle);
        const ret = wasm.write_smiles(mol.__wbg_ptr);
        deferred1_0 = ret[0];
        deferred1_1 = ret[1];
        return getStringFromWasm0(ret[0], ret[1]);
    } finally {
        wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
    }
}

/**
 * XLogP3 partition coefficient (alternative to Crippen LogP).
 * Returns JSON: `{"xlogp3": float}`.
 * @param {MolHandle} mol
 * @returns {string}
 */
export function xlogp3_json(mol) {
    let deferred1_0;
    let deferred1_1;
    try {
        _assertClass(mol, MolHandle);
        const ret = wasm.xlogp3_json(mol.__wbg_ptr);
        deferred1_0 = ret[0];
        deferred1_1 = ret[1];
        return getStringFromWasm0(ret[0], ret[1]);
    } finally {
        wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
    }
}

/**
 * Per-atom XLogP3 contributions.
 * Returns JSON array of floats (one per heavy atom).
 * @param {MolHandle} mol
 * @returns {string}
 */
export function xlogp3_per_atom_json(mol) {
    let deferred1_0;
    let deferred1_1;
    try {
        _assertClass(mol, MolHandle);
        const ret = wasm.xlogp3_per_atom_json(mol.__wbg_ptr);
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
        __wbg___wbindgen_is_undefined_721f8decd50c87a3: function(arg0) {
            const ret = arg0 === undefined;
            return ret;
        },
        __wbg___wbindgen_string_get_71bb4348194e31f0: function(arg0, arg1) {
            const obj = arg1;
            const ret = typeof(obj) === 'string' ? obj : undefined;
            var ptr1 = isLikeNone(ret) ? 0 : passStringToWasm0(ret, wasm.__wbindgen_malloc, wasm.__wbindgen_realloc);
            var len1 = WASM_VECTOR_LEN;
            getDataViewMemory0().setInt32(arg0 + 4 * 1, len1, true);
            getDataViewMemory0().setInt32(arg0 + 4 * 0, ptr1, true);
        },
        __wbg___wbindgen_throw_ea4887a5f8f9a9db: function(arg0, arg1) {
            throw new Error(getStringFromWasm0(arg0, arg1));
        },
        __wbg_error_933f449d72fef598: function(arg0) {
            console.error(arg0);
        },
        __wbg_new_from_slice_3f5658f83d8d0725: function(arg0, arg1) {
            const ret = new Uint32Array(getArrayU32FromWasm0(arg0, arg1));
            return ret;
        },
        __wbg_new_from_slice_98e57cb2fe2e6a5d: function(arg0, arg1) {
            const ret = new Float64Array(getArrayF64FromWasm0(arg0, arg1));
            return ret;
        },
        __wbg_now_e7c6795a7f81e10f: function(arg0) {
            const ret = arg0.now();
            return ret;
        },
        __wbg_performance_3fcf6e32a7e1ed0a: function(arg0) {
            const ret = arg0.performance;
            return ret;
        },
        __wbg_static_accessor_GLOBAL_THIS_2fee5048bcca5938: function() {
            const ret = typeof globalThis === 'undefined' ? null : globalThis;
            return isLikeNone(ret) ? 0 : addToExternrefTable0(ret);
        },
        __wbg_static_accessor_GLOBAL_ce44e66a4935da8c: function() {
            const ret = typeof global === 'undefined' ? null : global;
            return isLikeNone(ret) ? 0 : addToExternrefTable0(ret);
        },
        __wbg_static_accessor_SELF_44f6e0cb5e67cdad: function() {
            const ret = typeof self === 'undefined' ? null : self;
            return isLikeNone(ret) ? 0 : addToExternrefTable0(ret);
        },
        __wbg_static_accessor_WINDOW_168f178805d978fe: function() {
            const ret = typeof window === 'undefined' ? null : window;
            return isLikeNone(ret) ? 0 : addToExternrefTable0(ret);
        },
        __wbindgen_cast_0000000000000001: function(arg0) {
            // Cast intrinsic for `F64 -> Externref`.
            const ret = arg0;
            return ret;
        },
        __wbindgen_cast_0000000000000002: function(arg0, arg1) {
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

const ConformerHandleFinalization = (typeof FinalizationRegistry === 'undefined')
    ? { register: () => {}, unregister: () => {} }
    : new FinalizationRegistry(ptr => wasm.__wbg_conformerhandle_free(ptr, 1));
const DepictOptionsFinalization = (typeof FinalizationRegistry === 'undefined')
    ? { register: () => {}, unregister: () => {} }
    : new FinalizationRegistry(ptr => wasm.__wbg_depictoptions_free(ptr, 1));
const MhfpLshHandleFinalization = (typeof FinalizationRegistry === 'undefined')
    ? { register: () => {}, unregister: () => {} }
    : new FinalizationRegistry(ptr => wasm.__wbg_mhfplshhandle_free(ptr, 1));
const MolHandleFinalization = (typeof FinalizationRegistry === 'undefined')
    ? { register: () => {}, unregister: () => {} }
    : new FinalizationRegistry(ptr => wasm.__wbg_molhandle_free(ptr, 1));

function addToExternrefTable0(obj) {
    const idx = wasm.__externref_table_alloc();
    wasm.__wbindgen_externrefs.set(idx, obj);
    return idx;
}

function _assertClass(instance, klass) {
    if (!(instance instanceof klass)) {
        throw new Error(`expected instance of ${klass.name}`);
    }
}

function getArrayF64FromWasm0(ptr, len) {
    ptr = ptr >>> 0;
    return getFloat64ArrayMemory0().subarray(ptr / 8, ptr / 8 + len);
}

function getArrayU32FromWasm0(ptr, len) {
    ptr = ptr >>> 0;
    return getUint32ArrayMemory0().subarray(ptr / 4, ptr / 4 + len);
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

let cachedFloat64ArrayMemory0 = null;
function getFloat64ArrayMemory0() {
    if (cachedFloat64ArrayMemory0 === null || cachedFloat64ArrayMemory0.byteLength === 0) {
        cachedFloat64ArrayMemory0 = new Float64Array(wasm.memory.buffer);
    }
    return cachedFloat64ArrayMemory0;
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

function isLikeNone(x) {
    return x === undefined || x === null;
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
    cachedFloat64ArrayMemory0 = null;
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
