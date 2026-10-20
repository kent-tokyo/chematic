import { T, STATIC_IDS } from './translations.js';

// ============================================================
// i18n
// ============================================================


let lang = 'en';
function t(key) { return (T[lang] && T[lang][key]) || T.en[key] || key; }



function applyLang(newLang) {
  lang = newLang;
  document.documentElement.lang = lang === 'ja' ? 'ja' : lang === 'zh' ? 'zh-Hans' : 'en';
  localStorage.setItem('chematic-lang', lang);

  document.querySelectorAll('[data-lang]').forEach(btn => {
    btn.classList.toggle('active', btn.dataset.lang === lang);
  });

  for (const [id, key] of Object.entries(STATIC_IDS)) {
    const el = document.getElementById(id);
    if (el) el.textContent = t(key);
  }
  document.getElementById('loading-text').textContent = t('app.loading');
  document.getElementById('smiles-input').placeholder = t('desc.placeholder');
  document.getElementById('sim-a').placeholder = t('desc.placeholder');
  document.getElementById('sim-b').placeholder = t('desc.placeholder');
  document.getElementById('sdf-input').placeholder = t('sdf.placeholder');
  document.getElementById('report-smiles').placeholder = t('report.placeholder');
  document.getElementById('h-report-title').textContent = t('report.title');
  const btnReport = document.getElementById('btn-report');
  if (btnReport) btnReport.textContent = t('report.button');

  renderRuleCards(lastMol);
  renderSimTable(lastSimResults);
  renderGallery();
}

// ============================================================
// Tab switching
// ============================================================
let activeTab = 'tab-2d';
let activeCategory = 'cat-structure';
let viewer3dInited = false;
let viewer3dMolDirty = false;

function switchTab(tabId) {
  activeTab = tabId;
  document.querySelectorAll('.tab-btn').forEach(btn => {
    btn.classList.toggle('active', btn.dataset.tab === tabId);
  });
  document.querySelectorAll('.tab-panel').forEach(panel => {
    panel.classList.toggle('active', panel.id === tabId);
  });

  // Keep the category row in sync whenever a tab is opened directly
  // (the gallery's hardcoded switchTab('tab-2d') call, page load, etc.)
  // rather than only through a category click.
  const btnEl = document.querySelector(`.tab-bar-inner .tab-btn[data-tab="${tabId}"]`);
  const cat = btnEl ? btnEl.dataset.cat : null;
  if (cat && cat !== activeCategory) {
    activeCategory = cat;
    document.querySelectorAll('.tab-cat-btn').forEach(b => {
      b.classList.toggle('active', b.dataset.cat === cat);
    });
    document.querySelectorAll('.tab-bar-inner .tab-btn').forEach(b => {
      b.style.display = b.dataset.cat === cat ? '' : 'none';
    });
  }

  if (tabId === 'tab-3d' && lastMol) {
    if (!viewer3dInited) {
      init3dViewer();
    } else if (viewer3dMolDirty) {
      render3d(lastMol);
    }
  }
}

function setCategory(catId) {
  activeCategory = catId;
  document.querySelectorAll('.tab-cat-btn').forEach(btn => {
    btn.classList.toggle('active', btn.dataset.cat === catId);
  });
  document.querySelectorAll('.tab-bar-inner .tab-btn').forEach(btn => {
    btn.style.display = btn.dataset.cat === catId ? '' : 'none';
  });
  const firstBtn = document.querySelector(`.tab-bar-inner .tab-btn[data-cat="${catId}"]`);
  if (firstBtn) switchTab(firstBtn.dataset.tab);
}

// A small set of task-oriented entry points makes the Playground useful to
// people who do not yet know which chemistry API they need.
function launchWorkflow(workflow) {
  const aspirin = 'CC(=O)Oc1ccccc1C(=O)O';
  const caffeine = 'Cn1cnc2c1c(=O)n(C)c(=O)n2C';
  switch (workflow) {
    case 'molecule':
      setCategory('cat-structure');
      document.getElementById('smiles-input').value = aspirin;
      calculate();
      break;
    case 'similarity':
      setCategory('cat-analysis');
      switchTab('tab-sim');
      document.getElementById('sim-a').value = aspirin;
      document.getElementById('sim-b').value = caffeine;
      compare();
      break;
    case 'report':
      setCategory('cat-analysis');
      switchTab('tab-report');
      document.getElementById('report-smiles').value = aspirin;
      document.getElementById('btn-report').click();
      break;
    case 'reaction':
      setCategory('cat-reactions');
      switchTab('tab-rxn');
      document.getElementById('rxn-esterify')?.click();
      break;
    case 'formats':
      setCategory('cat-data');
      switchTab('tab-formats');
      document.getElementById('btn-formats-example')?.click();
      break;
  }
}

function shareCurrentMolecule() {
  const smiles = document.getElementById('smiles-input').value.trim();
  if (!smiles) return;
  const url = `${location.origin}${location.pathname}#smiles=${encodeURIComponent(smiles)}`;
  history.replaceState(null, '', url);
  const status = document.getElementById('share-molecule-status');
  const done = () => {
    status.textContent = t('share.copied');
    setTimeout(() => { status.textContent = ''; }, 2200);
  };
  if (navigator.clipboard?.writeText) {
    navigator.clipboard.writeText(url).then(done).catch(() => { status.textContent = url; });
  } else {
    status.textContent = url;
  }
}

function sharedMoleculeFromUrl() {
  const prefix = '#smiles=';
  if (!location.hash.startsWith(prefix)) return null;
  try {
    const smiles = decodeURIComponent(location.hash.slice(prefix.length));
    return smiles.length > 0 && smiles.length <= 20000 ? smiles : null;
  } catch (_) { return null; }
}

// ============================================================
// State
// ============================================================
let wasmReady = false;
let parseSmilesWasm, tanimotoEcfp4, tanimotoFcfp4, tanimotoAtomPair, tanimotoTorsion, bricsFrag;
let enumerateStereoIsomers, mmpPairsJson, rgroupDecomposeJson, depictSvgOpts;
let DepictOptions, runReactants, depictSvgGrid, smartsMatchAtoms, generate3dPdb, depictReactionSvg;
let molFromSdfBlock, sdfToSmilesJson, estateIndicesJson, tanimotoTopoPath;
let identifyFunctionalGroups, gasteigerChargesJson, saScore, slogpVsaJson, smrVsaJson, peoeVsaJson;
let detectFunctionalGroups, getAtomInfo, getBondBetween, molBlockFromSmiles;
let moleculeReportJson, compareMultipleJson, screenSmilesJson, generate3dOptimizedPdb;
// Formats & Materials tab
let molFromCube, cubeGridJson, cubeValuesF64, cubeShapeU32;
let opendxGridJson, opendxValuesF64, opendxShapeU32;
let molFromMmcif, mmcifToJson;
let molFromPqr, pqrToJson;
let molFromQcschemaMolecule, qcschemaMoleculeCoordsJson;
let molFromOrcaInput, orcaInputToJson, orcaOutputToJson;
let lammpsDataToJson, lammpsTrajectoryToJson;
let lastMol = null;
let lastSimResults = null;
let lastReportJson = null;
let highlightedAtoms = [];
let highlightColor = '#4f8ef7';
let showIndices = false;
let namedGroupsData = [];
let activeGroupName = null;
let groupColorMap = {};

// ============================================================
// WASM init
// ============================================================
// Keep the shipped WASM assets cache-busted when the binary changes. Without
// this, a browser can retain an older wasm-bindgen JS/wasm pair after a Pages
// deployment and display stale chemistry results until its cache expires.
const WASM_ASSET_VERSION = '1.2.0';

// Formats a byte count as a human-readable "X.X MB" / "XXX KB" string.
function formatBytes(n) {
  if (n >= 1024 * 1024) return (n / (1024 * 1024)).toFixed(1) + ' MB';
  return Math.round(n / 1024) + ' KB';
}

async function initWasm() {
  try {
    const mod = await import(`./pkg/chematic_wasm.js?v=${WASM_ASSET_VERSION}`);
    // Fetch the WASM binary ourselves (rather than letting init() do it
    // internally) so we can measure its real, current size and display it
    // -- the actual bytes served, never a hardcoded number that goes stale
    // across releases. init() accepts a Response directly, so this adds no
    // extra network round-trip.
    const wasmResp = await fetch(`./pkg/chematic_wasm_bg.wasm?v=${WASM_ASSET_VERSION}`);
    const wasmBytes = (await wasmResp.clone().arrayBuffer()).byteLength;
    await mod.default(wasmResp);
    // chematic_version() may be absent on an older locally-built ./pkg (e.g.
    // a fresh git clone before CI has rebuilt it) -- guarded the same way
    // depict_svg_opts is a few lines down, so a stale pkg degrades this one
    // badge instead of failing WASM init entirely.
    const versionEl = document.getElementById('version-badge');
    if (versionEl && typeof mod.chematic_version === 'function') {
      versionEl.textContent = 'v' + mod.chematic_version();
      versionEl.style.display = '';
    }
    const wasmSizeEl = document.getElementById('wasm-size-badge');
    if (wasmSizeEl) wasmSizeEl.textContent = formatBytes(wasmBytes);
    parseSmilesWasm  = mod.parse_smiles;
    tanimotoEcfp4    = mod.tanimoto_ecfp4;
    tanimotoFcfp4    = mod.tanimoto_fcfp4;
    tanimotoAtomPair = mod.tanimoto_atom_pair;
    tanimotoTorsion  = mod.tanimoto_torsion;
    bricsFrag        = mod.brics_fragment_count;
    DepictOptions    = mod.DepictOptions;
    runReactants     = mod.run_reactants;
    depictSvgGrid    = mod.depict_svg_grid;
    smartsMatchAtoms = mod.smarts_match_atoms;
    generate3dPdb       = mod.generate_3d_pdb;
    depictReactionSvg   = mod.depict_reaction_svg;
    molFromSdfBlock          = mod.mol_from_sdf_block;
    sdfToSmilesJson          = mod.sdf_to_smiles_json;
    if (new URLSearchParams(location.search).has('browser-smoke')) {
      window.__browserSmoke = window.__browserSmoke || {};
      window.__browserSmoke.sdfToSmilesJson = sdfToSmilesJson;
    }
    estateIndicesJson        = mod.estate_indices_json;
    tanimotoTopoPath         = mod.tanimoto_topo_path;
    identifyFunctionalGroups = mod.identify_functional_groups;
    gasteigerChargesJson     = mod.gasteiger_charges_json;
    saScore                  = mod.sa_score;
    slogpVsaJson             = mod.slogp_vsa_json;
    smrVsaJson               = mod.smr_vsa_json;
    peoeVsaJson              = mod.peoe_vsa_json;
    detectFunctionalGroups   = mod.detect_functional_groups;
    getAtomInfo              = mod.get_atom_info;
    getBondBetween           = mod.get_bond_between;
    molBlockFromSmiles       = mod.mol_block_from_smiles;
    // v0.1.20 new APIs
    enumerateStereoIsomers   = mod.enumerate_stereo_isomers_json;
    mmpPairsJson             = mod.mmp_pairs_json;
    rgroupDecomposeJson      = mod.rgroup_decompose_json;
    depictSvgOpts            = mod.depict_svg_opts ?? mod.depict_svg;
    // Workflow APIs (Phase 3 report & compare)
    moleculeReportJson       = mod.molecule_report_json;
    compareMultipleJson      = mod.compare_molecules_json;
    screenSmilesJson         = mod.screen_smiles_json;
    generate3dOptimizedPdb   = mod.generate_3d_optimized_pdb;
    // Formats & Materials tab (v0.17.0/v0.18.0 materials-science formats)
    molFromCube              = mod.mol_from_cube;
    cubeGridJson             = mod.cube_grid_json;
    cubeValuesF64            = mod.cube_values_f64;
    cubeShapeU32             = mod.cube_shape_u32;
    opendxGridJson           = mod.opendx_grid_json;
    opendxValuesF64          = mod.opendx_values_f64;
    opendxShapeU32           = mod.opendx_shape_u32;
    molFromMmcif             = mod.mol_from_mmcif;
    mmcifToJson              = mod.mmcif_to_json;
    molFromPqr               = mod.mol_from_pqr;
    pqrToJson                = mod.pqr_to_json;
    molFromQcschemaMolecule  = mod.mol_from_qcschema_molecule;
    qcschemaMoleculeCoordsJson = mod.qcschema_molecule_coords_json;
    molFromOrcaInput         = mod.mol_from_orca_input;
    orcaInputToJson          = mod.orca_input_to_json;
    orcaOutputToJson         = mod.orca_output_to_json;
    lammpsDataToJson         = mod.lammps_data_to_json;
    lammpsTrajectoryToJson   = mod.lammps_trajectory_to_json;
    // Store in window namespace for easy access in dynamic handlers
    window.chemWasm = mod;
    wasmReady = true;
  } catch (e) {
    document.getElementById('loading-text').textContent = 'WASM load failed: ' + String(e);
    throw e;
  }
}

// ============================================================
// 3D Viewer (3Dmol.js)
// ============================================================
let viewer3d = null;
let current3dStyle = 'stick';

function init3dViewer() {
  if (viewer3dInited) return;
  viewer3dInited = true;

  const ph = document.getElementById('viewer-placeholder');
  const phText = document.getElementById('viewer-placeholder-text');
  phText.textContent = t('3d.generating');
  const spin = document.createElement('div');
  spin.className = 'spinner';
  ph.prepend(spin);

  if (!generate3dPdb) {
    phText.textContent = 'WASM build missing generate_3d_pdb — please redeploy';
    spin.remove();
    return;
  }

  const script = document.createElement('script');
  script.src = 'https://3dmol.org/build/3Dmol-min.js';
  script.onload = () => {
    // Use rAF so the tab panel is fully laid out before measuring dimensions
    requestAnimationFrame(() => {
      const el = document.getElementById('viewer-3d');
      const w = el.offsetWidth  || 800;
      const h = el.offsetHeight || 480;
      viewer3d = $3Dmol.createViewer(el, {
        backgroundColor: '#080a0f',
        width:  w,
        height: h,
      });
      ph.style.display = 'none';
      if (lastMol) render3d(lastMol);
    });
  };
  script.onerror = () => {
    phText.textContent = '3Dmol.js failed to load. Check network connection.';
    spin.remove();
  };
  document.head.appendChild(script);
}

function render3d(mol) {
  viewer3dMolDirty = false;
  if (!viewer3d || !generate3dPdb) return;
  clearError('error-3d');
  try {
    const pdb = generate3dPdb(mol);
    viewer3d.clear();
    viewer3d.addModel(pdb, 'pdb');
    applyStyleOnly();
    viewer3d.zoomTo();
    viewer3d.render();
  } catch (e) {
    showError('error-3d', '3D render error: ' + String(e));
  }
}

// Apply style without calling render (render3d calls render at the end)
function applyStyleOnly() {
  if (!viewer3d) return;
  if (current3dStyle === 'ball') {
    viewer3d.setStyle({}, {
      stick:  { colorscheme: 'Jmol', radius: 0.12 },
      sphere: { colorscheme: 'Jmol', scale: 0.3 },
    });
  } else if (current3dStyle === 'surface') {
    viewer3d.setStyle({}, { stick: { colorscheme: 'Jmol', radius: 0.1 } });
    viewer3d.addSurface($3Dmol.SurfaceType.VDW, { opacity: 0.7, colorscheme: 'Jmol' });
  } else {
    viewer3d.setStyle({}, { stick: { colorscheme: 'Jmol', radius: 0.18 } });
  }
}

function apply3dStyle() {
  if (!viewer3d) return;
  // Re-render from scratch when style button is clicked
  if (lastMol) render3d(lastMol);
}

// ============================================================
// Descriptor definitions
// ============================================================
const DESC_DEFS = [
  { key: 'desc.formula',     fn: m => m.formula(),                    isSmiles: false },
  { key: 'desc.mw',          fn: m => m.molecular_weight().toFixed(2)              },
  { key: 'desc.exact_mass',  fn: m => m.exact_mass().toFixed(4)                    },
  { key: 'desc.tpsa',        fn: m => m.tpsa().toFixed(2)                          },
  { key: 'desc.logp',        fn: m => m.logp_crippen().toFixed(3)                  },
  { key: 'desc.fsp3',        fn: m => m.fsp3().toFixed(3)                          },
  { key: 'desc.qed',         fn: m => m.qed().toFixed(3)                           },
  { key: 'desc.mr',          fn: m => m.molar_refractivity().toFixed(2)            },
  { key: 'desc.heavy_atoms', fn: m => String(m.heavy_atom_count())                 },
  { key: 'desc.atoms',       fn: m => String(m.atom_count())                       },
  { key: 'desc.bonds',       fn: m => String(m.bond_count())                       },
  { key: 'desc.hbd',         fn: m => String(m.hbd_count())                        },
  { key: 'desc.hba',         fn: m => String(m.hba_count())                        },
  { key: 'desc.rotbonds',    fn: m => String(m.rotatable_bond_count())             },
  { key: 'desc.rings',       fn: m => String(m.ring_count())                       },
  { key: 'desc.arom_rings',  fn: m => String(m.aromatic_ring_count())              },
  { key: 'desc.heteroatoms', fn: m => String(m.num_heteroatoms())                  },
  { key: 'desc.stereo',      fn: m => String(m.num_stereocenters())                },
  { key: 'desc.charge',      fn: m => String(m.formal_charge_sum())                },
  { key: 'desc.brics',       fn: m => String(bricsFrag(m))                         },
  { key: 'desc.canonical',   fn: m => m.canonical_smiles(), isSmiles: true         },
  { key: 'desc.arom_het',    fn: m => String(m.num_aromatic_heterocycles())        },
  { key: 'desc.aliph_het',   fn: m => String(m.num_aliphatic_heterocycles())       },
  { key: 'desc.spiro',       fn: m => String(m.num_spiro_atoms())                  },
  { key: 'desc.bridge',      fn: m => String(m.num_bridgehead_atoms())             },
  { key: 'desc.estate_sum',  fn: m => m.sum_estate().toFixed(3)                   },
  { key: 'desc.estate_max',  fn: m => m.max_estate().toFixed(3)                   },
  { key: 'desc.estate_min',  fn: m => m.min_estate().toFixed(3)                   },
  { key: 'desc.sa_score',    fn: m => saScore ? saScore(m).toFixed(2) : '—'       },
  { key: 'desc.labute_asa',  fn: m => m.labute_asa().toFixed(2)                   },
];

const RULE_DEFS = [
  { nameKey: 'rules.lipinski', criteriaKey: 'rules.lipinski_c', fn: m => m.lipinski_passes() },
  { nameKey: 'rules.veber',    criteriaKey: 'rules.veber_c',    fn: m => m.veber_passes()    },
  { nameKey: 'rules.egan',     criteriaKey: 'rules.egan_c',     fn: m => m.egan_passes()     },
  { nameKey: 'rules.ghose',    criteriaKey: 'rules.ghose_c',    fn: m => m.ghose_passes()    },
  { nameKey: 'rules.reos',     criteriaKey: 'rules.reos_c',     fn: m => m.reos_passes()     },
  { nameKey: 'rules.pains',    criteriaKey: 'rules.pains_c',    fn: m => m.pains_passes()    },
];

const SIM_DEFS = [
  { nameKey: 'sim.ecfp4',    descKey: 'sim.ecfp4_d',    fn: (a,b) => tanimotoEcfp4(a,b)    },
  { nameKey: 'sim.fcfp4',    descKey: 'sim.fcfp4_d',    fn: (a,b) => tanimotoFcfp4(a,b)    },
  { nameKey: 'sim.atompair', descKey: 'sim.atompair_d', fn: (a,b) => tanimotoAtomPair(a,b) },
  { nameKey: 'sim.torsion',  descKey: 'sim.torsion_d',  fn: (a,b) => tanimotoTorsion(a,b)  },
  { nameKey: 'sim.topo_path',descKey: 'sim.topo_path_d',fn: (a,b) => tanimotoTopoPath(a,b) },
];

const GALLERY = [
  { catKey: 'gallery.nsaids', drugs: [
    { name: 'Aspirin',      smiles: 'CC(=O)Oc1ccccc1C(=O)O' },
    { name: 'Ibuprofen',    smiles: 'CC(C)Cc1ccc(cc1)C(C)C(=O)O' },
    { name: 'Naproxen',     smiles: 'CC(C(=O)O)c1ccc2cccc(OC)c2c1' },
    { name: 'Celecoxib',    smiles: 'Cc1ccc(-c2cc(C(F)(F)F)nn2-c2ccc(cc2)S(N)(=O)=O)cc1' },
    { name: 'Diclofenac',   smiles: 'OC(=O)Cc1ccccc1Nc1c(Cl)cccc1Cl' },
  ]},
  { catKey: 'gallery.antibiotics', drugs: [
    { name: 'Amoxicillin',   smiles: 'CC1(C(N2C(S1)C(C2=O)NC(=O)C(c3ccc(cc3)O)N)C(=O)O)C' },
    { name: 'Ciprofloxacin', smiles: 'C1CC1N2C=C(C(=O)C3=CC(=C(C=C23)N4CCNCC4)F)C(=O)O' },
    { name: 'Metronidazole', smiles: 'CC1=NC=C(N1CCO)[N+](=O)[O-]' },
  ]},
  { catKey: 'gallery.cardio', drugs: [
    { name: 'Atorvastatin', smiles: 'CC(C)c1c(C(=O)Nc2ccccc2)c(-c2ccccc2)c(-c2ccc(F)cc2)n1CC[C@@H](O)C[C@@H](O)CC(=O)O' },
    { name: 'Amlodipine',   smiles: 'CCOC(=O)C1=C(COCCN)NC(=C(C1c1ccccc1Cl)C(=O)OCC)C' },
    { name: 'Metoprolol',   smiles: 'CC(C)NCC(O)COc1ccc(CCOC)cc1' },
    { name: 'Warfarin',     smiles: 'CC(=O)CC(c1ccccc1)c1c(O)c2ccccc2oc1=O' },
  ]},
  { catKey: 'gallery.cns', drugs: [
    { name: 'Caffeine',     smiles: 'Cn1cnc2c1c(=O)n(c(=O)n2C)C' },
    { name: 'Fluoxetine',   smiles: 'CNCCC(c1ccccc1)Oc1ccc(cc1)C(F)(F)F' },
    { name: 'Morphine',     smiles: 'CN1CCC23C4C1CC(C2OC5=C(O)C=CC(=C35)C4O)O' },
    { name: 'Sertraline',   smiles: 'CN[C@@H]1CC[C@H](c2ccc(Cl)c(Cl)c2)c2ccccc21' },
  ]},
  { catKey: 'gallery.oncology', drugs: [
    { name: 'Imatinib',     smiles: 'Cc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc1nccc(-c2cccnc2)n1' },
    { name: 'Taxol',        smiles: 'CC1=C2C(C(=O)C3(C(CC4C(C3C(C(=C2OC1=O)C)(C)C)OC(=O)C5=CC=CC=C5)OC(=O)C4)O)(C)C' },
    { name: 'Erlotinib',    smiles: 'C#Cc1cccc(Nc2ncnc3cc(OCCO)c(OCCO)cc23)c1' },
  ]},
  { catKey: 'gallery.nutrients', drugs: [
    { name: 'Glucose',      smiles: 'OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O' },
    { name: 'Cholesterol',  smiles: 'CC(C)CCCC(C)C1CCC2C1(CCC3C2CC=C4C3(CCC(C4)O)C)C' },
    { name: 'Vitamin C',    smiles: 'OC[C@H](O)[C@H]1OC(=O)C(O)=C1O' },
    { name: 'Riboflavin',   smiles: 'Cc1cc2nc3c(=O)[nH]c(=O)nc3n(C[C@H](O)[C@H](O)[C@H](O)CO)c2cc1C' },
  ]},
];

// ============================================================
// Render: radar chart (pure SVG)
// ============================================================
function renderRadar(mol) {
  const wrap = document.getElementById('radar-wrap');
  wrap.textContent = '';
  if (!mol) return;

  const axes = [
    { label: 'MW',       val: mol.molecular_weight() / 500 },
    { label: 'LogP',     val: (mol.logp_crippen() + 5) / 10 },
    { label: 'HBD',      val: mol.hbd_count() / 5 },
    { label: 'HBA',      val: mol.hba_count() / 10 },
    { label: 'TPSA',     val: mol.tpsa() / 140 },
    { label: 'RotBonds', val: mol.rotatable_bond_count() / 10 },
  ];

  const size = 300, cx = 150, cy = 150, r = 95;
  const n = axes.length;
  const NS = 'http://www.w3.org/2000/svg';
  const svg = document.createElementNS(NS, 'svg');
  svg.setAttribute('viewBox', `0 0 ${size} ${size}`);
  svg.setAttribute('aria-label', 'Property radar chart');

  function angle(i) { return (i / n) * 2 * Math.PI - Math.PI / 2; }
  function pt(frac, i) {
    const a = angle(i);
    return [cx + r * frac * Math.cos(a), cy + r * frac * Math.sin(a)];
  }

  for (const frac of [0.25, 0.5, 0.75, 1.0]) {
    const pts = Array.from({ length: n }, (_, i) => pt(frac, i));
    const poly = document.createElementNS(NS, 'polygon');
    poly.setAttribute('points', pts.map(p => p.join(',')).join(' '));
    poly.setAttribute('fill', 'none');
    poly.setAttribute('stroke', frac === 1.0 ? '#4f8ef788' : '#2e3246');
    poly.setAttribute('stroke-width', frac === 1.0 ? '1.5' : '0.7');
    if (frac === 1.0) poly.setAttribute('stroke-dasharray', '4 3');
    svg.appendChild(poly);
  }

  Array.from({ length: n }, (_, i) => {
    const [x, y] = pt(1, i);
    const line = document.createElementNS(NS, 'line');
    line.setAttribute('x1', cx); line.setAttribute('y1', cy);
    line.setAttribute('x2', x);  line.setAttribute('y2', y);
    line.setAttribute('stroke', '#2e3246');
    line.setAttribute('stroke-width', '0.7');
    svg.appendChild(line);
  });

  const clamped = axes.map(a => Math.min(Math.max(a.val, 0), 1.5));
  const allInRange = clamped.every(v => v <= 1.0);
  const dataPts = clamped.map((v, i) => pt(v, i));
  const poly = document.createElementNS(NS, 'polygon');
  poly.setAttribute('points', dataPts.map(p => p.join(',')).join(' '));
  poly.setAttribute('fill', allInRange ? 'rgba(52,211,153,.22)' : 'rgba(248,113,113,.18)');
  poly.setAttribute('stroke', allInRange ? '#34d399' : '#f87171');
  poly.setAttribute('stroke-width', '2');
  poly.setAttribute('stroke-linejoin', 'round');
  svg.appendChild(poly);

  dataPts.forEach(([x, y]) => {
    const dot = document.createElementNS(NS, 'circle');
    dot.setAttribute('cx', x); dot.setAttribute('cy', y); dot.setAttribute('r', '3');
    dot.setAttribute('fill', allInRange ? '#34d399' : '#f87171');
    svg.appendChild(dot);
  });

  axes.forEach((axis, i) => {
    const [x, y] = pt(1.32, i);
    const text = document.createElementNS(NS, 'text');
    text.setAttribute('x', x);
    text.setAttribute('y', y);
    text.setAttribute('text-anchor', 'middle');
    text.setAttribute('dominant-baseline', 'middle');
    text.setAttribute('font-size', '11');
    text.setAttribute('fill', '#8b92a9');
    text.setAttribute('font-family', 'var(--font-ui, sans-serif)');
    text.textContent = axis.label;
    svg.appendChild(text);
  });

  wrap.appendChild(svg);
}

function renderDescriptors(mol) {
  const tbody = document.getElementById('desc-tbody');
  tbody.textContent = '';
  for (const d of DESC_DEFS) {
    const val = d.fn(mol);
    const tr = document.createElement('tr');
    const tdName = document.createElement('td');
    tdName.className = 'td-name';
    tdName.textContent = t(d.key);
    const tdVal = document.createElement('td');
    tdVal.className = d.isSmiles ? 'td-smiles' : 'td-value';
    tdVal.textContent = val;
    tr.appendChild(tdName);
    tr.appendChild(tdVal);
    tbody.appendChild(tr);
  }
  const note = document.getElementById('rdkit-note');
  if (note) note.style.display = '';
}

function updateCodeSnippets(smiles) {
  const wrap = document.getElementById('code-wrap');
  if (!wrap) return;
  const s = smiles.replace(/\\/g, '\\\\').replace(/"/g, '\\"');
  wrap._pySnippet = `import chematic\n\nmol = chematic.from_smiles("${s}")\nprint(f"MW: {mol.mw:.2f}  LogP: {mol.logp:.3f}  TPSA: {mol.tpsa:.2f}")\nprint(mol.descriptors())`;
  wrap._jsSnippet = `import init, { parse_smiles } from '@kent-tokyo/chematic';\n\nawait init();\nconst mol = parse_smiles("${s}");\nconsole.log('MW:', mol.molecular_weight().toFixed(2));\nconsole.log('LogP:', mol.logp_crippen().toFixed(3));`;
  wrap.style.display = '';
}

function flashCopy(text) {
  navigator.clipboard.writeText(text).catch(() => {});
  const el = document.getElementById('copy-ok');
  if (!el) return;
  el.textContent = 'Copied!';
  el.style.opacity = '1';
  setTimeout(() => { el.style.opacity = '0'; }, 2000);
}

const NAMED_GROUP_COLORS = {
  hydroxyl:        '#22c55e', carbonyl:        '#f97316',
  carboxyl:        '#ef4444', aldehyde:        '#fb923c',
  ketone:          '#f59e0b', ester:           '#84cc16',
  amide:           '#34d399', primary_amine:   '#06b6d4',
  secondary_amine: '#0ea5e9', tertiary_amine:  '#3b82f6',
  nitro:           '#dc2626', nitrile:         '#8b5cf6',
  thiol:           '#eab308', sulfide:         '#ca8a04',
  sulfoxide:       '#d97706', sulfone:         '#b45309',
  halogen:         '#f43f5e', ether:           '#2dd4bf',
  epoxide:         '#a78bfa', aromatic:        '#60a5fa',
};
const _FALLBACK_PALETTE = ['#ef4444','#f97316','#eab308','#22c55e','#06b6d4','#3b82f6','#a855f7','#ec4899','#14b8a6'];

function renderNamedGroups(mol) {
  const out = document.getElementById('ifg-output');
  if (!out) return;
  out.textContent = '';
  namedGroupsData = [];
  activeGroupName = null;
  groupColorMap = {};

  if (!detectFunctionalGroups) {
    if (identifyFunctionalGroups) {
      try {
        const g = JSON.parse(identifyFunctionalGroups(mol));
        out.style.fontFamily = 'var(--font-mono)';
        out.style.fontSize = '.82rem';
        out.style.color = 'var(--muted)';
        out.textContent = g.length ? g.map(x => `[${x.atoms.join(',')}] ${x.type}`).join('  ·  ') : 'None';
      } catch(_) { out.textContent = '—'; }
    }
    return;
  }
  out.style.fontFamily = '';
  out.style.fontSize = '';
  out.style.color = '';

  try {
    const raw = JSON.parse(detectFunctionalGroups(mol));
    if (!raw.length) {
      const none = document.createElement('span');
      none.style.color = 'var(--muted)';
      none.style.fontSize = '.82rem';
      none.textContent = lang === 'ja' ? 'なし' : lang === 'zh' ? '无' : 'None';
      out.appendChild(none);
      return;
    }
    const byName = new Map();
    for (const g of raw) {
      if (!byName.has(g.name)) byName.set(g.name, new Set());
      g.atoms.forEach(a => byName.get(g.name).add(a));
    }
    let ci = 0;
    for (const [name, atomSet] of byName) {
      const color = NAMED_GROUP_COLORS[name] || _FALLBACK_PALETTE[ci++ % _FALLBACK_PALETTE.length];
      namedGroupsData.push({ name, atoms: [...atomSet], color });
    }
    const fgDiv = document.createElement('div');
    fgDiv.className = 'fg-chips';
    for (const g of namedGroupsData) {
      const chip = document.createElement('button');
      chip.className = 'fg-chip';
      chip.textContent = g.name.replace(/_/g, ' ');
      chip.style.borderColor = g.color;
      chip.style.color = g.color;
      chip.style.background = g.color + '22';
      chip.dataset.name = g.name;
      chip.addEventListener('click', () => toggleGroupHighlight(g.name));
      fgDiv.appendChild(chip);
    }
    out.appendChild(fgDiv);
  } catch(_) { out.textContent = '—'; }
}

function toggleGroupHighlight(name) {
  const g = namedGroupsData.find(x => x.name === name);
  if (!g) return;
  if (activeGroupName === name) {
    activeGroupName = null; groupColorMap = {}; highlightedAtoms = [];
  } else {
    activeGroupName = name; groupColorMap = {}; highlightedAtoms = [];
    for (const idx of g.atoms) groupColorMap[idx] = g.color;
  }
  document.querySelectorAll('.fg-chip').forEach(c => c.classList.toggle('active', c.dataset.name === activeGroupName));
  clearAtomInfo();
  if (lastMol) renderSvg(lastMol);
}

function _makeAiSpan(key, val) {
  const span = document.createElement('span');
  const k = document.createElement('span');
  k.className = 'ai-key';
  k.textContent = key;
  const v = document.createElement('span');
  v.className = 'ai-val';
  v.textContent = String(val);
  span.appendChild(k);
  span.appendChild(v);
  return span;
}

function showAtomInfo(idx) {
  const bar = document.getElementById('atom-info-bar');
  if (!bar || !getAtomInfo || !lastMol) return;
  try {
    const info = JSON.parse(getAtomInfo(lastMol, idx));
    if (!info) { clearAtomInfo(); return; }
    bar.classList.add('has-atom');
    bar.textContent = '';
    bar.appendChild(_makeAiSpan('#', idx));
    bar.appendChild(_makeAiSpan('El', info.element));
    bar.appendChild(_makeAiSpan('Hyb', info.hybridization));
    bar.appendChild(_makeAiSpan('H', info.totalHydrogens));
    bar.appendChild(_makeAiSpan('q', (info.charge >= 0 ? '+' : '') + info.charge));
    if (info.isAromatic) {
      const badge = document.createElement('span');
      badge.className = 'ai-badge';
      badge.textContent = 'aromatic';
      bar.appendChild(badge);
    }
  } catch(_) { clearAtomInfo(); }
}

function clearAtomInfo() {
  const bar = document.getElementById('atom-info-bar');
  if (!bar) return;
  bar.classList.remove('has-atom');
  bar.textContent = '';
  const hint = document.createElement('span');
  hint.id = 'atom-info-hint';
  hint.style.color = 'var(--muted)';
  hint.style.fontSize = '.8rem';
  hint.textContent = t('atom.hint');
  bar.appendChild(hint);
}

function downloadMol() {
  if (!lastMol || !molBlockFromSmiles) return;
  const smiles = document.getElementById('smiles-input').value.trim();
  if (!smiles) return;
  try {
    const molBlock = molBlockFromSmiles(smiles);
    const blob = new Blob([molBlock], { type: 'chemical/x-mdl-molfile' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url; a.download = 'molecule.mol'; a.click();
    URL.revokeObjectURL(url);
  } catch(_) {}
}

function renderRuleCards(mol) {
  const grid = document.getElementById('rules-grid');
  grid.textContent = '';
  for (const r of RULE_DEFS) {
    const passes = mol ? r.fn(mol) : null;

    const card = document.createElement('div');
    card.className = 'rule-card';

    const header = document.createElement('div');
    header.className = 'rule-header';

    const name = document.createElement('span');
    name.className = 'rule-name';
    name.textContent = t(r.nameKey);

    const badge = document.createElement('span');
    if (mol === null) {
      badge.className = 'badge badge-pending';
      badge.textContent = '—';
    } else {
      badge.className = 'badge ' + (passes ? 'badge-pass' : 'badge-fail');
      badge.textContent = passes ? t('rules.pass') : t('rules.fail');
    }

    header.appendChild(name);
    header.appendChild(badge);

    const criteria = document.createElement('div');
    criteria.className = 'rule-criteria';
    criteria.textContent = t(r.criteriaKey);

    card.appendChild(header);
    card.appendChild(criteria);
    grid.appendChild(card);
  }
}

function simColor(v) {
  if (v >= 0.7) return 'var(--pass)';
  if (v >= 0.3) return 'var(--warn)';
  return 'var(--fail)';
}

function renderSimTable(results) {
  lastSimResults = results;
  const tbody = document.getElementById('sim-tbody');
  tbody.textContent = '';
  for (let i = 0; i < SIM_DEFS.length; i++) {
    const d = SIM_DEFS[i];
    const v = results ? results[i] : null;

    const tr = document.createElement('tr');

    const tdMethod = document.createElement('td');
    tdMethod.className = 'td-name';
    tdMethod.style.fontWeight = '600';
    tdMethod.textContent = t(d.nameKey);

    const tdDesc = document.createElement('td');
    tdDesc.className = 'td-name';
    tdDesc.textContent = t(d.descKey);

    const tdVal = document.createElement('td');
    tdVal.className = 'sim-val';
    tdVal.style.color = v !== null ? simColor(v) : 'var(--muted)';
    tdVal.textContent = v !== null ? v.toFixed(3) : '—';

    const tdBar = document.createElement('td');
    tdBar.style.cssText = 'width:140px;vertical-align:middle;padding-right:1rem';
    const wrap = document.createElement('div');
    wrap.className = 'sim-bar-wrap';
    const bar = document.createElement('div');
    bar.className = 'sim-bar';
    bar.style.width = v !== null ? (v * 100).toFixed(1) + '%' : '0%';
    bar.style.background = v !== null ? simColor(v) : 'transparent';
    wrap.appendChild(bar);
    tdBar.appendChild(wrap);

    tr.appendChild(tdMethod);
    tr.appendChild(tdDesc);
    tr.appendChild(tdVal);
    tr.appendChild(tdBar);
    tbody.appendChild(tr);
  }
}

function renderSimSvgs(mA, mB) {
  function insertSvg(mol, boxEl) {
    boxEl.textContent = '';
    try {
      const opts = makeDepictOpts();
      const str = mol.depict_svg_opts(opts);
      opts.free();
      const parser = new DOMParser();
      const doc = parser.parseFromString(str, 'image/svg+xml');
      if (doc.querySelector('parsererror')) throw new Error('bad svg');
      boxEl.appendChild(doc.documentElement);
    } catch (_) {
      const ph = document.createElement('span');
      ph.style.cssText = 'color:#aaa;font-size:.85rem;padding:1rem';
      ph.textContent = '—';
      boxEl.appendChild(ph);
    }
  }

  insertSvg(mA, document.getElementById('sim-svg-a'));
  insertSvg(mB, document.getElementById('sim-svg-b'));
  document.getElementById('sim-svgs').style.display = '';
}

function renderGallery() {
  const container = document.getElementById('gallery-content');
  container.textContent = '';
  for (const cat of GALLERY) {
    const catDiv = document.createElement('div');
    catDiv.className = 'gallery-category';

    const title = document.createElement('div');
    title.className = 'gallery-cat-title';
    title.textContent = t(cat.catKey);
    catDiv.appendChild(title);

    const grid = document.createElement('div');
    grid.className = 'gallery-grid';

    for (const drug of cat.drugs) {
      const card = document.createElement('button');
      card.className = 'gallery-card';
      card.textContent = drug.name;
      card.dataset.smiles = drug.smiles;
      card.addEventListener('click', () => {
        document.querySelectorAll('.gallery-card.active').forEach(c => c.classList.remove('active'));
        card.classList.add('active');
        document.getElementById('smiles-input').value = drug.smiles;
        switchTab('tab-2d');
        calculate();
      });
      grid.appendChild(card);
    }

    catDiv.appendChild(grid);
    container.appendChild(catDiv);
  }
}

// ============================================================
// DepictOptions helper
// ============================================================
function makeDepictOpts() {
  const opts = new DepictOptions();
  opts.set_dark(true);
  opts.set_background('transparent');
  opts.set_atom_ids(true);
  if (showIndices) opts.set_show_atom_indices(true);
  return opts;
}

// ============================================================
// SVG pan + zoom
// ============================================================
function enableSvgPanZoom(area, svgEl) {
  let scale = 1, tx = 0, ty = 0;
  let dragging = false, startX = 0, startY = 0;

  function apply() {
    svgEl.style.transformOrigin = '50% 50%';
    svgEl.style.transform = `translate(${tx}px,${ty}px) scale(${scale})`;
  }

  area.addEventListener('wheel', e => {
    e.preventDefault();
    const f = e.deltaY < 0 ? 1.12 : 1 / 1.12;
    scale = Math.min(Math.max(scale * f, 0.15), 12);
    apply();
  }, { passive: false });

  area.addEventListener('pointerdown', e => {
    if (e.target.classList.contains('svg-zoom-reset')) return;
    if (e.target.classList.contains('btn-svg-ctrl')) return;
    dragging = true;
    startX = e.clientX - tx;
    startY = e.clientY - ty;
    area.setPointerCapture(e.pointerId);
    area.style.cursor = 'grabbing';
  });
  area.addEventListener('pointermove', e => {
    if (!dragging) return;
    tx = e.clientX - startX;
    ty = e.clientY - startY;
    apply();
  });
  area.addEventListener('pointerup', () => { dragging = false; area.style.cursor = 'grab'; });
  area.addEventListener('dblclick', reset);

  const btn = document.getElementById('svg-zoom-reset');
  if (btn) btn.onclick = reset;

  function reset() { scale = 1; tx = 0; ty = 0; apply(); }
  area.style.cursor = 'grab';
  apply();
}

// ============================================================
// SVG depiction
// ============================================================
function renderSvg(mol) {
  const area = document.getElementById('svg-area');
  try {
    const opts = makeDepictOpts();
    if (highlightedAtoms.length > 0) {
      opts.set_highlight_atoms(new Uint32Array(highlightedAtoms));
      opts.set_highlight_color(highlightColor);
    }
    for (const [idx, color] of Object.entries(groupColorMap)) {
      opts.set_atom_color(parseInt(idx), color);
    }
    const svgStr = mol.depict_svg_opts(opts);
    opts.free();
    const parser = new DOMParser();
    const doc = parser.parseFromString(svgStr, 'image/svg+xml');
    if (doc.querySelector('parsererror')) throw new Error('invalid SVG');

    const svgEl = doc.documentElement;
    svgEl.removeAttribute('width');
    svgEl.removeAttribute('height');
    svgEl.style.cssText = 'max-width:100%;max-height:360px;display:block;';

    area.querySelectorAll('svg').forEach(el => el.remove());
    const ph = document.getElementById('svg-placeholder');
    if (ph) ph.textContent = '';
    area.appendChild(svgEl);
    enableSvgPanZoom(area, svgEl);
    enableAtomClick(svgEl);
  } catch (_) {
    area.querySelectorAll('svg').forEach(el => el.remove());
    const ph = document.getElementById('svg-placeholder');
    if (ph) { ph.style.color = 'var(--muted)'; ph.textContent = t('desc.svg_ph'); }
  }
}

// ============================================================
// Atom click — highlight on click, deselect on re-click
// ============================================================
function enableAtomClick(svgEl) {
  svgEl.addEventListener('click', e => {
    const atomEl = e.target.closest('[data-atom-idx]');
    if (atomEl) {
      const idx = parseInt(atomEl.getAttribute('data-atom-idx'), 10);
      if (highlightedAtoms.length === 1 && highlightedAtoms[0] === idx) {
        highlightedAtoms = [];
        clearAtomInfo();
      } else {
        highlightedAtoms = [idx];
        highlightColor = '#4f8ef7';
        showAtomInfo(idx);
      }
      clearSmartsUI();
      renderSvg(lastMol);
    } else {
      if (highlightedAtoms.length > 0 || Object.keys(groupColorMap).length > 0) {
        highlightedAtoms = [];
        groupColorMap = {};
        activeGroupName = null;
        document.querySelectorAll('.fg-chip').forEach(c => c.classList.remove('active'));
        clearAtomInfo();
        renderSvg(lastMol);
      }
    }
  });
}

// ============================================================
// SMARTS substructure highlight
// ============================================================
function applySmarts() {
  if (!lastMol || !smartsMatchAtoms) return;
  const pat = document.getElementById('smarts-input').value.trim();
  if (!pat) return;
  clearError('error-desc');
  try {
    const json = smartsMatchAtoms(pat, lastMol);
    const matches = JSON.parse(json);
    const idxs = [...new Set(matches.flat())].map(Number);
    const badge = document.getElementById('smarts-match-badge');
    const nomatch = document.getElementById('smarts-nomatch-msg');
    if (idxs.length > 0) {
      highlightedAtoms = idxs;
      highlightColor = '#34d399';
      badge.textContent = idxs.length + ' ' + t('smarts.matches');
      badge.classList.add('visible');
      nomatch.classList.remove('visible');
    } else {
      highlightedAtoms = [];
      badge.classList.remove('visible');
      nomatch.classList.add('visible');
    }
    renderSvg(lastMol);
  } catch (e) {
    showError('error-desc', 'SMARTS: ' + String(e));
  }
}

function clearSmarts() {
  clearSmartsUI();
  highlightedAtoms = [];
  if (lastMol) renderSvg(lastMol);
}

function clearSmartsUI() {
  const badge = document.getElementById('smarts-match-badge');
  const nomatch = document.getElementById('smarts-nomatch-msg');
  if (badge) badge.classList.remove('visible');
  if (nomatch) nomatch.classList.remove('visible');
}

// ============================================================
// Reaction Scheme
// ============================================================
function runReaction() {
  clearError('error-rxn');
  const smirks = document.getElementById('rxn-smirks').value.trim();
  const reactants = document.getElementById('rxn-reactants').value.trim();
  if (!smirks || !reactants) { showError('error-rxn', t('error.empty')); return; }
  if (!wasmReady) { showError('error-rxn', t('error.wasm')); return; }
  try {
    const productsJson = runReactants(smirks, reactants);
    const products = JSON.parse(productsJson);
    const reactantList = reactants.split('|').map(s => s.trim()).filter(Boolean);
    const productList = products.length > 0 ? products[0] : [];
    const allSmiles = [...reactantList, ...productList];
    if (allSmiles.length === 0) {
      showError('error-rxn', 'No products generated.');
      return;
    }
    const cols = Math.max(reactantList.length, productList.length, 1);
    const gridSvg = depictSvgGrid(allSmiles.join('\n'), cols);
    const parser = new DOMParser();
    const doc = parser.parseFromString(gridSvg, 'image/svg+xml');
    if (doc.querySelector('parsererror')) throw new Error('Invalid SVG from grid renderer');
    const svgEl = doc.documentElement;
    svgEl.style.cssText = 'max-width:100%;height:auto;display:block;';
    const wrap = document.getElementById('rxn-svg-wrap');
    const ph = wrap.querySelector('.rxn-placeholder');
    if (ph) ph.style.display = 'none';
    wrap.querySelectorAll('svg').forEach(el => el.remove());
    wrap.appendChild(svgEl);
  } catch (e) {
    showError('error-rxn', String(e));
  }
}

// ============================================================
// Reaction Equation Viewer
// ============================================================
function runReactionEquation() {
  clearError('error-rxn-eq');
  const smi = document.getElementById('rxn-eq-input').value.trim();
  if (!smi) { showError('error-rxn-eq', t('error.empty')); return; }
  if (!wasmReady) { showError('error-rxn-eq', t('error.wasm')); return; }
  try {
    const svg = depictReactionSvg(smi);
    const parser = new DOMParser();
    const doc = parser.parseFromString(svg, 'image/svg+xml');
    if (doc.querySelector('parsererror')) throw new Error('Invalid SVG from renderer');
    const svgEl = doc.documentElement;
    svgEl.style.cssText = 'max-width:100%;height:auto;display:block;';
    const wrap = document.getElementById('rxn-eq-svg-wrap');
    const ph = wrap.querySelector('.rxn-placeholder');
    if (ph) ph.style.display = 'none';
    wrap.querySelectorAll('svg').forEach(el => el.remove());
    wrap.appendChild(svgEl);
  } catch (e) {
    showError('error-rxn-eq', String(e));
  }
}

// ============================================================
// Error helpers
// ============================================================
function showError(id, msg) {
  const el = document.getElementById(id);
  if (!el) return;
  el.textContent = msg;
  el.classList.remove('hidden');
}
function clearError(id) {
  const el = document.getElementById(id);
  if (el) el.classList.add('hidden');
}

// ============================================================
// Calculate descriptors
// ============================================================
function calculate() {
  clearError('error-desc');
  const smiles = document.getElementById('smiles-input').value.trim();
  if (!smiles) { showError('error-desc', t('error.empty')); return; }
  if (!wasmReady) { showError('error-desc', t('error.wasm')); return; }

  let mol;
  try { mol = parseSmilesWasm(smiles); }
  catch (e) { showError('error-desc', t('error.parse') + String(e)); return; }

  if (lastMol) { try { lastMol.free(); } catch (_) {} }
  lastMol = mol;
  highlightedAtoms = [];
  groupColorMap = {};
  activeGroupName = null;
  clearSmartsUI();
  clearAtomInfo();

  // Critical path: SVG + radar + primary descriptors render immediately.
  renderSvg(mol);
  renderRadar(mol);
  renderDescriptors(mol);
  updateCodeSnippets(smiles);

  // Secondary: rule cards + functional groups are visually below the fold;
  // defer to keep the main-thread responsive after each keystroke.
  const rIC = typeof requestIdleCallback === 'function' ? requestIdleCallback : (fn) => setTimeout(fn, 0);
  rIC(() => { if (lastMol === mol) renderRuleCards(mol); });
  rIC(() => { if (lastMol === mol) renderNamedGroups(mol); });
  // Clear stale stereo isomers when molecule changes.
  document.getElementById('stereo-grid').textContent = '';

  // Mark 3D dirty — update when user switches to 3D tab
  if (activeTab === 'tab-3d' && viewer3dInited && viewer3d) {
    render3d(mol);
  } else {
    viewer3dMolDirty = true;
  }
}

// ============================================================
// Compare similarity
// ============================================================
function compare() {
  clearError('error-sim');
  const sA = document.getElementById('sim-a').value.trim();
  const sB = document.getElementById('sim-b').value.trim();
  if (!sA || !sB) { showError('error-sim', t('error.empty')); return; }
  if (!wasmReady) { showError('error-sim', t('error.wasm')); return; }

  let mA, mB;
  try { mA = parseSmilesWasm(sA); } catch(e) { showError('error-sim', 'A: ' + t('error.parse') + String(e)); return; }
  try { mB = parseSmilesWasm(sB); } catch(e) { mA.free(); showError('error-sim', 'B: ' + t('error.parse') + String(e)); return; }

  renderSimSvgs(mA, mB);
  const results = SIM_DEFS.map(d => d.fn(mA, mB));
  mA.free(); mB.free();
  renderSimTable(results);
}

// ============================================================
// SDF / MOL load
// ============================================================
function loadSdf() {
  clearError('error-sdf');
  const text = document.getElementById('sdf-input').value.trim();
  if (!text) { showError('error-sdf', t('error.empty')); return; }
  if (!wasmReady) { showError('error-sdf', t('error.wasm')); return; }

  const out = document.getElementById('sdf-grid-output');
  out.textContent = '';

  let smilesList = [];
  try {
    const json = sdfToSmilesJson(text);
    const parsed = JSON.parse(json);
    smilesList = parsed.filter(s => s !== null);
  } catch (e) {
    showError('error-sdf', String(e));
    return;
  }

  if (smilesList.length === 0) {
    showError('error-sdf', 'No valid molecules found in input.');
    return;
  }

  // Load first molecule into the descriptor panel
  try {
    const firstMol = parseSmilesWasm(smilesList[0]);
    lastMol = firstMol;
    groupColorMap = {};
    activeGroupName = null;
    clearAtomInfo();
    renderSvg(firstMol);
    renderDescriptors(firstMol);
    renderRuleCards(firstMol);
    renderNamedGroups(firstMol);
    viewer3dMolDirty = true;
  } catch (_) {}

  // Render molecule grid using DOMParser (same pattern as runReaction/runReactionEquation)
  try {
    const gridSvg = depictSvgGrid(smilesList.join('\n'), 4);
    const parser = new DOMParser();
    const doc = parser.parseFromString(gridSvg, 'image/svg+xml');
    if (!doc.querySelector('parsererror')) {
      const svgEl = doc.documentElement;
      svgEl.style.cssText = 'max-width:100%;height:auto;display:block;';
      out.appendChild(svgEl);
    }
  } catch (_) {}
}

// ============================================================
// Stereo isomers
// ============================================================
function showStereoIsomers() {
  const grid = document.getElementById('stereo-grid');
  grid.textContent = '';
  clearError('error-stereo');
  if (!lastMol) return;
  if (!enumerateStereoIsomers) { showError('error-stereo', 'enumerate_stereo_isomers_json not available'); return; }
  try {
    const result = enumerateStereoIsomers(lastMol);
    const isomers = JSON.parse(result);
    if (!isomers.length) {
      const msg = document.createElement('span');
      msg.style.cssText = 'color:var(--muted);font-size:.82rem';
      msg.textContent = 'No unspecified stereocenters found.';
      grid.appendChild(msg);
      return;
    }
    for (const smi of isomers) {
      const cell = document.createElement('div');
      cell.className = 'stereo-cell';
      try {
        const svgHtml = depictSvgGrid(smi, 1);
        const parser = new DOMParser();
        const doc = parser.parseFromString(svgHtml, 'image/svg+xml');
        if (!doc.querySelector('parsererror')) cell.appendChild(doc.documentElement);
      } catch (_) {}
      const smiLabel = document.createElement('div');
      smiLabel.className = 'stereo-smi';
      smiLabel.textContent = smi;
      cell.appendChild(smiLabel);
      grid.appendChild(cell);
    }
  } catch (e) {
    showError('error-stereo', String(e));
  }
}

// ============================================================
// MMP analysis
// ============================================================
const MMP_PRESETS = {
  alkylbenzenes: 'CCc1ccccc1\nCCCc1ccccc1\nCCCCc1ccccc1\nCc1ccccc1C\nCC(C)c1ccccc1',
  nsaids: 'CC(=O)Oc1ccccc1C(=O)O\nCC(C)Cc1ccc(cc1)C(C)C(=O)O\nCC(C(=O)O)c1ccc2cccc(OC)c2c1',
};

function makeCodeCell(text) {
  const td = document.createElement('td');
  const code = document.createElement('code');
  code.className = 'mmp-smi-code';
  code.textContent = text;
  td.appendChild(code);
  return td;
}

function runMmp() {
  const tbody = document.getElementById('mmp-tbody');
  const table = document.getElementById('mmp-table');
  tbody.textContent = '';
  table.style.display = 'none';
  clearError('error-mmp');
  if (!mmpPairsJson) { showError('error-mmp', 'mmp_pairs_json not available'); return; }

  const lines = document.getElementById('mmp-input').value
    .split('\n').map(s => s.trim()).filter(Boolean);
  if (lines.length < 2) { showError('error-mmp', 'Enter at least 2 SMILES (one per line).'); return; }

  try {
    const pairs = JSON.parse(mmpPairsJson(JSON.stringify(lines)));
    if (!pairs.length) { showError('error-mmp', 'No matched molecular pairs found.'); return; }
    for (const p of pairs) {
      const tr = document.createElement('tr');
      for (const s of [p.mol_a, p.mol_b, p.core, p.fragment_a, p.fragment_b]) {
        tr.appendChild(makeCodeCell(s));
      }
      tbody.appendChild(tr);
    }
    table.style.display = '';
  } catch (e) {
    showError('error-mmp', String(e));
  }
}

// ============================================================
// R-group decomposition
// ============================================================
const RGROUP_PRESETS = {
  benzene: { smiles: 'CCc1ccccc1\nCCCc1ccccc1\nFc1ccccc1\nClc1ccccc1\nCOc1ccccc1', core: 'c1ccc(*)cc1' },
  pyridine: { smiles: 'CCc1ccncc1\nFc1ccncc1\nClc1ccncc1\nCOc1ccncc1', core: 'c1cc(*)ccn1' },
};

function runRgroup() {
  const grid = document.getElementById('rgroup-grid');
  grid.textContent = '';
  clearError('error-rgroup');
  if (!rgroupDecomposeJson) { showError('error-rgroup', 'rgroup_decompose_json not available'); return; }

  const lines = document.getElementById('rgroup-input').value
    .split('\n').map(s => s.trim()).filter(Boolean);
  const core = document.getElementById('rgroup-core').value.trim();
  if (!lines.length) { showError('error-rgroup', 'Enter at least 1 SMILES.'); return; }
  if (!core) { showError('error-rgroup', 'Enter a core SMARTS pattern.'); return; }

  try {
    const results = JSON.parse(rgroupDecomposeJson(JSON.stringify(lines), core));
    for (let i = 0; i < results.length; i++) {
      const r = results[i];
      const smi = lines[i];
      const cell = document.createElement('div');
      cell.className = 'rgroup-cell';

      try {
        const svgHtml = depictSvgGrid(smi, 1);
        const parser = new DOMParser();
        const doc = parser.parseFromString(svgHtml, 'image/svg+xml');
        if (!doc.querySelector('parsererror')) cell.appendChild(doc.documentElement);
      } catch (_) {}

      const info = document.createElement('div');
      info.className = 'rgroup-info';

      const status = document.createElement('span');
      status.className = r.matched ? 'matched' : 'unmatched';
      status.textContent = r.matched ? '✓ matched' : '✗ no match';
      info.appendChild(status);

      if (r.matched) {
        const rkeys = Object.keys(r).filter(k => k !== 'matched');
        for (const k of rkeys) {
          const br = document.createElement('br');
          info.appendChild(br);
          const line = document.createElement('span');
          line.textContent = k + ': ';
          const code = document.createElement('code');
          code.textContent = r[k];
          info.appendChild(line);
          info.appendChild(code);
        }
      }
      cell.appendChild(info);
      grid.appendChild(cell);
    }
  } catch (e) {
    showError('error-rgroup', String(e));
  }
}

// ============================================================
// Formats & Materials tab
// ============================================================
let formatsMode = 'cube';

const FORMATS_EXAMPLES = {
  cube: `Water density
Generated for chematic examples
1    0.000000    0.000000    0.000000
2    1.000000    0.000000    0.000000
2    0.000000    1.000000    0.000000
2    0.000000    0.000000    1.000000
8    8.000000    0.500000    0.500000    0.500000
0.0 1.0 2.0 3.0
4.0 5.0 6.0 7.0
`,
  opendx: `object 1 class gridpositions counts 2 2 2
origin -1.0 -1.0 -1.0
delta 0.5 0.0 0.0
delta 0.0 0.5 0.0
delta 0.0 0.0 0.5
object 2 class gridconnections counts 2 2 2
object 3 class array type double rank 0 items 8 data follows
0.0 1.0 2.0
3.0 4.0 5.0
6.0 7.0
attribute "dep" string "positions"
object "regular positions regular connections" class field
component "positions" value 1
component "connections" value 2
component "data" value 3
`,
  mmcif: `data_EXAMPLE
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
ATOM   1  O  O1  HOH A 1.000 2.000 3.000
ATOM   2  H  H1  HOH A 1.500 2.500 3.500
`,
  pqr: `ATOM      1  N   ALA     1     -0.966   1.523   1.412 -0.400  1.500
ATOM      2  CA  ALA     1      0.257   0.679   1.911  0.100  1.700
HETATM    3  ZN  ZN      2      5.000   5.000   5.000  2.000  1.090
`,
  qcschema: JSON.stringify({
    schema_name: 'qcschema_molecule',
    schema_version: 1,
    symbols: ['O', 'H', 'H'],
    geometry: [
      0.0, 0.0, 0.0,
      0.0, 1.431544636036637, 1.109419726497757,
      0.0, -1.431544636036637, 1.109419726497757,
    ],
    molecular_charge: 0.0,
    molecular_multiplicity: 1,
    connectivity: [[0, 1, 1.0], [0, 2, 1.0]],
  }, null, 2),
  orca_input: `! B3LYP def2-SVP
! Opt Freq

* xyz 0 1
O 0 0 0
*
`,
  orca_output: `
 Total Charge           Charge          ....    0
 Multiplicity           Mult            ....    1

-------------------------   --------------------
FINAL SINGLE POINT ENERGY       -76.025678123456
-------------------------   --------------------

****ORCA TERMINATED NORMALLY****
`,
  lammps_data: `comment

1 atoms
1 atom types

0.0 10.0 xlo xhi
0.0 10.0 ylo yhi
0.0 10.0 zlo zhi

Atoms # atomic

1 1 5.0 5.0 5.0
`,
  lammps_dump:
    'ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n2\nITEM: BOX BOUNDS pp pp pp\n' +
    '0.0 10.0\n0.0 10.0\n0.0 10.0\nITEM: ATOMS id x y z\n1 1.0 2.0 3.0\n2 4.0 5.0 6.0\n',
};

const FORMATS_MODES = {
  cube: {
    parse(text) {
      const mol = molFromCube(text);
      const atomCount = mol.atom_count();
      mol.free();
      const grid = JSON.parse(cubeGridJson(text));
      const values = cubeValuesF64(text);
      return {
        summary: [
          [t('formats.f.atoms'), atomCount],
          [t('formats.f.shape'), grid.shape.join(' × ')],
          [t('formats.f.units'), grid.units],
          [t('formats.f.values'), values.length],
        ],
        raw: grid,
        histogram: values,
      };
    },
  },
  opendx: {
    parse(text) {
      const grid = JSON.parse(opendxGridJson(text));
      const values = opendxValuesF64(text);
      return {
        summary: [
          [t('formats.f.shape'), grid.shape.join(' × ')],
          [t('formats.f.units'), grid.units],
          [t('formats.f.values'), values.length],
        ],
        raw: grid,
        histogram: values,
      };
    },
  },
  mmcif: {
    parse(text) {
      const mol = molFromMmcif(text);
      const atomCount = mol.atom_count();
      mol.free();
      const result = JSON.parse(mmcifToJson(text));
      const summary = [[t('formats.f.atoms'), atomCount]];
      if (result.cell) {
        const c = result.cell;
        summary.push([t('formats.f.cell'),
          `a=${c.a} b=${c.b} c=${c.c} α=${c.alpha} β=${c.beta} γ=${c.gamma}`]);
      }
      if (result.space_group) summary.push([t('formats.f.space_group'), result.space_group]);
      if (result.unhandled_columns && result.unhandled_columns.length) {
        summary.push([t('formats.f.warnings'), result.unhandled_columns.join(', ')]);
      }
      return { summary, raw: result, histogram: null };
    },
  },
  pqr: {
    parse(text) {
      const mol = molFromPqr(text);
      const atomCount = mol.atom_count();
      mol.free();
      const result = JSON.parse(pqrToJson(text));
      return { summary: [[t('formats.f.atoms'), atomCount]], raw: result, histogram: null };
    },
  },
  qcschema: {
    parse(text) {
      const mol = molFromQcschemaMolecule(text);
      const atomCount = mol.atom_count();
      mol.free();
      const result = JSON.parse(qcschemaMoleculeCoordsJson(text));
      return {
        summary: [
          [t('formats.f.atoms'), atomCount],
          [t('formats.f.charge'), result.molecular_charge],
          [t('formats.f.multiplicity'), result.molecular_multiplicity],
        ],
        raw: result,
        histogram: null,
      };
    },
  },
  orca_input: {
    parse(text) {
      let atomCount = null;
      try {
        const mol = molFromOrcaInput(text);
        atomCount = mol.atom_count();
        mol.free();
      } catch (_) {
        // xyzfile/gzmtfile/int (Z-matrix) coordinate block, or none at all --
        // orca_input_to_json below still describes it, just without an atom count.
      }
      const result = JSON.parse(orcaInputToJson(text));
      const summary = [
        [t('formats.f.keywords'), result.keywords.join(' ')],
        [t('formats.f.blocks'), result.blocks.length],
      ];
      if (atomCount !== null) summary.unshift([t('formats.f.atoms'), atomCount]);
      if (result.coords) summary.push([t('formats.f.coord_kind'), result.coords.type]);
      return { summary, raw: result, histogram: null };
    },
  },
  orca_output: {
    parse(text) {
      const result = JSON.parse(orcaOutputToJson(text));
      const summary = [
        [t('formats.f.energy'), result.final_energy_hartree],
        [t('formats.f.termination'), result.termination.kind],
        [t('formats.f.convergence'), result.optimization_convergence],
        [t('formats.f.frames'), result.trajectory.length],
      ];
      return { summary, raw: result, histogram: null };
    },
  },
  lammps_data: {
    needsAtomStyle: true,
    parse(text, atomStyle) {
      const result = JSON.parse(lammpsDataToJson(text, atomStyle));
      const summary = [
        [t('formats.f.atom_style'), result.atom_style],
        [t('formats.f.atoms'), result.atoms.length],
        [t('formats.f.bonds'), result.bonds.length],
      ];
      if (result.unparsed_sections.length) {
        summary.push([t('formats.f.warnings'), result.unparsed_sections.map(s => s[0]).join(', ')]);
      }
      return { summary, raw: result, histogram: null };
    },
  },
  lammps_dump: {
    parse(text) {
      const frames = JSON.parse(lammpsTrajectoryToJson(text));
      const summary = [[t('formats.f.frames'), frames.length]];
      if (frames[0]) {
        summary.push([t('formats.f.atoms'), frames[0].num_atoms]);
        summary.push([t('formats.f.timestep'), frames[0].timestep]);
      }
      return { summary, raw: frames, histogram: null };
    },
  },
};

function formatsExampleMatches(value, mode) {
  return value.trim() === (FORMATS_EXAMPLES[mode] || '').trim();
}

function setFormatsMode(mode) {
  const input = document.getElementById('formats-input');
  const previousMode = formatsMode;
  // Keep an uploaded or hand-written document intact. The only safe automatic
  // replacement is the built-in example from the mode the user just left.
  // Otherwise switching the parser would make the visible example invalid.
  if (!input.value.trim() || formatsExampleMatches(input.value, previousMode)) {
    input.value = FORMATS_EXAMPLES[mode] || '';
  }
  formatsMode = mode;
  document.querySelectorAll('.format-mode-btn').forEach(btn => {
    btn.classList.toggle('active', btn.dataset.mode === mode);
  });
  document.getElementById('formats-atom-style-row').style.display =
    FORMATS_MODES[mode].needsAtomStyle ? '' : 'none';
  document.getElementById('formats-output').style.display = 'none';
  clearError('error-formats');
}

function loadFormatsExample() {
  document.getElementById('formats-input').value = FORMATS_EXAMPLES[formatsMode] || '';
  clearError('error-formats');
}

function readFormatsFile(file) {
  const reader = new FileReader();
  reader.onload = () => {
    document.getElementById('formats-input').value = String(reader.result);
    clearError('error-formats');
  };
  reader.onerror = () => showError('error-formats', String(reader.error));
  reader.readAsText(file);
}

function renderFormatsSummary(rows) {
  const tbody = document.getElementById('formats-summary-tbody');
  tbody.textContent = '';
  for (const [name, value] of rows) {
    const tr = document.createElement('tr');
    const tdName = document.createElement('td');
    tdName.className = 'td-name';
    tdName.textContent = name;
    const tdVal = document.createElement('td');
    tdVal.className = 'td-value';
    tdVal.textContent = String(value);
    tr.appendChild(tdName);
    tr.appendChild(tdVal);
    tbody.appendChild(tr);
  }
}

// Dependency-free bar histogram of a scalar field (Cube/OpenDX grid values).
// Min/max/bin-count computed via plain loops, not Math.min(...values) /
// Math.max(...values) -- spreading a real (potentially large) grid's values
// into call arguments risks a stack-size blowup that a loop never does.
function renderFormatsHistogram(values) {
  const canvas = document.getElementById('formats-histogram');
  if (!values || values.length === 0) { canvas.style.display = 'none'; return; }
  canvas.style.display = '';

  let min = Infinity, max = -Infinity;
  for (const v of values) { if (v < min) min = v; if (v > max) max = v; }
  const binCount = 20;
  const bins = new Array(binCount).fill(0);
  const range = (max - min) || 1;
  for (const v of values) {
    let idx = Math.floor((v - min) / range * binCount);
    if (idx >= binCount) idx = binCount - 1;
    if (idx < 0) idx = 0;
    bins[idx]++;
  }
  let maxCount = 0;
  for (const c of bins) if (c > maxCount) maxCount = c;
  maxCount = maxCount || 1;

  const ctx = canvas.getContext('2d');
  const w = canvas.width, h = canvas.height;
  ctx.clearRect(0, 0, w, h);
  const accentVar = getComputedStyle(document.documentElement).getPropertyValue('--accent').trim();
  ctx.fillStyle = accentVar || '#4f8ef7';
  const barW = w / binCount;
  for (let i = 0; i < binCount; i++) {
    const barH = (bins[i] / maxCount) * (h - 4);
    ctx.fillRect(i * barW + 1, h - barH, barW - 2, barH);
  }
}

function runFormatsParse() {
  clearError('error-formats');
  const text = document.getElementById('formats-input').value.trim();
  if (!text) { showError('error-formats', t('error.empty')); return; }
  if (!wasmReady) { showError('error-formats', t('error.wasm')); return; }

  const atomStyle = document.getElementById('formats-atom-style').value;
  let result;
  try {
    result = FORMATS_MODES[formatsMode].parse(text, atomStyle);
  } catch (e) {
    document.getElementById('formats-output').style.display = 'none';
    showError('error-formats', String(e));
    return;
  }

  renderFormatsSummary(result.summary);
  renderFormatsHistogram(result.histogram);
  document.getElementById('formats-raw-json').textContent = JSON.stringify(result.raw, null, 2);
  document.getElementById('formats-output').style.display = '';
}

function initFormatsTab() {
  document.querySelectorAll('.format-mode-btn').forEach(btn => {
    btn.addEventListener('click', () => setFormatsMode(btn.dataset.mode));
  });
  document.getElementById('btn-formats-example').addEventListener('click', loadFormatsExample);
  document.getElementById('btn-formats-parse').addEventListener('click', runFormatsParse);

  const fileInput = document.getElementById('formats-file-input');
  document.getElementById('btn-formats-browse').addEventListener('click', () => fileInput.click());
  fileInput.addEventListener('change', () => {
    if (fileInput.files[0]) readFormatsFile(fileInput.files[0]);
    fileInput.value = '';
  });

  const dropzone = document.getElementById('formats-dropzone');
  dropzone.addEventListener('dragover', (e) => { e.preventDefault(); dropzone.classList.add('drag-over'); });
  dropzone.addEventListener('dragleave', () => dropzone.classList.remove('drag-over'));
  dropzone.addEventListener('drop', (e) => {
    e.preventDefault();
    dropzone.classList.remove('drag-over');
    const file = e.dataTransfer.files[0];
    if (file) readFormatsFile(file);
  });

  setFormatsMode('cube');
}

// ============================================================
// Boot
// ============================================================
(async () => {
  const saved = localStorage.getItem('chematic-lang');
  const browser = navigator.language.toLowerCase();
  const detected = saved ?? (browser.startsWith('zh') ? 'zh' : browser.startsWith('ja') ? 'ja' : 'en');
  applyLang(detected);
  clearAtomInfo();

  renderRuleCards(null);
  renderSimTable(null);
  renderGallery();

  await initWasm();
  document.getElementById('loading-overlay').classList.add('hidden');

  // Auto-load the shared molecule when present; otherwise show a useful first
  // example without making the user understand SMILES before seeing a result.
  document.getElementById('smiles-input').value = sharedMoleculeFromUrl() || 'CC(=O)Oc1ccccc1C(=O)O';
  calculate();

  // Batch Report: generate downloadable HTML from multiple SMILES
  document.getElementById('btn-batch-report').addEventListener('click', () => {
    clearError('error-batch-report');
    const raw = document.getElementById('batch-report-smiles').value.trim();
    if (!raw) { showError('error-batch-report', 'Enter at least one SMILES'); return; }
    if (!wasmReady || typeof window.chemWasm.batch_report_html !== 'function') {
      showError('error-batch-report', 'batch_report_html not available'); return;
    }
    try {
      const html = window.chemWasm.batch_report_html(raw);
      const blob = new Blob([html], { type: 'text/html' });
      const url = URL.createObjectURL(blob);
      const a = document.getElementById('batch-report-download');
      a.href = url;
      a.style.display = '';
    } catch (e) {
      showError('error-batch-report', 'Report generation failed: ' + String(e));
    }
  });

  // Report Tab: Populate descriptors table from workflow JSON
  function populateReportDescriptors(report) {
    const tbody = document.getElementById('report-desc-tbody');
    tbody.textContent = '';
    if (!report.descriptors) return;

    const desc = report.descriptors;
    const rows = [
      ['Canonical SMILES', report.canonical_smiles],
      ['Molecular Formula', report.formula],
      ['Molecular Weight (Da)', desc.molecular_weight.toFixed(2)],
      ['Exact Mass (Da)', desc.exact_mass.toFixed(4)],
      ['TPSA (Ų)', desc.tpsa.toFixed(2)],
      ['LogP (Crippen)', desc.logp.toFixed(3)],
      ['Fsp3', desc.fsp3.toFixed(3)],
      ['QED', desc.qed.toFixed(3)],
      ['Molar Refractivity', desc.molar_refractivity.toFixed(2)],
      ['Heavy Atoms', String(desc.heavy_atom_count)],
      ['HBD', String(desc.hbd)],
      ['HBA', String(desc.hba)],
      ['Rotatable Bonds', String(desc.rotatable_bonds)],
      ['Ring Count (SSSR)', String(desc.ring_count)],
      ['Heteroatoms', String(desc.num_heteroatoms)],
      ['Stereocenters', String(desc.num_stereocenters)],
      ['SA Score', desc.sa_score.toFixed(2)],
      ['Bertz CT', desc.bertz_ct.toFixed(2)],
      ['Wiener Index', desc.wiener_index.toFixed(2)],
    ];

    for (const [name, val] of rows) {
      const tr = document.createElement('tr');
      const tdName = document.createElement('td');
      tdName.className = 'td-name';
      tdName.textContent = name;
      const tdVal = document.createElement('td');
      tdVal.className = 'td-value';
      tdVal.textContent = String(val);
      tr.appendChild(tdName);
      tr.appendChild(tdVal);
      tbody.appendChild(tr);
    }
  }

  // Report Tab: Populate drug-likeness filters grid from workflow JSON
  function populateReportFilters(report) {
    const grid = document.getElementById('report-filters-grid');
    grid.textContent = '';
    if (!report.filters) return;

    const filters = report.filters;
    const filterList = [
      { name: 'Lipinski (Ro5)',      passes: filters.lipinski_passes },
      { name: 'Veber',               passes: filters.veber_passes },
      { name: 'Egan',                passes: filters.egan_passes },
      { name: 'Ghose',               passes: filters.ghose_passes },
      { name: 'REOS',                passes: filters.reos_passes },
      { name: 'PAINS',               passes: filters.pains_passes },
    ];

    for (const rule of filterList) {
      const card = document.createElement('div');
      card.className = 'rule-card';

      const header = document.createElement('div');
      header.className = 'rule-header';

      const name = document.createElement('span');
      name.className = 'rule-name';
      name.textContent = rule.name;

      const badge = document.createElement('span');
      badge.className = 'badge ' + (rule.passes ? 'badge-pass' : 'badge-fail');
      badge.textContent = rule.passes ? t('rules.pass') : t('rules.fail');

      header.appendChild(name);
      header.appendChild(badge);
      card.appendChild(header);
      grid.appendChild(card);
    }

    // Show PAINS alerts if any
    if (filters.pains_alerts && filters.pains_alerts.length > 0) {
      const alertCard = document.createElement('div');
      alertCard.className = 'rule-card';
      const header = document.createElement('div');
      header.className = 'rule-header';
      const name = document.createElement('span');
      name.className = 'rule-name';
      name.textContent = 'PAINS Alerts';
      header.appendChild(name);
      alertCard.appendChild(header);
      const alerts = document.createElement('div');
      alerts.className = 'rule-criteria';
      alerts.textContent = filters.pains_alerts.join(', ');
      alertCard.appendChild(alerts);
      grid.appendChild(alertCard);
    }
  }

  // Tab switching
  document.querySelectorAll('.tab-btn').forEach(btn => {
    btn.addEventListener('click', () => switchTab(btn.dataset.tab));
  });
  document.querySelectorAll('.tab-cat-btn').forEach(btn => {
    btn.addEventListener('click', () => setCategory(btn.dataset.cat));
  });

  // Language switcher
  document.querySelectorAll('[data-lang]').forEach(btn => {
    btn.addEventListener('click', () => applyLang(btn.dataset.lang));
  });

  // Calculate (button + Enter + debounced input)
  document.getElementById('btn-calc').addEventListener('click', calculate);
  document.getElementById('smiles-input').addEventListener('keydown', e => {
    if (e.key === 'Enter') calculate();
  });
  let calcTimer;
  document.getElementById('smiles-input').addEventListener('input', () => {
    clearTimeout(calcTimer);
    calcTimer = setTimeout(() => { if (wasmReady) calculate(); }, 350);
  });

  // Descriptor presets
  document.querySelectorAll('.preset-btn').forEach(btn => {
    btn.addEventListener('click', () => {
      document.getElementById('smiles-input').value = btn.dataset.smiles;
      calculate();
    });
  });

  // Similarity presets
  document.querySelectorAll('.preset-sim').forEach(btn => {
    btn.addEventListener('click', () => {
      document.getElementById(btn.dataset.target === 'a' ? 'sim-a' : 'sim-b').value = btn.dataset.smiles;
    });
  });

  // Compare (button + Enter)
  document.getElementById('btn-sim').addEventListener('click', compare);
  ['sim-a', 'sim-b'].forEach(id => {
    document.getElementById(id).addEventListener('keydown', e => {
      if (e.key === 'Enter') compare();
    });
  });

  // MOL download
  document.getElementById('btn-download-mol').addEventListener('click', downloadMol);
  document.getElementById('btn-share-molecule').addEventListener('click', shareCurrentMolecule);

  // Guided workflow cards
  document.querySelectorAll('.quickstart-card').forEach(btn => {
    btn.addEventListener('click', () => launchWorkflow(btn.dataset.workflow));
  });

  // Atom index toggle
  document.getElementById('toggle-indices').addEventListener('click', () => {
    showIndices = !showIndices;
    document.getElementById('toggle-indices').classList.toggle('active', showIndices);
    if (lastMol) renderSvg(lastMol);
  });

  // SMARTS
  document.getElementById('btn-smarts-hi').addEventListener('click', applySmarts);
  document.getElementById('btn-smarts-clear').addEventListener('click', clearSmarts);
  document.getElementById('smarts-input').addEventListener('keydown', e => {
    if (e.key === 'Enter') applySmarts();
  });

  // Reaction presets
  document.querySelectorAll('.preset-rxn').forEach(btn => {
    btn.addEventListener('click', () => {
      document.getElementById('rxn-smirks').value = btn.dataset.smirks;
      document.getElementById('rxn-reactants').value = btn.dataset.reactants;
    });
  });

  // Run reaction (SMIRKS transformer)
  document.getElementById('btn-rxn').addEventListener('click', runReaction);

  // Reaction equation presets
  document.querySelectorAll('.preset-rxn-eq').forEach(btn => {
    btn.addEventListener('click', () => {
      document.getElementById('rxn-eq-input').value = btn.dataset.rxn;
    });
  });

  // Draw reaction equation
  document.getElementById('btn-rxn-eq').addEventListener('click', runReactionEquation);
  document.getElementById('rxn-eq-input').addEventListener('keydown', e => {
    if (e.key === 'Enter') runReactionEquation();
  });

  // Stereo isomers
  document.getElementById('btn-stereo').addEventListener('click', showStereoIsomers);

  // MMP
  document.getElementById('btn-mmp').addEventListener('click', runMmp);
  document.querySelectorAll('.preset-mmp').forEach(btn => {
    btn.addEventListener('click', () => {
      if (btn.id === 'mmp-preset-alkylbenzenes') {
        document.getElementById('mmp-input').value = MMP_PRESETS.alkylbenzenes;
      } else if (btn.id === 'mmp-preset-nsaids') {
        document.getElementById('mmp-input').value = MMP_PRESETS.nsaids;
      }
    });
  });

  // R-group
  document.getElementById('btn-rgroup').addEventListener('click', runRgroup);
  document.querySelectorAll('.preset-rgroup').forEach(btn => {
    btn.addEventListener('click', () => {
      const preset = btn.id === 'rgroup-preset-benzene' ? RGROUP_PRESETS.benzene : RGROUP_PRESETS.pyridine;
      document.getElementById('rgroup-input').value = preset.smiles;
      document.getElementById('rgroup-core').value = preset.core;
    });
  });

  // SDF / MOL upload
  document.getElementById('btn-sdf-load').addEventListener('click', loadSdf);

  // Formats & Materials
  initFormatsTab();

  // 3D style buttons
  document.querySelectorAll('[data-style]').forEach(btn => {
    btn.addEventListener('click', () => {
      document.querySelectorAll('[data-style]').forEach(b => b.classList.remove('active'));
      btn.classList.add('active');
      current3dStyle = btn.dataset.style;
      if (viewer3d && lastMol) apply3dStyle();
    });
  });

  // Dynamics: Coulomb energy
  document.getElementById('btn-coulomb').addEventListener('click', async () => {
    if (!lastMol) { showError('coulomb-result', 'No molecule loaded'); return; }
    try {
      const result = window.chemWasm.coulomb_energy_json(lastMol);
      const data = JSON.parse(result);
      document.getElementById('coulomb-output').textContent =
        `Coulomb Energy: ${data.coulomb_energy.toFixed(4)} ${data.unit}`;
      document.getElementById('coulomb-result').style.display = 'block';

      // Show energy breakdown
      showEnergyBreakdown({
        coulomb: data.coulomb_energy,
        total: data.coulomb_energy,
        type: 'Coulomb (Direct)'
      });
    } catch (e) {
      showError('coulomb-result', `Error: ${e.message}`);
    }
  });

  // Dynamics: MD simulation
  document.getElementById('btn-md').addEventListener('click', async () => {
    if (!lastMol) { showError('md-result', 'No molecule loaded'); return; }
    try {
      const steps = parseInt(document.getElementById('md-steps').value) || 100;
      const result = window.chemWasm.run_md_json(lastMol, steps, 300.0);
      const data = JSON.parse(result);
      const frames = data.frames || [];
      document.getElementById('md-frames').textContent = frames.length;

      let output = `MD Simulation Results (${frames.length} frames):\n\n`;
      frames.forEach((f, i) => {
        output += `Frame ${i+1} (Step ${f.step}): E_pot=${f.potential.toFixed(2)}, E_kin=${f.kinetic.toFixed(2)}, T=${f.temp.toFixed(1)}K\n`;
      });
      document.getElementById('md-output').textContent = output;
      document.getElementById('md-result').style.display = 'block';

      // Show energy breakdown for final frame
      if (frames.length > 0) {
        const finalFrame = frames[frames.length - 1];
        showEnergyBreakdown({
          potential: finalFrame.potential,
          kinetic: finalFrame.kinetic,
          total: finalFrame.potential + finalFrame.kinetic,
          type: 'MD (Final Frame)'
        });
      }
    } catch (e) {
      showError('md-result', `Error: ${e.message}`);
    }
  });

  // Helper: Display energy breakdown chart and table
  function showEnergyBreakdown(energyData) {
    const energyChart = document.getElementById('energy-chart');
    const energyHint = document.getElementById('energy-hint');
    const tableBody = document.getElementById('energy-table-body');

    // Build component list from energyData
    const components = [];
    if (energyData.potential !== undefined) components.push(['Potential', energyData.potential]);
    if (energyData.kinetic !== undefined) components.push(['Kinetic', energyData.kinetic]);
    if (energyData.coulomb !== undefined) components.push(['Coulomb', energyData.coulomb]);
    if (energyData.bond !== undefined) components.push(['Bond', energyData.bond]);
    if (energyData.angle !== undefined) components.push(['Angle', energyData.angle]);
    if (energyData.vdw !== undefined) components.push(['VDW', energyData.vdw]);

    // Populate table
    tableBody.textContent = '';
    const total = energyData.total || Math.abs(components.reduce((s, [_, e]) => s + e, 0));
    components.forEach(([name, energy]) => {
      const row = document.createElement('tr');
      row.style.borderBottom = '1px solid var(--border)';
      const pct = Math.abs(energy / total) * 100;
      row.innerHTML = `<td style="padding: 0.5rem;">${name}</td><td style="text-align: right; padding: 0.5rem;">${energy.toFixed(4)}</td><td style="text-align: right; padding: 0.5rem;">${pct.toFixed(1)}%</td>`;
      tableBody.appendChild(row);
    });
    const totalRow = document.createElement('tr');
    totalRow.style.fontWeight = '600';
    totalRow.innerHTML = `<td style="padding: 0.5rem;">Total</td><td style="text-align: right; padding: 0.5rem;">${total.toFixed(4)}</td><td style="text-align: right; padding: 0.5rem;">100%</td>`;
    tableBody.appendChild(totalRow);

    energyChart.style.display = 'block';
    energyHint.style.display = 'none';

    // Draw simple bar chart on canvas
    const canvas = document.getElementById('energy-canvas');
    if (canvas && canvas.getContext) {
      drawEnergyChart(canvas, components, total);
    }
  }

  function drawEnergyChart(canvas, components, total) {
    const ctx = canvas.getContext('2d');
    if (!ctx) return;

    const padding = 40;
    const width = canvas.width - 2 * padding;
    const height = canvas.height - 2 * padding;
    const barWidth = width / components.length / 1.5;
    const maxEnergy = Math.max(...components.map(([_, e]) => Math.abs(e)));

    ctx.clearRect(0, 0, canvas.width, canvas.height);
    ctx.fillStyle = 'var(--text)';
    ctx.font = '12px var(--font-mono)';
    ctx.textAlign = 'center';

    const colors = ['#4f8ef7', '#34d399', '#f87171', '#fbbf24', '#a78bfa', '#06b6d4'];
    components.forEach(([name, energy], i) => {
      const barHeight = (Math.abs(energy) / maxEnergy) * height;
      const x = padding + i * (width / components.length) + (width / components.length - barWidth) / 2;
      const y = padding + height - barHeight;

      ctx.fillStyle = colors[i % colors.length];
      ctx.fillRect(x, y, barWidth, barHeight);

      ctx.fillStyle = 'var(--muted)';
      ctx.fillText(name, x + barWidth / 2, canvas.height - padding + 20);

      ctx.fillStyle = 'var(--text)';
      ctx.fillText(energy.toFixed(2), x + barWidth / 2, y - 5);
    });

    // Axes
    ctx.strokeStyle = 'var(--border)';
    ctx.beginPath();
    ctx.moveTo(padding, padding);
    ctx.lineTo(padding, padding + height);
    ctx.lineTo(padding + width, padding + height);
    ctx.stroke();
  }

  // Dynamics: DREIDING minimize
  document.getElementById('btn-minimize').addEventListener('click', async () => {
    if (!lastMol) { showError('minimize-result', 'No molecule loaded'); return; }
    try {
      const pdbStr = window.chemWasm.minimize_dreiding_json(lastMol);
      const resultDiv = document.getElementById('minimize-result');
      resultDiv.textContent = '';
      const p1 = document.createElement('p');
      p1.textContent = '✓ Optimization complete! Geometry minimized with DREIDING force field.';
      const p2 = document.createElement('p');
      p2.textContent = 'Minimized structure (PDB format):';
      const pre = document.createElement('pre');
      pre.style.maxHeight = '200px';
      pre.style.overflowY = 'auto';
      pre.style.fontSize = '0.8rem';
      pre.textContent = pdbStr.split('\n').slice(0, 10).join('\n') + '\n...';
      resultDiv.appendChild(p1);
      resultDiv.appendChild(p2);
      resultDiv.appendChild(pre);
      resultDiv.style.display = 'block';
    } catch (e) {
      showError('minimize-result', `Error: ${e.message}`);
    }
  });

  // ── MMFF94 Force Field Optimizer ─────────────────────────────────────────────
  function makeTd(text, style) {
    const td = document.createElement('td');
    td.textContent = text;
    if (style) Object.assign(td.style, style);
    return td;
  }
  function makeBar(pct, color) {
    const td = document.createElement('td');
    const bar = document.createElement('div');
    Object.assign(bar.style, {background:color,height:'10px',borderRadius:'3px',width:Math.max(0,pct)+'%',minWidth:'2px'});
    td.appendChild(bar);
    return td;
  }

  document.getElementById('btn-mmff94-opt').addEventListener('click', async () => {
    if (!lastMol) { alert('No molecule loaded'); return; }
    if (!window.chemWasm) { alert('WASM not ready'); return; }
    const algo = document.getElementById('mmff94-algo').value;
    const maxIter = parseInt(document.getElementById('mmff94-max-iter').value, 10) || 300;
    try {
      const fn = algo === 'lbfgs' ? 'minimize_mmff94_lbfgs_json' : 'minimize_mmff94_json';
      const r = JSON.parse(window.chemWasm[fn](lastMol, maxIter));
      if (r.error) { alert('MMFF94 error: ' + r.error); return; }
      document.getElementById('mmff94-energy').textContent = r.energy.toFixed(4);
      document.getElementById('mmff94-rmsd').textContent   = r.rmsd.toFixed(3);
      document.getElementById('mmff94-iters').textContent  = r.iterations;
      document.getElementById('mmff94-conv').textContent   = r.converged ? 'Yes' : 'No';

      const bd = JSON.parse(window.chemWasm.mmff94_energy_breakdown_json(lastMol));
      const tbody = document.getElementById('mmff94-breakdown-body');
      tbody.textContent = '';
      const terms = [['Bond', bd.bond], ['Angle', bd.angle], ['Torsion', bd.torsion], ['vdW', bd.vdw], ['Electrostatic', bd.elec]];
      const posSum = terms.reduce((s, [, v]) => s + Math.max(0, v), 0) || 1;
      const colors = ['#4c9be8','#f59e42','#6bcb77','#e84c4c','#a78bfa'];
      terms.forEach(([name, val], idx) => {
        const tr = document.createElement('tr');
        tr.style.borderBottom = '1px solid var(--border)';
        const pct = (Math.abs(val) / Math.abs(bd.total || 1) * 100).toFixed(1);
        tr.appendChild(makeTd(name, {padding:'.35rem .6rem'}));
        tr.appendChild(makeTd(val.toFixed(4), {textAlign:'right',padding:'.35rem .6rem'}));
        tr.appendChild(makeTd(pct+'%', {textAlign:'right',padding:'.35rem .6rem'}));
        tr.appendChild(makeBar(Math.max(0, val / posSum * 100), colors[idx]));
        tbody.appendChild(tr);
      });
      const tr2 = document.createElement('tr');
      tr2.style.fontWeight = '600';
      tr2.appendChild(makeTd('Total', {padding:'.35rem .6rem'}));
      tr2.appendChild(makeTd(bd.total.toFixed(4), {textAlign:'right',padding:'.35rem .6rem'}));
      tr2.appendChild(document.createElement('td'));
      tr2.appendChild(document.createElement('td'));
      tbody.appendChild(tr2);
      document.getElementById('mmff94-opt-result').style.display = 'block';
    } catch (e) {
      alert('Error: ' + e.message);
    }
  });

  // ── Force Field Comparison ────────────────────────────────────────────────────
  document.getElementById('btn-ff-compare').addEventListener('click', async () => {
    if (!lastMol) { alert('No molecule loaded'); return; }
    if (!window.chemWasm) { alert('WASM not ready'); return; }
    try {
      const t0 = performance.now();
      window.chemWasm.minimize_dreiding_json(lastMol);
      const drMs = (performance.now() - t0).toFixed(0);

      const t1 = performance.now();
      const mr = JSON.parse(window.chemWasm.minimize_mmff94_lbfgs_json(lastMol, 300));
      const mmMs = (performance.now() - t1).toFixed(0);

      const tbody = document.getElementById('ff-compare-body');
      tbody.textContent = '';
      [[`DREIDING`, `—`, `—`, `~${drMs} ms`],
       [`MMFF94 L-BFGS`, mr.error ? 'error' : mr.energy.toFixed(4), mr.error ? '—' : mr.rmsd.toFixed(3), `${mr.iterations||'—'} iters / ${mmMs} ms`]
      ].forEach(([ff, e, rmsd, extra]) => {
        const tr = document.createElement('tr');
        tr.style.borderBottom = '1px solid var(--border)';
        const b = document.createElement('strong'); b.textContent = ff;
        const td0 = document.createElement('td'); td0.style.padding = '.4rem .6rem'; td0.appendChild(b);
        tr.appendChild(td0);
        tr.appendChild(makeTd(e,    {textAlign:'right',padding:'.4rem .6rem'}));
        tr.appendChild(makeTd(rmsd, {textAlign:'right',padding:'.4rem .6rem'}));
        tr.appendChild(makeTd(extra,{textAlign:'right',padding:'.4rem .6rem',color:'var(--muted)',fontSize:'.82rem'}));
        tbody.appendChild(tr);
      });
      document.getElementById('ff-compare-result').style.display = 'block';
    } catch (e) {
      alert('Error: ' + e.message);
    }
  });

  // ── Torsion Scan ─────────────────────────────────────────────────────────────
  document.getElementById('btn-torsion-scan').addEventListener('click', async () => {
    if (!lastMol) { alert('No molecule loaded'); return; }
    if (!window.chemWasm) { alert('WASM not ready'); return; }
    const si = parseInt(document.getElementById('scan-i').value, 10);
    const sj = parseInt(document.getElementById('scan-j').value, 10);
    const sk = parseInt(document.getElementById('scan-k').value, 10);
    const sl = parseInt(document.getElementById('scan-l').value, 10);
    const steps = parseInt(document.getElementById('scan-steps').value, 10) || 36;
    try {
      const result = JSON.parse(window.chemWasm.torsion_scan_json(lastMol, si, sj, sk, sl, steps));
      if (result.error) { alert('Scan error: ' + result.error); return; }
      const energies = result.map(p => p.energy);
      const minE = Math.min(...energies), maxE = Math.max(...energies);
      const range = (maxE - minE) || 1;
      const W = 420, H = 180, pl = 55, pt = 15, pr = 15, pb = 30;
      const cW = W - pl - pr, cH = H - pt - pb;
      const toX = a => pl + (a / 360) * cW;
      const toY = e => pt + (1 - (e - minE) / range) * cH;
      const pts = result.map(p => `${toX(p.angle).toFixed(1)},${toY(p.energy).toFixed(1)}`).join(' ');
      const minIdx = energies.indexOf(minE), maxIdx = energies.indexOf(maxE);
      document.getElementById('scan-min').textContent = minE.toFixed(3);
      document.getElementById('scan-min-angle').textContent = result[minIdx].angle.toFixed(1);
      document.getElementById('scan-max').textContent = maxE.toFixed(3);
      document.getElementById('scan-max-angle').textContent = result[maxIdx].angle.toFixed(1);
      const ns = 'http://www.w3.org/2000/svg';
      const svg = document.getElementById('torsion-plot');
      while (svg.firstChild) svg.removeChild(svg.firstChild);
      function svgEl(tag, attrs) {
        const el = document.createElementNS(ns, tag);
        Object.entries(attrs).forEach(([k, v]) => el.setAttribute(k, v));
        return el;
      }
      svg.appendChild(svgEl('line', {x1:pl,y1:pt,x2:pl,y2:pt+cH,stroke:'var(--border)','stroke-width':1}));
      svg.appendChild(svgEl('line', {x1:pl,y1:pt+cH,x2:pl+cW,y2:pt+cH,stroke:'var(--border)','stroke-width':1}));
      [0,90,180,270,360].forEach(a => {
        const t = svgEl('text', {x:toX(a).toFixed(1),y:pt+cH+12,'text-anchor':'middle','font-size':9,fill:'var(--muted)'}); t.textContent = a; svg.appendChild(t);
      });
      const tmin = svgEl('text',{x:pl-4,y:pt+cH+4,'text-anchor':'end','font-size':9,fill:'var(--muted)'}); tmin.textContent = minE.toFixed(1); svg.appendChild(tmin);
      const tmax = svgEl('text',{x:pl-4,y:pt+4,'text-anchor':'end','font-size':9,fill:'var(--muted)'}); tmax.textContent = maxE.toFixed(1); svg.appendChild(tmax);
      svg.appendChild(svgEl('polyline', {points:pts,fill:'none',stroke:'#4c9be8','stroke-width':2,'stroke-linejoin':'round'}));
      svg.appendChild(svgEl('circle', {cx:toX(result[minIdx].angle).toFixed(1),cy:toY(minE).toFixed(1),r:4,fill:'#6bcb77'}));
      svg.appendChild(svgEl('circle', {cx:toX(result[maxIdx].angle).toFixed(1),cy:toY(maxE).toFixed(1),r:4,fill:'#e84c4c'}));
      document.getElementById('torsion-scan-result').style.display = 'block';
    } catch (e) {
      alert('Error: ' + e.message);
    }
  });

  // ── MMFF94 Charge Map (2D tab) ────────────────────────────────────────────────
  let chargeMapActive = false;
  document.getElementById('btn-charge-map').addEventListener('click', async () => {
    if (!lastMol) { alert('No molecule loaded'); return; }
    if (!window.chemWasm) { alert('WASM not ready'); return; }
    chargeMapActive = !chargeMapActive;
    const btn = document.getElementById('btn-charge-map');
    btn.style.background = chargeMapActive ? '#4c9be8' : '';
    btn.style.color = chargeMapActive ? '#fff' : '';
    if (!chargeMapActive) { if (lastSmiles) document.getElementById('btn-calc').click(); return; }
    try {
      const qJson = JSON.parse(window.chemWasm.mmff94_partial_charges_json(lastMol));
      if (qJson.error) { alert('Charge error: ' + qJson.error); chargeMapActive = false; btn.style.background = ''; return; }
      const charges = qJson.charges;
      function chargeColor(q) {
        const v = Math.max(-1, Math.min(1, q / 0.5));
        if (v < 0) {
          const g = Math.round(220 * (1 + v));
          return 'rgb(220,' + g + ',' + g + ')';
        }
        const rb = Math.round(220 * (1 - v));
        return 'rgb(' + rb + ',' + rb + ',220)';
      }
      // Build atom_color_map JSON and request highlighted SVG
      const colorMap = {};
      charges.forEach((q, idx) => { colorMap[idx] = chargeColor(q); });
      const smiles = lastSmiles || document.getElementById('smiles-input').value;
      const colorMapJson = JSON.stringify(colorMap);
      let svgStr;
      if (window.chemWasm.depict_svg_with_atom_colors_json) {
        svgStr = window.chemWasm.depict_svg_with_atom_colors_json(lastMol, colorMapJson);
      } else {
        const opts = JSON.stringify({atomColorMap: colorMap});
        svgStr = window.chemWasm.depict_svg_json
          ? window.chemWasm.depict_svg_json(lastMol, opts)
          : window.chemWasm.depict_svg(lastMol);
      }
      const svgArea = document.getElementById('svg-area');
      const existing = svgArea.querySelector('svg');
      if (existing && svgStr && svgStr.startsWith('<svg')) {
        const tmp = document.createElement('div');
        tmp.textContent = svgStr;
        // Parse SVG safely via DOMParser
        const parser = new DOMParser();
        const parsed = parser.parseFromString(svgStr, 'image/svg+xml');
        const newSvg = parsed.querySelector('svg');
        if (newSvg) svgArea.replaceChild(document.adoptNode(newSvg), existing);
      }
    } catch (e) {
      chargeMapActive = false;
      document.getElementById('btn-charge-map').style.background = '';
      alert('Error: ' + e.message);
    }
  });

  // Report tab: workflow API integration
  document.getElementById('btn-report').addEventListener('click', async () => {
    const smiles = document.getElementById('report-smiles').value.trim();
    if (!smiles) { showError('error-report', 'SMILES required'); return; }
    if (!window.chemWasm) { showError('error-report', 'WASM not ready'); return; }

    try {
      const json = window.chemWasm.molecule_report_json(smiles);
      const report = JSON.parse(json);
      lastReportJson = json;

      document.getElementById('report-canonical').textContent = `Canonical: ${report.canonical_smiles}`;
      populateReportDescriptors(report.descriptors);
      populateReportFilters(report.filters);
      document.getElementById('report-output').style.display = 'block';
      document.getElementById('error-report').classList.add('hidden');
    } catch (e) {
      showError('error-report', `Error: ${e.message}`);
    }
  });

  document.getElementById('btn-report-json').addEventListener('click', () => {
    if (!lastReportJson) return;
    const status = document.getElementById('report-json-status');
    const copied = () => {
      status.textContent = t('report.copied');
      setTimeout(() => { status.textContent = ''; }, 2200);
    };
    const writePromise = navigator.clipboard?.writeText(JSON.stringify(JSON.parse(lastReportJson), null, 2));
    if (writePromise) {
      writePromise.then(copied).catch(() => { status.textContent = lastReportJson; });
    } else {
      status.textContent = lastReportJson;
    }
  });

  // ── Copy code buttons ────────────────────────────────────────────────────
  document.getElementById('btn-copy-py').addEventListener('click', () => {
    flashCopy(document.getElementById('code-wrap')._pySnippet || '');
  });
  document.getElementById('btn-copy-js').addEventListener('click', () => {
    flashCopy(document.getElementById('code-wrap')._jsSnippet || '');
  });

  function populateReportDescriptors(descriptors) {
    const tbody = document.getElementById('report-desc-tbody');
    tbody.textContent = '';
    const rows = [
      ['MW', descriptors.molecular_weight.toFixed(2)],
      ['LogP', descriptors.logp.toFixed(2)],
      ['TPSA', descriptors.tpsa.toFixed(2)],
      ['HBA', descriptors.hba],
      ['HBD', descriptors.hbd],
      ['Rot.Bonds', descriptors.rotatable_bonds],
      ['QED', descriptors.qed.toFixed(3)],
      ['SA Score', descriptors.sa_score.toFixed(2)],
    ];
    rows.forEach(([prop, val]) => {
      const tr = document.createElement('tr');
      const td1 = document.createElement('td');
      td1.textContent = prop;
      const td2 = document.createElement('td');
      td2.textContent = val;
      td2.style.textAlign = 'right';
      tr.appendChild(td1);
      tr.appendChild(td2);
      tbody.appendChild(tr);
    });
  }

  function populateReportFilters(filters) {
    const grid = document.getElementById('report-filters-grid');
    grid.textContent = '';
    const filterList = [
      ['Lipinski', filters.lipinski_passes],
      ['Veber', filters.veber_passes],
      ['Egan', filters.egan_passes],
      ['Ghose', filters.ghose_passes],
      ['REOS', filters.reos_passes],
      ['PAINS', filters.pains_passes],
    ];
    filterList.forEach(([name, passes]) => {
      const div = document.createElement('div');
      div.style.padding = '0.75rem';
      div.style.border = `1px solid ${passes ? 'var(--pass)' : 'var(--fail)'}`;
      div.style.borderRadius = 'var(--radius)';
      div.style.backgroundColor = passes ? 'color-mix(in srgb,var(--pass) 8%,transparent)' : 'color-mix(in srgb,var(--fail) 8%,transparent)';
      div.style.color = passes ? 'var(--pass)' : 'var(--fail)';
      div.style.textAlign = 'center';
      div.style.fontSize = '0.85rem';
      const strong = document.createElement('strong');
      strong.textContent = name;
      div.appendChild(strong);
      div.appendChild(document.createElement('br'));
      const status = document.createTextNode(passes ? '✓ Pass' : '✗ Fail');
      div.appendChild(status);
      grid.appendChild(div);
    });
  }
})();
