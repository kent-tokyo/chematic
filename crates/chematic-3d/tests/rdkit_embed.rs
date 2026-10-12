//! `rdkit_embed_molecule` against RDKit 2026.09.1 ETKDGv3.
//!
//! The historical 2026.03 coordinate arrays below retain atom-count coverage;
//! current numeric parity is pinned by the simple 2026.09 oracle and the
//! versioned corpus-scale A6 comparison lane.

use chematic_3d::rdkit_embed::{
    RdkitEmbedError, RdkitEmbedOptions, rdkit_embed_initial_molecule, rdkit_embed_molecule,
    rdkit_embed_multiple_confs, rdkit_embed_multiple_confs_mmff,
    rdkit_embed_multiple_confs_mmff_adaptive, rdkit_embed_multiple_confs_mmff_diagnostics,
    rdkit_embed_multiple_confs_mmff_gradient_rescue, rdkit_embed_multiple_confs_mmff_rescue,
};

#[allow(clippy::approx_constant, clippy::excessive_precision)]
const CASES: &[(&str, &[[f64; 3]])] = &[
    (
        "CCO",
        &[
            [
                -0.9534551842756148,
                0.04780421609619026,
                0.04249937663911896,
            ],
            [
                0.4879186876831028,
                -0.32153931807429664,
                -0.18904973635963182,
            ],
            [1.276977526346871, 0.3247836784511813, 0.7376732425087271],
            [-1.3570126123500212, 0.7339751954162034, -0.7332225613847179],
            [
                -1.5933209403261652,
                -0.8625698312047165,
                0.021759475665728743,
            ],
            [-1.0710312763466439, 0.5964393608027803, 1.011162760505495],
            [
                0.7874060542693289,
                -0.14060868171000282,
                -1.2357743239018328,
            ],
            [
                0.5900278787458185,
                -1.4261096608429211,
                -0.04403048576654172,
            ],
            [1.832489866253323, 1.047825041065579, 0.38898225209365683],
        ],
    ),
    (
        "C[C@H](N)C(=O)O",
        &[
            [
                -0.8592407634599415,
                1.1514661472637542,
                -0.16522593775029595,
            ],
            [
                -0.1853754145711437,
                -0.15987068859904846,
                -0.5099590208150512,
            ],
            [-0.802105522995388, -1.1734690014909304, 0.28318342058733675],
            [
                1.2541597265949056,
                -0.13057261963056138,
                -0.1652195934721842,
            ],
            [1.972011877866091, 0.7167166775677084, -0.7762230953089931],
            [1.821236640219406, -0.9882817051797496, 0.7878718150743558],
            [-0.2750021524916736, 1.5837702964536673, 0.6895632730001544],
            [-0.9005780761432235, 1.8285273117575287, -1.0365760821155205],
            [-1.8862461561959238, 0.996256849997651, 0.2550251212858728],
            [
                -0.33699872662889885,
                -0.34379839275665236,
                -1.583884215673667,
            ],
            [-1.3606582832409595, -0.8755071898300736, 1.0831481019512659],
            [
                -1.0289647171341654,
                -2.0490815471739134,
                -0.18946831775330764,
            ],
            [2.587761568180817, -0.5561561383794941, 1.3277645309899992],
        ],
    ),
    (
        "C/C=C/C(=O)NC",
        &[
            [
                2.9500828488785844,
                -0.08627105255554705,
                -0.3036799472277232,
            ],
            [1.47502409975286, -0.05564872623557736, -0.298038277940883],
            [0.7462542835276481, 0.10666780396736188, 0.7863367537462322],
            [-0.6993261572608513, 0.1396219382349811, 0.8134940563823534],
            [-1.335975703119055, 0.30178013891356376, 1.9103455748880995],
            [
                -1.4920149575296606,
                -0.0021786907326252565,
                -0.3463814166473156,
            ],
            [
                -2.937443464600719,
                0.05061140387758834,
                -0.17326275093631568,
            ],
            [3.2508618255964032, -1.0397077261928445, -0.8298914774544935],
            [3.4037463378015125, 0.8020278093879532, -0.8213283732031329],
            [3.3952986808220786, -0.18139819140794777, 0.6948143472885984],
            [0.9730111636912291, -0.1748948678226098, -1.2488998364636617],
            [1.287474088190356, 0.22280894203568685, 1.7203052793556681],
            [-1.0657691238282314, -0.137560338817556, -1.2779628579355562],
            [
                -3.2384953489368313,
                -0.24704448339444945,
                0.8556874221799673,
            ],
            [
                -3.3577188132776814,
                -0.7330824783561386,
                -0.8426667399030011,
            ],
            [
                -3.3550097597077464,
                1.0342685190981658,
                -0.48458400701503496,
            ],
        ],
    ),
    (
        "c1ccc2[nH]ccc2c1",
        &[
            [
                -2.0578332763082767,
                -0.7035875397948492,
                -0.11842283437971346,
            ],
            [-2.076928263828804, 0.6597464844272636, 0.15496810694940372],
            [-0.8782412090611019, 1.3598627343420218, 0.2821205895378317],
            [0.31445796519695973, 0.6828345290456674, 0.13338152126831698],
            [1.5903548196076058, 1.1051462053315486, 0.20400873747041373],
            [
                2.398254675305616,
                0.035174186192287606,
                -0.01925387086388418,
            ],
            [
                1.6356548445385486,
                -1.0946597765706132,
                -0.23727987633595646,
            ],
            [0.3157543296318745, -0.6692260944787213, -0.1375587074188902],
            [
                -0.8501238421555126,
                -1.384313127972235,
                -0.26806609740714565,
            ],
            [-2.991786689972323, -1.257367988810752, -0.2191293473513957],
            [-3.0157985547530464, 1.1870258419120763, 0.270904559999922],
            [-0.9482951310254392, 2.417916505189476, 0.4949104802526355],
            [1.8910231586324955, 2.0967752942338813, 0.39941527965355417],
            [3.484430532656664, 0.0877707682895754, -0.020629946980938042],
            [2.021820305262349, -2.0835353564198065, -0.4396569386144659],
            [
                -0.8327436637276501,
                -2.4395626649168154,
                -0.4797116557796753,
            ],
        ],
    ),
    (
        "OC1CCC2CCCCC2C1",
        &[
            [3.2954199085382236, -0.9510233257776786, 1.1165690314026173],
            [2.3275410277274395, -0.4469885977256517, 0.272409883144957],
            [2.195442103942001, 1.0443150807151513, 0.32750083997390317],
            [0.8151597007067561, 1.5584009545964457, 0.32625800784537423],
            [-0.176198173839366, 0.7904135998132864, -0.5610052469885649],
            [-1.5312405000441, 1.24445369963706, -0.16694937007345398],
            [-2.608704767356766, 0.22417464108649968, -0.1273737317083551],
            [
                -2.3781778278248997,
                -1.0066056251798055,
                -0.9418567925495308,
            ],
            [
                -1.0603539380738596,
                -1.5847446138988595,
                -0.5364019820703144,
            ],
            [
                0.08931106563223179,
                -0.6573581361710933,
                -0.5991076415353017,
            ],
            [0.9752926024055132, -1.0115687675617018, 0.5999058029669032],
            [4.124867257556764, -0.41938650603415534, 0.9505135875924029],
            [2.5495175505975727, -0.7198253724268069, -0.7818426617925159],
            [2.8118925076346275, 1.5754512894769754, -0.44722699010636097],
            [2.659499134931002, 1.349728381740679, 1.3007703679268323],
            [0.3259865672714306, 1.604311817609698, 1.342608200218282],
            [
                0.8229710643221996,
                2.6260808259935313,
                -0.018653508459217063,
            ],
            [0.04069020565451612, 1.210742376986471, -1.5895729628038169],
            [-1.8904284119541452, 2.090016133107567, -0.8283724942949106],
            [-1.541842629642985, 1.7310931568225298, 0.8436317758254865],
            [
                -2.8283268743661742,
                -0.11858609412208881,
                0.9192932059587355,
            ],
            [-3.587533229767922, 0.6479295050851678, -0.4788010412156998],
            [
                -3.1562395079828183,
                -1.7508506578039986,
                -0.6261940677753618,
            ],
            [-2.496404435380898, -0.840428280539962, -2.0230845224798473],
            [-1.2104079642451766, -2.0435674372542865, 0.4783497311279739],
            [
                -0.8898689617342375,
                -2.4748155771134392,
                -1.2001332064395536,
            ],
            [0.7306190031583832, -0.9323981631382868, -1.4796782716134183],
            [0.5408766454456073, -0.6266311330925628, 1.5401947509640248],
            [1.0506408766890507, -2.112333174830769, 0.6465492781986695],
        ],
    ),
];

// First-atom snapshot from RDKit 2026.09.1 ETKDGv3 with randomSeed=42.
// Broader parity is gated by the versioned A6 bounds/quality corpus because a
// last-bit bounds change can legitimately move a flexible molecule to another
// local minimum while preserving its chemical constraints.
const ETHANOL_FIRST_ATOM_2026_09_1: [f64; 3] = [
    -0.9639577069927618,
    0.015102231353043739,
    0.0410206427653118,
];

const ETHANOL_INITIAL_2026_09_1: [[f64; 3]; 9] = [
    [-0.7500043088651186, 0.4084390159458829, 0.06848496479484886],
    [
        0.44564780848200714,
        -0.30357044432024505,
        -0.2925815152031927,
    ],
    [1.1754601191726461, -0.458924916623023, 1.1664313728171822],
    [-1.2354609528915097, 1.0831648725953975, -0.5316473904156442],
    [
        -1.5724182037685812,
        -1.0383486994621525,
        0.45351405445069587,
    ],
    [
        -1.2999968573900929,
        0.40241491422299236,
        0.40181982062149446,
    ],
    [0.5525517516037974, -0.5577005629983048, -0.9232674186573715],
    [0.9120715852990083, -0.80779822184343, -0.5889154542504689],
    [1.7721490583578425, 1.2723240424828839, 0.24616156584245455],
];

#[test]
fn initial_metric_embedding_matches_rdkit_2026_09_boundary() {
    let mh = chematic_chem::add_hydrogens(&chematic_smiles::parse("CCO").unwrap());
    let got = rdkit_embed_initial_molecule(&mh, &RdkitEmbedOptions::default()).unwrap();
    assert_eq!(got.len(), ETHANOL_INITIAL_2026_09_1.len());
    for (actual, expected) in got.iter().zip(ETHANOL_INITIAL_2026_09_1) {
        for coordinate in 0..3 {
            assert!((actual[coordinate] - expected[coordinate]).abs() <= 5e-4);
        }
    }
    assert_eq!(
        got,
        rdkit_embed_initial_molecule(&mh, &RdkitEmbedOptions::default()).unwrap()
    );
}

/// Bit-identical on the recorded Linux x86-64 comparison lane; elsewhere the
/// platform libm's last-bit differences are allowed to accumulate through
/// the minimizations, so coordinates are compared within 5e-4 Å. The wider
/// bound is required on macOS arm64, where the accumulated libm difference in
/// this fixture reaches about 1.3e-4 Å.
#[test]
#[allow(unexpected_cfgs)]
fn etkdg_2026_09_seed_is_deterministic_and_ethanol_matches_oracle() {
    // LLVM coverage instrumentation perturbs the accumulated floating-point
    // result by a few ULPs even on the pinned Linux x86-64 lane. Keep exact
    // comparison for ordinary release/test builds and use the documented
    // cross-platform tolerance for instrumented coverage builds.
    let exact = cfg!(all(
        target_os = "linux",
        target_arch = "x86_64",
        not(coverage)
    ));
    for &(smi, legacy_expected) in CASES {
        let mh = chematic_chem::add_hydrogens(&chematic_smiles::parse(smi).unwrap());
        let got = rdkit_embed_molecule(&mh, &RdkitEmbedOptions::default())
            .unwrap_or_else(|e| panic!("{smi}: {e}"));
        assert_eq!(got.len(), legacy_expected.len(), "{smi}");
        let repeated = rdkit_embed_molecule(&mh, &RdkitEmbedOptions::default())
            .unwrap_or_else(|e| panic!("{smi}: {e}"));
        assert_eq!(
            got, repeated,
            "{smi}: fixed-seed embedding changed between calls"
        );
        if smi == "CCO" {
            for coordinate in 0..3 {
                if exact {
                    assert_eq!(got[0][coordinate], ETHANOL_FIRST_ATOM_2026_09_1[coordinate]);
                } else {
                    assert!(
                        (got[0][coordinate] - ETHANOL_FIRST_ATOM_2026_09_1[coordinate]).abs()
                            <= 5e-4,
                        "{smi}: {:?} vs RDKit first atom {ETHANOL_FIRST_ATOM_2026_09_1:?}",
                        got[0]
                    );
                }
            }
        }
    }
}

#[test]
fn negative_seed_is_refused_and_fragments_are_embedded() {
    let mh = chematic_chem::add_hydrogens(&chematic_smiles::parse("CCO").unwrap());
    let opts = RdkitEmbedOptions {
        random_seed: -1,
        ..Default::default()
    };
    assert!(matches!(
        rdkit_embed_molecule(&mh, &opts),
        Err(RdkitEmbedError::Unsupported(_))
    ));
    // Multi-fragment molecules are embedded fragment by fragment, as
    // RDKit's EmbedMultipleConfs does with embedFragmentsSeparately.
    let salt = chematic_chem::add_hydrogens(&chematic_smiles::parse("CC(=O)[O-].[Na+]").unwrap());
    let coords = rdkit_embed_molecule(&salt, &RdkitEmbedOptions::default())
        .unwrap_or_else(|e| panic!("salt: {e}"));
    assert_eq!(coords.len(), salt.atom_count());
}

#[test]
fn multiple_conformers_match_rdkit_seed_progression() {
    let mh = chematic_chem::add_hydrogens(&chematic_smiles::parse("CCO").unwrap());
    let options = RdkitEmbedOptions::default();
    let conformers = rdkit_embed_multiple_confs(&mh, 3, &options).unwrap();
    assert_eq!(conformers.len(), 3);
    assert_eq!(conformers[0], rdkit_embed_molecule(&mh, &options).unwrap());

    let expected_first_atoms = [
        [
            -0.9639577069927618,
            0.015102231353043739,
            0.0410206427653118,
        ],
        [
            -0.946941193051227,
            0.20466995457287168,
            0.007322645786689314,
        ],
        [
            -0.8949265887204274,
            0.21063088448005707,
            -0.017107180966725415,
        ],
    ];
    for (conformer, expected) in conformers.iter().zip(expected_first_atoms) {
        for coordinate in 0..3 {
            assert!((conformer[0][coordinate] - expected[coordinate]).abs() <= 5.0e-4);
        }
    }

    assert!(
        rdkit_embed_multiple_confs(&mh, 0, &options)
            .unwrap()
            .is_empty()
    );
}

#[test]
fn multiple_conformer_mmff_is_deterministic_and_reuses_exact_field() {
    let mh = chematic_chem::add_hydrogens(&chematic_smiles::parse("CCO").unwrap());
    let options = RdkitEmbedOptions::default();
    let embedded = rdkit_embed_multiple_confs(&mh, 3, &options).unwrap();
    let optimized = rdkit_embed_multiple_confs_mmff(&mh, 3, &options, 200).unwrap();
    let repeated = rdkit_embed_multiple_confs_mmff(&mh, 3, &options, 200).unwrap();
    assert_eq!(optimized.len(), 3);
    for (left, right) in optimized.iter().zip(&repeated) {
        assert_eq!(left.coords, right.coords);
        assert_eq!(left.status, right.status);
        assert_eq!(left.energy, right.energy);
        assert!(left.energy.is_finite());
    }

    let field =
        chematic_ff::rdkit_mmff::RdkitMmffField::new(&mh, &embedded[0], 100.0, true).unwrap();
    let mut flat: Vec<f64> = embedded[0].iter().flatten().copied().collect();
    let (status, energy) = field.optimize(&mut flat, 200);
    assert_eq!(optimized[0].status, status);
    assert_eq!(optimized[0].energy, energy);
    assert_eq!(optimized[0].coords.as_flattened(), flat);
}

#[test]
fn multiple_conformer_mmff_diagnostics_preserve_compact_results() {
    let mh = chematic_chem::add_hydrogens(&chematic_smiles::parse("CCCO").unwrap());
    let options = RdkitEmbedOptions::default();
    let compact = rdkit_embed_multiple_confs_mmff(&mh, 4, &options, 2).unwrap();
    let detailed = rdkit_embed_multiple_confs_mmff_diagnostics(&mh, 4, &options, 2).unwrap();

    assert_eq!(compact.len(), detailed.len());
    for (compact, detailed) in compact.iter().zip(detailed) {
        assert_eq!(compact.coords, detailed.conformer.coords);
        assert_eq!(compact.status, detailed.conformer.status);
        assert_eq!(compact.energy, detailed.conformer.energy);
        assert_eq!(detailed.iterations, 2);
        assert!(detailed.final_gradient_test.is_some_and(f64::is_finite));
    }
}

#[test]
fn adaptive_mmff_retries_only_low_energy_nonconverged_candidates() {
    let mh = chematic_chem::add_hydrogens(&chematic_smiles::parse("CCCO").unwrap());
    let options = RdkitEmbedOptions::default();
    let initial = rdkit_embed_multiple_confs_mmff(&mh, 4, &options, 1).unwrap();
    let adaptive = rdkit_embed_multiple_confs_mmff_adaptive(&mh, 4, &options, 1, 200, 1).unwrap();
    let repeated = rdkit_embed_multiple_confs_mmff_adaptive(&mh, 4, &options, 1, 200, 1).unwrap();

    assert_eq!(adaptive.len(), initial.len());
    assert_eq!(adaptive.len(), repeated.len());
    assert!(initial.iter().all(|conformer| conformer.status != 0));
    assert_eq!(adaptive.iter().filter(|c| c.status == 0).count(), 1);
    assert_eq!(
        adaptive
            .iter()
            .zip(&initial)
            .filter(|(after, before)| after.coords != before.coords)
            .count(),
        1
    );
    for ((after, before), again) in adaptive.iter().zip(&initial).zip(&repeated) {
        assert!(after.energy <= before.energy);
        assert_eq!(after.coords, again.coords);
        assert_eq!(after.status, again.status);
        assert_eq!(after.energy, again.energy);
    }
}

#[test]
fn rescue_mmff_skips_retry_when_any_candidate_converged() {
    let mh = chematic_chem::add_hydrogens(&chematic_smiles::parse("CCCO").unwrap());
    let options = RdkitEmbedOptions::default();
    let initial = rdkit_embed_multiple_confs_mmff(&mh, 4, &options, 200).unwrap();
    assert!(initial.iter().any(|conformer| conformer.status == 0));
    let rescued = rdkit_embed_multiple_confs_mmff_rescue(&mh, 4, &options, 200, 1000, 1).unwrap();
    assert_eq!(rescued.len(), initial.len());
    for (after, before) in rescued.iter().zip(initial) {
        assert_eq!(after.coords, before.coords);
        assert_eq!(after.status, before.status);
        assert_eq!(after.energy, before.energy);
    }
}

#[test]
fn gradient_rescue_retries_only_candidates_above_threshold() {
    let mh = chematic_chem::add_hydrogens(&chematic_smiles::parse("CCCO").unwrap());
    let options = RdkitEmbedOptions::default();
    let initial = rdkit_embed_multiple_confs_mmff(&mh, 4, &options, 1).unwrap();
    let rescued =
        rdkit_embed_multiple_confs_mmff_gradient_rescue(&mh, 4, &options, 1, 200, 1, 0.0).unwrap();
    let skipped =
        rdkit_embed_multiple_confs_mmff_gradient_rescue(&mh, 4, &options, 1, 200, 1, f64::MAX)
            .unwrap();

    assert_eq!(rescued.iter().filter(|row| row.status == 0).count(), 1);
    assert_eq!(
        rescued
            .iter()
            .zip(&initial)
            .filter(|(after, before)| after.coords != before.coords)
            .count(),
        1
    );
    for (after, before) in skipped.iter().zip(initial) {
        assert_eq!(after.coords, before.coords);
        assert_eq!(after.status, before.status);
        assert_eq!(after.energy, before.energy);
    }
}
