fn main() {
    #[cfg(feature = "native-inchi")]
    build_native_inchi();
}

#[cfg(feature = "native-inchi")]
fn build_native_inchi() {
    let target_family = std::env::var("CARGO_CFG_TARGET_FAMILY").unwrap_or_default();
    if target_family == "wasm" {
        panic!(
            "The `native-inchi` feature cannot be enabled for WASM targets \
             (wasm32-unknown-unknown). Remove the feature flag when building for WASM."
        );
    }

    let manifest = std::path::PathBuf::from(std::env::var("CARGO_MANIFEST_DIR").unwrap());
    let vendor = manifest.join("vendor/inchi-src");

    if !vendor.exists() {
        panic!(
            "Vendor directory not found at {vendor:?}.\n\
             Download the IUPAC InChI v1.07.5 source and place it under\n\
             crates/chematic-inchi/vendor/inchi-src/"
        );
    }

    let base = vendor.join("INCHI_BASE/src");
    let lib_src = vendor.join("INCHI_API/libinchi/src");
    let ixa = vendor.join("INCHI_API/libinchi/src/ixa");

    let base_files = [
        "bcf_s.c",
        "ichi_bns.c",
        "ichi_io.c",
        "ichican2.c",
        "ichicano.c",
        "ichicans.c",
        "ichierr.c",
        "ichiisot.c",
        "ichimak2.c",
        "ichimake.c",
        "ichimap1.c",
        "ichimap2.c",
        "ichimap4.c",
        "ichinorm.c",
        "ichiparm.c",
        "ichiprt1.c",
        "ichiprt2.c",
        "ichiprt3.c",
        "ichiqueu.c",
        "ichiread.c",
        "ichiring.c",
        "ichirvr1.c",
        "ichirvr2.c",
        "ichirvr3.c",
        "ichirvr4.c",
        "ichirvr5.c",
        "ichirvr6.c",
        "ichirvr7.c",
        "ichisort.c",
        "ichister.c",
        "ichitaut.c",
        "ikey_base26.c",
        "ikey_dll.c",
        "inchi_gui.c",
        "mol2atom.c",
        "mol_fmt1.c",
        "mol_fmt2.c",
        "mol_fmt3.c",
        "mol_fmt4.c",
        "permutation_util.c",
        "readinch.c",
        "runichi.c",
        "runichi2.c",
        "runichi3.c",
        "runichi4.c",
        "sha2.c",
        "strutil.c",
        "util.c",
    ];

    let lib_files = [
        "ichilnct.c",
        "inchi_dll.c",
        "inchi_dll_a.c",
        "inchi_dll_a2.c",
        "inchi_dll_b.c",
        "inchi_dll_main.c",
    ];

    let ixa_files = [
        "ixa_builder.c",
        "ixa_inchikey_builder.c",
        "ixa_mol.c",
        "ixa_read_inchi.c",
        "ixa_read_mol.c",
        "ixa_status.c",
    ];

    let mut build = cc::Build::new();
    build
        .include(&base)
        .include(&lib_src)
        .include(&ixa)
        .define("TARGET_API_LIB", None)
        .define("COMPILE_ANSI_ONLY", None)
        .flag_if_supported("-fno-strict-aliasing")
        .flag_if_supported("-Wno-all")
        .flag_if_supported("-Wno-everything")
        .warnings(false);

    for f in &base_files {
        build.file(base.join(f));
    }
    for f in &lib_files {
        build.file(lib_src.join(f));
    }
    for f in &ixa_files {
        build.file(ixa.join(f));
    }

    build.compile("inchi_native");

    // Math library required on Linux
    let target_os = std::env::var("CARGO_CFG_TARGET_OS").unwrap_or_default();
    if target_os == "linux" {
        println!("cargo:rustc-link-lib=m");
    }
}
