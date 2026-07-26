//! Standard IUPAC InChI generation via the vendored IUPAC InChI C library (v1.07.5).
//!
//! Enabled only with the `native-inchi` feature. Not available for WASM targets.
//!
//! # Stereo note
//! Tetrahedral stereo (/t, /m, /s) and E/Z double bond stereo (/b) are included.
//! Formula, connectivity, hydrogen, charge and isotope layers are bit-exact with
//! the IUPAC reference.

mod convert;
pub(crate) mod ffi;

use chematic_core::Molecule;
use std::ffi::{CStr, CString};
use std::os::raw::c_char;

use convert::ConvertError;
pub(crate) use convert::{
    has_unrepresentable_multi_h_stereocenter, has_unresolved_specified_tetrahedral_stereo,
    unresolved_specified_tetrahedral_stereo_atoms,
};

use ffi::{FreeStdINCHI, GetStdINCHI, GetStdINCHIKeyFromStdINCHI, InchiInput, InchiOutput};

/// The IUPAC InChI C library's command-line option prefix character is fixed
/// at the *library's own* compile time (`INCHI_OPTION_PREFX` in `mode.h`:
/// `/` on Windows, `-` everywhere else) -- `cfg(windows)` matches it here
/// since the vendored C library and this Rust code are always built for the
/// same target.
#[cfg(windows)]
const SNON_OPTION: &CStr = c"/SNon";
#[cfg(not(windows))]
const SNON_OPTION: &CStr = c"-SNon";

/// Error returned by [`standard_inchi`] and [`standard_inchi_key`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum InchiError {
    /// The C library returned a null InChI pointer.
    NullOutput,
    /// The C library returned an error or fatal code.
    LibError(String),
    /// The input string contained an interior NUL byte.
    InvalidInput(String),
    /// Aromatic bonds in the molecule could not be Kekulized.
    ///
    /// Occurs for some complex fused ring systems. The IUPAC C library
    /// requires explicit single/double/triple bonds; aromatic bonds must be
    /// converted first.
    KekulizationFailed(String),
}

impl std::fmt::Display for InchiError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::NullOutput => write!(f, "InChI library returned null output"),
            Self::LibError(msg) => write!(f, "InChI library error: {msg}"),
            Self::InvalidInput(msg) => write!(f, "invalid InChI input: {msg}"),
            Self::KekulizationFailed(msg) => write!(f, "Kekulization failed: {msg}"),
        }
    }
}

impl std::error::Error for InchiError {}

/// Generate a standard IUPAC InChI string using the vendored InChI C library.
///
/// Layers included: formula, connectivity (/c), hydrogen (/h), charge (/q if non-zero),
/// isotope (/i if present), tetrahedral stereo (/t, /m, /s) and E/Z double bond stereo (/b).
///
/// Returns an `Err` only if the C library signals an error. Warnings are accepted
/// as success (return code 1).
pub fn standard_inchi(mol: &Molecule) -> Result<String, InchiError> {
    generate_standard_inchi(mol, std::ptr::null())
}

/// Generate a standard IUPAC InChI string with stereo perception suppressed
/// at generation time, via the InChI C API's own `SNon` option ("exclude
/// stereo" -- one of the options explicitly documented as compatible with
/// `GetStdINCHI`'s standard-format enforcement, unlike e.g. `/SRel`). This
/// clears the library's internal `REQ_MODE_STEREO` bit before any stereo
/// perception runs, so the returned string has no `/b`, `/t`, `/m`, or `/s`
/// layer at all -- **not** a full InChI string with those layers stripped
/// out afterward by this crate.
///
/// `pub(crate)`: used only by [`crate::dedup::IdentityPolicy::StereoIgnored`].
pub(crate) fn standard_inchi_no_stereo(mol: &Molecule) -> Result<String, InchiError> {
    generate_standard_inchi(mol, SNON_OPTION.as_ptr())
}

/// Shared worker behind [`standard_inchi`] and [`standard_inchi_no_stereo`].
/// `options` is either null (default: full stereo) or a NUL-terminated
/// options string understood by the InChI C library's option parser.
fn generate_standard_inchi(mol: &Molecule, options: *const c_char) -> Result<String, InchiError> {
    let (mut atoms, mut stereo) = convert::mol_to_inchi_atoms(mol).map_err(|e| match e {
        ConvertError::KekulizationFailed(msg) => InchiError::KekulizationFailed(msg),
    })?;

    // Guard: InChI is not defined for molecules with no heavy atoms (e.g. [H][H]).
    // The C library returns NullOutput for num_atoms=0; return a clean error instead.
    if atoms.is_empty() {
        return Err(InchiError::InvalidInput(
            "molecule has no heavy atoms; InChI is not defined for pure-hydrogen molecules"
                .to_string(),
        ));
    }

    // The InChI C library uses i16 (AT_NUM / signed short) for atom counts.
    // Guard against silent wrapping: i16::MAX = 32767, but the library rejects
    // anything above its internal NORMALLY_ALLOWED_INP_MAX_ATOMS (1024 by default)
    // anyway. We check here so the error is a clean Rust Err rather than UB.
    if atoms.len() > 32766 {
        return Err(InchiError::InvalidInput(format!(
            "molecule has {} heavy atoms; InChI C library maximum is 32766",
            atoms.len()
        )));
    }

    let mut input = InchiInput {
        atom: atoms.as_mut_ptr(), // non-empty guaranteed by the guard above
        stereo0d: if stereo.is_empty() {
            std::ptr::null_mut()
        } else {
            stereo.as_mut_ptr()
        },
        options,
        num_atoms: atoms.len() as i16,
        num_stereo0d: stereo.len() as i16,
    };

    let mut output = InchiOutput {
        sz_inchi: std::ptr::null_mut(),
        sz_aux_info: std::ptr::null_mut(),
        sz_message: std::ptr::null_mut(),
        sz_log: std::ptr::null_mut(),
    };

    // SAFETY: `input` and `output` are valid; we free `output` before returning.
    let ret = unsafe { GetStdINCHI(&mut input, &mut output) };

    // ret: 0=OK, 1=warning (both success), 2=error, 3=fatal
    if ret > 1 {
        let msg = if output.sz_message.is_null() {
            String::from("unknown error")
        } else {
            unsafe {
                CStr::from_ptr(output.sz_message)
                    .to_string_lossy()
                    .into_owned()
            }
        };
        unsafe { FreeStdINCHI(&mut output) };
        return Err(InchiError::LibError(msg));
    }

    if output.sz_inchi.is_null() {
        unsafe { FreeStdINCHI(&mut output) };
        return Err(InchiError::NullOutput);
    }

    let inchi = unsafe {
        CStr::from_ptr(output.sz_inchi)
            .to_string_lossy()
            .into_owned()
    };
    unsafe { FreeStdINCHI(&mut output) };
    Ok(inchi)
}

/// Compute the standard 27-character InChIKey from a standard InChI string.
pub fn standard_inchi_key(inchi_str: &str) -> Result<String, InchiError> {
    let c_inchi = CString::new(inchi_str).map_err(|e| InchiError::InvalidInput(e.to_string()))?;

    // 27 characters + null terminator = 28 bytes.
    let mut key_buf = vec![0 as c_char; 28];

    // SAFETY: `c_inchi` is valid; `key_buf` is 28 bytes as required.
    let ret = unsafe { GetStdINCHIKeyFromStdINCHI(c_inchi.as_ptr(), key_buf.as_mut_ptr()) };

    if ret != 0 {
        return Err(InchiError::LibError(format!(
            "GetStdINCHIKeyFromStdINCHI returned error code {ret}"
        )));
    }

    let key = unsafe {
        CStr::from_ptr(key_buf.as_ptr())
            .to_string_lossy()
            .into_owned()
    };
    Ok(key)
}

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;

    fn mol(s: &str) -> Molecule {
        parse(s).unwrap_or_else(|e| panic!("parse {s:?}: {e}"))
    }

    /// `standard_inchi_no_stereo` must produce a string with no `/b`, `/t`,
    /// `/m`, or `/s` layer at all (generation-time suppression), not the
    /// full string with those layers cut out afterward -- there is no
    /// string-splitting/substring-removal step anywhere in this function.
    #[test]
    fn no_stereo_generation_omits_all_stereo_layers() {
        // (E)-but-2-ene: standard_inchi has a /b layer.
        let ez = mol("C/C=C/C");
        let full = standard_inchi(&ez).unwrap();
        assert!(full.contains("/b"), "sanity: full InChI has /b: {full}");
        let no_stereo = standard_inchi_no_stereo(&ez).unwrap();
        for tag in ["/b", "/t", "/m", "/s"] {
            assert!(
                !no_stereo.contains(tag),
                "SNon output must have no {tag} layer: {no_stereo}"
            );
        }

        // L-alanine: standard_inchi has /t (+ /m/s).
        let ala = mol("N[C@@H](C)C(=O)O");
        let full = standard_inchi(&ala).unwrap();
        assert!(full.contains("/t"), "sanity: full InChI has /t: {full}");
        let no_stereo = standard_inchi_no_stereo(&ala).unwrap();
        for tag in ["/b", "/t", "/m", "/s"] {
            assert!(
                !no_stereo.contains(tag),
                "SNon output must have no {tag} layer: {no_stereo}"
            );
        }
    }

    /// An E/Z pair must produce byte-identical output under `SNon`
    /// generation (both entirely lack stereo now), even though their full
    /// `standard_inchi` output differs.
    #[test]
    fn no_stereo_generation_collapses_ez_pair() {
        let e = mol("C/C=C/C");
        let z = mol("C/C=C\\C");
        assert_ne!(standard_inchi(&e).unwrap(), standard_inchi(&z).unwrap());
        assert_eq!(
            standard_inchi_no_stereo(&e).unwrap(),
            standard_inchi_no_stereo(&z).unwrap()
        );
    }

    /// Same, for an enantiomer pair (tetrahedral stereo).
    #[test]
    fn no_stereo_generation_collapses_enantiomer_pair() {
        let l = mol("N[C@@H](C)C(=O)O");
        let d = mol("N[C@H](C)C(=O)O");
        assert_ne!(standard_inchi(&l).unwrap(), standard_inchi(&d).unwrap());
        assert_eq!(
            standard_inchi_no_stereo(&l).unwrap(),
            standard_inchi_no_stereo(&d).unwrap()
        );
    }

    /// Non-stereo layers (formula, connectivity, H, isotope) must be
    /// unaffected by `SNon` -- it suppresses stereo perception only.
    #[test]
    fn no_stereo_generation_preserves_non_stereo_layers() {
        let d4 = mol("[2H]C([2H])([2H])[2H]");
        let full = standard_inchi(&d4).unwrap();
        let no_stereo = standard_inchi_no_stereo(&d4).unwrap();
        assert!(full.contains("/i1D4"), "sanity: {full}");
        assert!(
            no_stereo.contains("/i1D4"),
            "isotope layer must survive SNon: {no_stereo}"
        );
    }
}
