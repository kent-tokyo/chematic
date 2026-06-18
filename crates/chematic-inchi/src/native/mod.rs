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

use std::ffi::{CStr, CString};
use std::os::raw::c_char;
use chematic_core::Molecule;


use convert::ConvertError;


use ffi::{FreeStdINCHI, GetStdINCHI, GetStdINCHIKeyFromStdINCHI, InchiInput, InchiOutput};

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
    let (mut atoms, mut stereo) = convert::mol_to_inchi_atoms(mol).map_err(|e| match e {
        ConvertError::KekulizationFailed(msg) => InchiError::KekulizationFailed(msg),
    })?;

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
        atom: if atoms.is_empty() {
            std::ptr::null_mut()
        } else {
            atoms.as_mut_ptr()
        },
        stereo0d: if stereo.is_empty() {
            std::ptr::null_mut()
        } else {
            stereo.as_mut_ptr()
        },
        options: std::ptr::null(),
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
            unsafe { CStr::from_ptr(output.sz_message).to_string_lossy().into_owned() }
        };
        unsafe { FreeStdINCHI(&mut output) };
        return Err(InchiError::LibError(msg));
    }

    if output.sz_inchi.is_null() {
        unsafe { FreeStdINCHI(&mut output) };
        return Err(InchiError::NullOutput);
    }

    let inchi = unsafe { CStr::from_ptr(output.sz_inchi).to_string_lossy().into_owned() };
    unsafe { FreeStdINCHI(&mut output) };
    Ok(inchi)
}

/// Compute the standard 27-character InChIKey from a standard InChI string.
pub fn standard_inchi_key(inchi_str: &str) -> Result<String, InchiError> {
    let c_inchi = CString::new(inchi_str)
        .map_err(|e| InchiError::InvalidInput(e.to_string()))?;

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
