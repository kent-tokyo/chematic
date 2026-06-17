use std::os::raw::{c_char, c_int};

pub const MAXVAL: usize = 20;
pub const ATOM_EL_LEN: usize = 6;
pub const NUM_H_ISOTOPES: usize = 3;

/// Matches `inchi_Atom` from inchi_api.h exactly (field order chosen to avoid padding).
///
/// Size: 3×f64 + 20×i16 + 20×i8 + 20×i8 + 6×i8 + i16 + 4×i8 + i16 + i8 + i8 = 120 bytes.
#[repr(C)]
#[derive(Clone)]
pub struct InchiAtom {
    pub x: f64,
    pub y: f64,
    pub z: f64,
    pub neighbor: [i16; MAXVAL],
    pub bond_type: [i8; MAXVAL],
    pub bond_stereo: [i8; MAXVAL],
    pub elname: [c_char; ATOM_EL_LEN],
    pub num_bonds: i16,
    pub num_iso_h: [i8; NUM_H_ISOTOPES + 1],
    pub isotopic_mass: i16,
    pub radical: i8,
    pub charge: i8,
}

impl Default for InchiAtom {
    fn default() -> Self {
        // SAFETY: all-zero is valid for this Plain-Old-Data struct.
        unsafe { std::mem::zeroed() }
    }
}

/// Matches `inchi_Stereo0D` from inchi_api.h.
#[repr(C)]
#[derive(Clone, Default)]
pub struct InchiStereo0D {
    pub neighbor: [i16; 4],
    pub central_atom: i16,
    pub stereo_type: i8,
    pub parity: i8,
}

/// Matches `inchi_Input` from inchi_api.h.
#[repr(C)]
pub struct InchiInput {
    pub atom: *mut InchiAtom,
    pub stereo0d: *mut InchiStereo0D,
    pub options: *const c_char,
    pub num_atoms: i16,
    pub num_stereo0d: i16,
}

/// Matches `inchi_Output` from inchi_api.h.
#[repr(C)]
pub struct InchiOutput {
    pub sz_inchi: *mut c_char,
    pub sz_aux_info: *mut c_char,
    pub sz_message: *mut c_char,
    pub sz_log: *mut c_char,
}

unsafe extern "C" {
    /// Generate standard InChI from an `inchi_Input`.
    /// Returns 0=OK, 1=warning, 2=error, 3=fatal.
    pub fn GetStdINCHI(inp: *mut InchiInput, out: *mut InchiOutput) -> c_int;

    /// Free memory allocated by `GetStdINCHI`.
    pub fn FreeStdINCHI(out: *mut InchiOutput);

    /// Compute standard InChIKey from a standard InChI string.
    /// `key` must point to a buffer of at least 28 bytes.
    /// Returns 0 on success.
    pub fn GetStdINCHIKeyFromStdINCHI(inchi: *const c_char, key: *mut c_char) -> c_int;
}

/// Compile-time size check for InchiAtom (must be 120 bytes on 64-bit).
const _: () = {
    assert!(std::mem::size_of::<InchiAtom>() == 120);
};
