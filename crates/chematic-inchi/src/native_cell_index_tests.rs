//! Boundary regression for the vendored InChI cell traversal.

use std::process::Command;
use std::time::{Duration, Instant};

#[repr(C)]
struct Partition {
    rank: *mut u16,
    at_number: *mut u16,
}

#[repr(C)]
struct Cell {
    first: i32,
    next: i32,
    prev: i32,
}

unsafe extern "C" {
    fn CellGetMinNode(p: *mut Partition, w: *mut Cell, v: u16, cd: *mut std::ffi::c_void) -> u16;
}

#[test]
fn cell_index_can_advance_past_u16_max() {
    const CHILD: &str = "CHEMATIC_CELL_INDEX_TEST_CHILD";
    if std::env::var_os(CHILD).is_some() {
        // Exercise the loop boundary with synthetic storage. This does not claim
        // that a molecule with this many atoms is accepted by InChI.
        let mut atoms = vec![0u16; usize::from(u16::MAX) + 1];
        let mut ranks = [0u16];
        let mut partition = Partition {
            rank: ranks.as_mut_ptr(),
            at_number: atoms.as_mut_ptr(),
        };
        let mut cell = Cell {
            first: i32::from(u16::MAX),
            next: i32::from(u16::MAX) + 1,
            prev: 0,
        };
        // Both arrays cover every index accessed by this synthetic cell.
        let result = unsafe { CellGetMinNode(&mut partition, &mut cell, 0, std::ptr::null_mut()) };
        assert_eq!(result, 1);
        return;
    }

    // A narrow loop counter wraps and never terminates; isolate it so a
    // regression fails instead of hanging the whole test suite.
    let mut child = Command::new(std::env::current_exe().unwrap())
        .args([
            "--exact",
            "native_cell_index_tests::cell_index_can_advance_past_u16_max",
        ])
        .env(CHILD, "1")
        .spawn()
        .unwrap();
    let deadline = Instant::now() + Duration::from_secs(10);
    loop {
        if let Some(status) = child.try_wait().unwrap() {
            assert!(status.success(), "boundary test failed: {status}");
            break;
        }
        if Instant::now() >= deadline {
            child.kill().unwrap();
            child.wait().unwrap();
            panic!("cell traversal did not terminate");
        }
        std::thread::sleep(Duration::from_millis(20));
    }
}
