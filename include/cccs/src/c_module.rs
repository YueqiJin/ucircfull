use libc::c_char;
use std::ffi::CString;
use std::os::raw::c_uint;

use crate::scan::find_concensus;

// Define a struct for the return value
#[repr(C)]
pub struct CcsResult {
    pub segment_ptr: *mut c_char,
    pub ccs_seq_ptr: *mut c_char,
}

// Define the Rust function with C-compatible types
#[no_mangle]
pub extern "C" fn ccs(seq: *const u8, seq_len: c_uint) -> *mut CcsResult {
    // Convert the input sequence from C to Rust types
    let seq_slice: &[u8];
    unsafe { 
        seq_slice = std::slice::from_raw_parts(seq, seq_len as usize);
    }
    //println!("seq_slice: {:?}", seq_slice);

    // Call the original Rust function
    let result = match find_concensus(seq_slice) {
        Err(_) => CcsResult {
            segment_ptr: std::ptr::null_mut(),
            ccs_seq_ptr: std::ptr::null_mut(),
        },
        Ok(x) => {
            let ccs_seq = x.1;
            let mut segments: Vec<String> = Vec::default();
            for (s, e) in x.0 {
                segments.push(String::from(vec![s.to_string(), e.to_string()].join("-")));
            }
            let segment_str = segments.join(";");
            let segment_cstr = CString::new(segment_str).expect("Failed to convert segment to CString").into_raw();
            let ccs_seq_cstr = CString::new(ccs_seq).expect("Failed to convert ccs_seq to CString").into_raw();
            CcsResult {
                segment_ptr: segment_cstr as *mut c_char,
                ccs_seq_ptr: ccs_seq_cstr as *mut c_char,
            }
        }
    };

    // Allocate memory for the result struct and return the pointer
    Box::into_raw(Box::new(result))
}

// Define a function to free the memory returned by the Rust function
#[no_mangle]
pub extern "C" fn free_ccs_result(result: *mut CcsResult) {
    if !result.is_null() {
        unsafe {
            if !(*result).segment_ptr.is_null() {
                let _ = CString::from_raw((*result).segment_ptr);
            }
            if !(*result).ccs_seq_ptr.is_null() {
                let _ = CString::from_raw((*result).ccs_seq_ptr);
            }
            let _ = Box::from_raw(result);
        }
    }
}