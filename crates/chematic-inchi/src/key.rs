// Minimal pure-software SHA-256 — replaces the sha2 crate (~56 KB WASM code).
// Saves ~15 KB gzip in the WASM build. Output is identical to sha2::Sha256.
#[rustfmt::skip]
const K: [u32; 64] = [
    0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
    0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
    0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
    0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7, 0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
    0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
    0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
    0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
    0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2,
];

fn sha256(data: &[u8]) -> [u8; 32] {
    let mut h: [u32; 8] = [
        0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a, 0x510e527f, 0x9b05688c, 0x1f83d9ab,
        0x5be0cd19,
    ];
    let bit_len = (data.len() as u64).wrapping_mul(8);
    let mut padded = data.to_vec();
    padded.push(0x80);
    while padded.len() % 64 != 56 {
        padded.push(0);
    }
    padded.extend_from_slice(&bit_len.to_be_bytes());
    for chunk in padded.as_chunks::<64>().0.iter() {
        let mut w = [0u32; 64];
        for i in 0..16 {
            w[i] = u32::from_be_bytes([
                chunk[4 * i],
                chunk[4 * i + 1],
                chunk[4 * i + 2],
                chunk[4 * i + 3],
            ]);
        }
        for i in 16..64 {
            let s0 = w[i - 15].rotate_right(7) ^ w[i - 15].rotate_right(18) ^ (w[i - 15] >> 3);
            let s1 = w[i - 2].rotate_right(17) ^ w[i - 2].rotate_right(19) ^ (w[i - 2] >> 10);
            w[i] = w[i - 16]
                .wrapping_add(s0)
                .wrapping_add(w[i - 7])
                .wrapping_add(s1);
        }
        let [mut a, mut b, mut c, mut d, mut e, mut f, mut g, mut hh] = h;
        for i in 0..64 {
            let s1 = e.rotate_right(6) ^ e.rotate_right(11) ^ e.rotate_right(25);
            let ch = (e & f) ^ (!e & g);
            let t1 = hh
                .wrapping_add(s1)
                .wrapping_add(ch)
                .wrapping_add(K[i])
                .wrapping_add(w[i]);
            let s0 = a.rotate_right(2) ^ a.rotate_right(13) ^ a.rotate_right(22);
            let maj = (a & b) ^ (a & c) ^ (b & c);
            let t2 = s0.wrapping_add(maj);
            (a, b, c, d, e, f, g, hh) = (t1.wrapping_add(t2), a, b, c, d.wrapping_add(t1), e, f, g);
        }
        h[0] = h[0].wrapping_add(a);
        h[1] = h[1].wrapping_add(b);
        h[2] = h[2].wrapping_add(c);
        h[3] = h[3].wrapping_add(d);
        h[4] = h[4].wrapping_add(e);
        h[5] = h[5].wrapping_add(f);
        h[6] = h[6].wrapping_add(g);
        h[7] = h[7].wrapping_add(hh);
    }
    let mut out = [0u8; 32];
    for (i, &hi) in h.iter().enumerate() {
        out[4 * i..4 * i + 4].copy_from_slice(&hi.to_be_bytes());
    }
    out
}

/// Generate InChIKey (27-character identifier) from an InChI string.
///
/// Format: `XXXXXXXXXXXXXX-XXXXXXXXXX-N` where N is the version flag.
pub fn inchi_key(inchi_str: &str) -> String {
    // Remove "InChI=1S/" prefix
    let inchi_content = if let Some(pos) = inchi_str.find("/") {
        &inchi_str[pos..]
    } else {
        inchi_str
    };

    // Split into parts: /c.../h... is connectivity+hydrogen, rest is charge/isotope
    let parts: Vec<&str> = inchi_content.split('/').collect();

    // First block: hash connectivity and hydrogen layers (/c and /h)
    let connectivity_input = if parts.len() >= 3 {
        // /c...../h..... format
        format!("{}/{}", parts[1], parts[2])
    } else if parts.len() >= 2 {
        // Just /c.....
        parts[1].to_string()
    } else {
        inchi_content.to_string()
    };

    // Second block: hash remaining layers (charge, isotope, stereo)
    let remaining_input = if parts.len() > 3 {
        parts[3..].join("/")
    } else {
        String::new()
    };

    // Compute SHA-256 hashes
    let hash1_bytes = sha256_hash(&connectivity_input);
    let hash2_bytes = sha256_hash(&remaining_input);

    // Convert to base-26 (A-Z) representation
    let block1 = bytes_to_base26(&hash1_bytes[..12]); // 12 bytes → ~14 chars
    let block2 = bytes_to_base26(&hash2_bytes[..9]); // 9 bytes → ~10 chars

    // Ensure block1 is exactly 14 chars and block2 is exactly 9 chars
    let block1_padded = format!("{:<14}", block1);
    let block2_padded = format!("{:<9}", block2);

    // Take only first 14 and 10 chars respectively, then truncate to exact size
    let block1_final = &block1_padded[..14];
    let block2_final = &block2_padded[..10];

    format!("{}-{}-N", block1_final, block2_final)
}

fn sha256_hash(input: &str) -> Vec<u8> {
    sha256(input.as_bytes()).to_vec()
}

fn bytes_to_base26(bytes: &[u8]) -> String {
    let mut result = String::new();
    for &byte in bytes {
        // Convert byte to 2 base-26 digits
        let d1 = (byte / 26) % 26;
        let d2 = byte % 26;
        result.push((b'A' + d1) as char);
        result.push((b'A' + d2) as char);
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sha256_abc() {
        // RFC 4634 test vector: SHA-256("abc")
        let got = sha256(b"abc");
        let expected: [u8; 32] = [
            0xba, 0x78, 0x16, 0xbf, 0x8f, 0x01, 0xcf, 0xea, 0x41, 0x41, 0x40, 0xde, 0x5d, 0xae,
            0x22, 0x23, 0xb0, 0x03, 0x61, 0xa3, 0x96, 0x17, 0x7a, 0x9c, 0xb4, 0x10, 0xff, 0x61,
            0xf2, 0x00, 0x15, 0xad,
        ];
        assert_eq!(got, expected);
    }

    #[test]
    fn sha256_empty() {
        // SHA-256("")
        let got = sha256(b"");
        let expected: [u8; 32] = [
            0xe3, 0xb0, 0xc4, 0x42, 0x98, 0xfc, 0x1c, 0x14, 0x9a, 0xfb, 0xf4, 0xc8, 0x99, 0x6f,
            0xb9, 0x24, 0x27, 0xae, 0x41, 0xe4, 0x64, 0x9b, 0x93, 0x4c, 0xa4, 0x95, 0x99, 0x1b,
            0x78, 0x52, 0xb8, 0x55,
        ];
        assert_eq!(got, expected);
    }

    #[test]
    fn test_inchi_key_format() {
        let inchi = "InChI=1S/C6H6/c1-2-3-4-5-6-1/h1-6H";
        let key = inchi_key(inchi);
        assert_eq!(key.len(), 27);
        assert_eq!(&key[14..15], "-");
        assert_eq!(&key[25..26], "-");
        assert_eq!(&key[26..27], "N");
    }

    #[test]
    fn test_inchi_key_deterministic() {
        let inchi = "InChI=1S/C6H6/c1-2-3-4-5-6-1/h1-6H";
        let key1 = inchi_key(inchi);
        let key2 = inchi_key(inchi);
        assert_eq!(key1, key2);
    }

    #[test]
    fn test_inchi_key_different_for_different_inchi() {
        let inchi1 = "InChI=1S/CH4/h1H4";
        let inchi2 = "InChI=1S/C2H6/c1-2/h1-2H3";
        let key1 = inchi_key(inchi1);
        let key2 = inchi_key(inchi2);
        assert_ne!(key1, key2);
    }
}
