use sha2::{Digest, Sha256};

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
    let mut hasher = Sha256::new();
    hasher.update(input.as_bytes());
    hasher.finalize().to_vec()
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
