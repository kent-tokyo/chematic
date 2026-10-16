use super::ChemDrawError as Error;
use super::binary::Reader;

fn number(data: &[u8], width: usize, offset: usize) -> Result<u32, Error> {
    if data.len() != width {
        return Err(Error::Invalid {
            offset,
            detail: "property width",
        });
    }
    Ok(data
        .iter()
        .enumerate()
        .fold(0u32, |v, (i, b)| v | (u32::from(*b) << (8 * i))))
}
pub(super) fn attribute(
    r: &mut Reader<'_>,
    owner: u16,
    tag: u16,
    offset: usize,
    start: usize,
    n: usize,
) -> Result<String, Error> {
    let data = &r.input[start..start + n];
    let u16value = || number(data, 2, offset);
    let attr = match (owner, tag) {
        (0x8004, 0x0400) if u16value()? == 1 => String::new(),
        (0x8004, 0x0402) => {
            let value = u16value()?;
            if value == 0 || value > 118 {
                return Err(Error::Invalid {
                    offset,
                    detail: "element atomic number",
                });
            }
            format!(" Element=\"{value}\"")
        }
        (0x8004, 0x0200) => {
            if n != 8 {
                return Err(Error::Invalid {
                    offset,
                    detail: "2D point width",
                });
            }
            let y = i32::from_le_bytes(data[..4].try_into().unwrap()) as f64 / 65536.;
            let x = i32::from_le_bytes(data[4..].try_into().unwrap()) as f64 / 65536.;
            format!(" p=\"{x} {y}\"")
        }
        (0x8004, 0x0420) => {
            let value = u16value()?;
            if value > i16::MAX as u32 {
                return Err(Error::Invalid {
                    offset,
                    detail: "negative isotope mass",
                });
            }
            format!(" Isotope=\"{value}\"")
        }
        (0x8004, 0x0421) => format!(" Charge=\"{}\"", number(data, 1, offset)? as u8 as i8),
        (0x8004, 0x042b) => {
            let value = u16value()?;
            if value > 255 {
                return Err(Error::Invalid {
                    offset,
                    detail: "hydrogen count range",
                });
            }
            format!(" NumHydrogens=\"{value}\"")
        }
        (0x8005, 0x0604 | 0x0605) => {
            let value = number(data, 4, offset)?;
            if value == 0 {
                return Err(Error::Invalid {
                    offset,
                    detail: "zero atom reference",
                });
            }
            format!(" {}=\"{value}\"", if tag == 0x0604 { "B" } else { "E" })
        }
        (0x8005, 0x0600) => {
            let value = match u16value()? {
                1 => "1",
                2 => "2",
                4 => "3",
                128 => "1.5",
                _ => return Err(Error::Unsupported { offset, tag }),
            };
            format!(" Order=\"{value}\"")
        }
        (0x8005, 0x0601 | 0x0602) => {
            if u16value()? != 0 {
                return Err(Error::Unsupported { offset, tag });
            }
            String::new()
        }
        // Known presentation-only properties. Unknown/custom tags and all
        // unhandled chemistry (queries, radicals, enhanced stereo, etc.) fail.
        (_, 0x0001..=0x000a | 0x000f..=0x0011 | 0x0204 | 0x0300..=0x0302 | 0x0800..=0x081f)
        | (0x8004, 0x0401 | 0x0439..=0x043c)
        | (0x8005, 0x0603 | 0x060c..=0x060f)
        | (0x8006, 0x0700..=0x0704) => {
            r.omit(tag, offset)?;
            String::new()
        }
        _ => return Err(Error::Unsupported { offset, tag }),
    };
    Ok(attr)
}
