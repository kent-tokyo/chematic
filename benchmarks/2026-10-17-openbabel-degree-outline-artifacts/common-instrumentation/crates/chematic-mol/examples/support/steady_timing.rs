use std::{error::Error, hint::black_box, time::Instant};

pub struct Timing {
    pub units: usize,
    pub repeats: usize,
    pub elapsed_ns: u128,
    pub warmup_repeats: usize,
    pub warmup_elapsed_ns: u128,
    pub minimum_warmup_ns: u64,
    pub minimum_timed_ns: u64,
}

fn number(name: &str) -> Result<u64, Box<dyn Error>> {
    let value = std::env::var(name).unwrap_or_else(|_| "0".into()).parse()?;
    if value > 10_000_000_000 {
        return Err("benchmark duration exceeds ten seconds".into());
    }
    Ok(value)
}

pub fn warmed_dump() -> bool {
    std::env::var("CHEMATIC_BENCH_WARM_DUMP").is_ok_and(|value| value == "1")
}

pub fn dump_repeats() -> Result<usize, Box<dyn Error>> {
    Ok(std::env::var("CHEMATIC_BENCH_MIN_REPEATS")
        .unwrap_or_else(|_| "3".into())
        .parse()?)
}

pub fn measure(
    minimum_repeats: usize,
    mut work: impl FnMut() -> Result<usize, Box<dyn Error>>,
) -> Result<Timing, Box<dyn Error>> {
    if minimum_repeats == 0 {
        return Err("positive minimum repetitions required".into());
    }
    let minimum_warmup_ns = number("CHEMATIC_BENCH_WARMUP_NS")?;
    let minimum_timed_ns = number("CHEMATIC_BENCH_TIMED_NS")?;
    let warmup_started = Instant::now();
    let mut warmup_repeats = 0;
    while warmup_started.elapsed().as_nanos() < u128::from(minimum_warmup_ns) {
        black_box(work()?);
        warmup_repeats += 1;
    }
    let warmup_elapsed_ns = warmup_started.elapsed().as_nanos();
    let (mut units, mut repeats) = (0usize, 0usize);
    let started = Instant::now();
    loop {
        units += work()?;
        repeats += 1;
        if repeats >= minimum_repeats
            && (minimum_timed_ns == 0
                || started.elapsed().as_nanos() >= u128::from(minimum_timed_ns))
        {
            break;
        }
    }
    Ok(Timing {
        units,
        repeats,
        elapsed_ns: started.elapsed().as_nanos(),
        warmup_repeats,
        warmup_elapsed_ns,
        minimum_warmup_ns,
        minimum_timed_ns,
    })
}
