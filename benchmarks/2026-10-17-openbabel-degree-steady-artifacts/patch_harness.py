from pathlib import Path
import re,shutil
O=Path(__file__).resolve().parent

def apply(root):
 root=Path(root);support=root/'crates/chematic-mol/examples/support';support.mkdir(exist_ok=True)
 shutil.copy2(O/'steady_timing.rs',support/'steady_timing.rs')
 shutil.copy2(O/'openbabel_steady_timing.hpp',root/'tools/openbabel_steady_timing.hpp')
 for kind in ['sdf','v3000','v3000_3d']:
  p=root/'crates/chematic-mol/examples'/('openbabel_'+kind+'_throughput.rs');s=p.read_text()
  s=s.replace(', time::Instant','').replace('    time::Instant,\n','')
  i=s.index('\nfn ');s=s[:i]+'\n#[path = "support/steady_timing.rs"]\nmod steady_timing;\n'+s[i:]
  if kind=='sdf':
   mode='mode.as_str()';rd='read(black_box(&input), expected)?';wr='write(&prepared, &mut buffer)';rt='write(&read(black_box(&input), expected)?, &mut buffer)'
  else:
   mode='args[2].as_str()';rd='read(black_box(&blocks))?';wr='write(&prepared)';rt='write(&read(black_box(&blocks))?)'
  needle='        let mut stdout = std::io::stdout().lock();'
  assert s.count(needle)==1
  warm=f'''        if steady_timing::warmed_dump() {{
            steady_timing::measure(steady_timing::dump_repeats()?, || Ok({wr}))?;
        }}
'''
  s=s.replace(needle,warm+needle)
  a=s.index('    let mut units = 0usize;');b=s.index('    println!(',a)
  s=s[:a]+f'''    let timing = steady_timing::measure(repeats, || {{
        Ok(match {mode} {{
            "parse" => {{
                let records = black_box({rd});
                records.iter().map(|r| r.mol.atom_count()).sum::<usize>()
            }}
            "write" => {wr},
            "roundtrip" => {rt},
            _ => return Err("unknown operation".into()),
        }})
    }})?;
'''+s[b:]
  s=s.replace('expected * repeats','expected * timing.repeats').replace('expected*repeats','expected*timing.repeats')
  s=s.replace('"repeats":repeats','"repeats":timing.repeats').replace('"output_units":units','"output_units":timing.units').replace('"elapsed_ns":started.elapsed().as_nanos()', '''"elapsed_ns":timing.elapsed_ns,
        "warmup_repeats":timing.warmup_repeats,"warmup_elapsed_ns":timing.warmup_elapsed_ns,
        "minimum_warmup_ns":timing.minimum_warmup_ns,"minimum_timed_ns":timing.minimum_timed_ns''')
  p.write_text(s)
 for kind in ['sdf','v3000']:
  p=root/'tools'/('openbabel_'+kind+'_throughput.cpp');s=p.read_text();s='#include "openbabel_steady_timing.hpp"\n'+s
  wr='write_records(prepared, writer)' if kind=='sdf' else 'write(prepared, writer)'
  rd='read_records(input, expected)' if kind=='sdf' else 'read(blocks)'
  write_read='write_records(read_records(input, expected), writer)' if kind=='sdf' else 'write(read(blocks), writer)'
  needle='    if (mode == "dump") {';assert s.count(needle)==1
  s=s.replace(needle,needle+f'''\n      if (chematic_benchmark::warmed_dump()) {{
        chematic_benchmark::measure(chematic_benchmark::dump_repeats(), [&]() {{ return {wr}; }});
      }}''')
  a=s.index('    std::size_t units = 0;');b=s.index('    std::cout << "{',a)
  s=s[:a]+f'''    const auto timing = chematic_benchmark::measure(repeats, [&]() -> std::uint64_t {{
      if (mode == "parse") {{
        const auto molecules = {rd};
        std::uint64_t units = 0;
        for (const auto& molecule : molecules) units += molecule->NumAtoms();
        return units;
      }}
      if (mode == "write") return {wr};
      if (mode == "roundtrip") return {write_read};
      throw std::runtime_error("unknown mode");
    }});
'''+s[b:]
  s=s.replace('expected * repeats','expected * timing.repeats').replace('<< repeats','<< timing.repeats').replace('<< units','<< timing.units').replace('<< elapsed','<< timing.elapsed_ns')
  s=s.replace('<< "}\\n";', '''<< ",\\"warmup_repeats\\":" << timing.warmup_repeats
              << ",\\"warmup_elapsed_ns\\":" << timing.warmup_elapsed_ns
              << ",\\"minimum_warmup_ns\\":" << timing.minimum_warmup_ns
              << ",\\"minimum_timed_ns\\":" << timing.minimum_timed_ns << "}\\n";''')
  p.write_text(s)

if __name__=='__main__':
 import sys
 for root in sys.argv[1:]:apply(root)
