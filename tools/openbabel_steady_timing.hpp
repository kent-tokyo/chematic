#pragma once
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <stdexcept>
#include <string>

namespace chematic_benchmark {
struct Timing {
  std::uint64_t units, repeats, elapsed_ns, warmup_repeats, warmup_elapsed_ns;
  std::uint64_t minimum_warmup_ns, minimum_timed_ns;
};
inline std::uint64_t number(const char* name) {
  const char* value = std::getenv(name);
  const auto n = value ? std::stoull(value) : 0;
  if (n > 10000000000ULL) throw std::runtime_error("benchmark duration exceeds ten seconds");
  return n;
}
inline bool warmed_dump() {
  const char* value = std::getenv("CHEMATIC_BENCH_WARM_DUMP");
  return value && std::string(value) == "1";
}
inline std::uint64_t dump_repeats() {
  const char* value = std::getenv("CHEMATIC_BENCH_MIN_REPEATS");
  return value ? std::stoull(value) : 3;
}
template<class Work> Timing measure(std::uint64_t minimum_repeats, Work work) {
  if (!minimum_repeats) throw std::runtime_error("positive minimum repetitions required");
  const auto warmup_target = number("CHEMATIC_BENCH_WARMUP_NS");
  const auto timed_target = number("CHEMATIC_BENCH_TIMED_NS");
  using Clock = std::chrono::steady_clock;
  auto ns_since = [](auto start) -> std::uint64_t {
    return std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now() - start).count();
  };
  const auto warmup_started = Clock::now();
  std::uint64_t warmup_repeats = 0, warmup_units = 0;
  while (ns_since(warmup_started) < warmup_target) {
    warmup_units += work();
    ++warmup_repeats;
  }
  // Keep the warm workload observable without charging its units to timing.
  volatile auto observed_warmup_units = warmup_units;
  (void)observed_warmup_units;
  const auto warmup_elapsed = ns_since(warmup_started);
  const auto started = Clock::now();
  std::uint64_t units = 0, repeats = 0;
  do { units += work(); ++repeats; }
  while (repeats < minimum_repeats || (timed_target && ns_since(started) < timed_target));
  return {units, repeats, ns_since(started), warmup_repeats, warmup_elapsed, warmup_target, timed_target};
}
}  // namespace chematic_benchmark
