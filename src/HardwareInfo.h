#ifndef FASTQNSN_HARDWARE_INFO_H
#define FASTQNSN_HARDWARE_INFO_H

#include <fstream>
#include <iostream>
#include <string>
#include <thread>
#include <vector>

#if defined(__linux__)
#include <unistd.h>
#elif defined(__APPLE__)
#include <sys/sysctl.h>
#elif defined(_WIN32)
#include <limits>
#include <windows.h>
#endif

#if defined(__x86_64__) || defined(_M_X64) || defined(__i386__) ||             \
    defined(_M_IX86)
#include <cpuid.h>
#endif

namespace fastqnsn_dynamic {

enum class SIMDLevel { None, SSE4_2, AVX2, AVX512, Neon };

struct HardwareInfo {
  size_t l1_cache_size = 32768;   // Default 32KB
  size_t l2_cache_size = 262144;  // Default 256KB
  size_t l3_cache_size = 4194304; // Default 4MB
  size_t cache_line_size = 64;
  size_t num_physical_cores = 1;
  size_t num_logical_cores = 1;
  SIMDLevel simd_level = SIMDLevel::None;

  HardwareInfo() { discover(); }

  void discover() {
    num_logical_cores = std::thread::hardware_concurrency();

#if defined(__linux__)
    discover_linux();
#elif defined(__APPLE__)
    discover_macos();
#elif defined(_WIN32)
    discover_windows();
#endif

    discover_simd();
  }

private:
  void discover_linux() {
    // Try sysfs first
    auto read_sysfs = [](const std::string &path) -> size_t {
      std::ifstream f(path);
      if (!f.is_open())
        return 0;
      size_t val;
      f >> val;
      return val;
    };

    // L2 cache size (usually index 2)
    size_t l2 = read_sysfs("/sys/devices/system/cpu/cpu0/cache/index2/size");
    if (l2 > 0)
      l2_cache_size = l2 * 1024; // If in KB

    // Cache line size
    size_t cls = read_sysfs(
        "/sys/devices/system/cpu/cpu0/cache/index0/coherency_line_size");
    if (cls > 0)
      cache_line_size = cls;

// Fallback to sysconf
#ifdef _SC_LEVEL2_CACHE_SIZE
    if (l2 == 0) {
      long l2_sc = sysconf(_SC_LEVEL2_CACHE_SIZE);
      if (l2_sc > 0)
        l2_cache_size = (size_t)l2_sc;
    }
#endif
  }

  void discover_macos() {
#ifdef __APPLE__
    size_t len = sizeof(size_t);
    sysctlbyname("hw.l2cachesize", &l2_cache_size, &len, NULL, 0);
    sysctlbyname("hw.cachelinesize", &cache_line_size, &len, NULL, 0);
    sysctlbyname("hw.physicalcpu", &num_physical_cores, &len, NULL, 0);
#endif
  }

  void discover_windows() {
#ifdef _WIN32
// Simplistic fallback for now, real implementation would use
// GetLogicalProcessorInformation
#endif
  }

  void discover_simd() {
#if defined(__x86_64__) || defined(_M_X64)
    if (__builtin_cpu_supports("avx512f"))
      simd_level = SIMDLevel::AVX512;
    else if (__builtin_cpu_supports("avx2"))
      simd_level = SIMDLevel::AVX2;
    else if (__builtin_cpu_supports("sse4.2"))
      simd_level = SIMDLevel::SSE4_2;
#elif defined(__aarch64__)
    simd_level = SIMDLevel::Neon; // Always on ARMv8
#endif
  }
};

} // namespace fastqnsn_dynamic

#endif // FASTQNSN_HARDWARE_INFO_H
