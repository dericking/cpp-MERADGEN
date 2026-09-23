#ifdef MERADGEN_PARITY_TRACE

#include "meradgen_parity_trace.hpp"

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>

static FILE* ptrc_fp() {
  static FILE* fp = nullptr;
  static bool once = false;
  if (!once) {
    once = true;
    const char* p = std::getenv("MERADGEN_PARITY_TRACE");
    if (p && p[0])
      fp = std::fopen(p, "w");
  }
  return fp;
}

static std::uint64_t ptrc_bits(double x) {
  std::uint64_t u = 0;
  std::memcpy(&u, &x, 8);
  return u;
}

void ptrc_event() {
  FILE* f = ptrc_fp();
  if (!f)
    return;
  static int n = 0;
  n++;
  std::fprintf(f, "EVENT %d\n", n);
}

void ptrc_d(const char* tag, double x) {
  FILE* f = ptrc_fp();
  if (!f)
    return;
  std::fprintf(f, "D %s 0x%016llX\n", tag, (unsigned long long)ptrc_bits(x));
}

void ptrc_i(const char* tag, int x) {
  FILE* f = ptrc_fp();
  if (!f)
    return;
  std::fprintf(f, "I %s %d\n", tag, x);
}

void ptrc_row(char kind, int idx, double a, double b, double c, double d) {
  FILE* f = ptrc_fp();
  if (!f)
    return;
  std::fprintf(f, "%c %d 0x%016llX 0x%016llX 0x%016llX 0x%016llX\n", kind, idx,
               (unsigned long long)ptrc_bits(a), (unsigned long long)ptrc_bits(b),
               (unsigned long long)ptrc_bits(c), (unsigned long long)ptrc_bits(d));
}

#endif
