#ifndef MERADGEN_PARITY_TRACE_HPP
#define MERADGEN_PARITY_TRACE_HPP

#ifdef MERADGEN_PARITY_TRACE
void ptrc_event();
void ptrc_d(const char* tag, double x);
void ptrc_i(const char* tag, int x);
void ptrc_row(char kind, int idx, double a, double b, double c, double d);
#else
inline void ptrc_event() {}
inline void ptrc_d(const char*, double) {}
inline void ptrc_i(const char*, int) {}
inline void ptrc_row(char, int, double, double, double, double) {}
#endif

#endif
