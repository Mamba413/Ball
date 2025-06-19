#ifndef Ball_OMP_H_
#define Ball_OMP_H_
#if defined(_OPENMP)
#include <omp.h>
#else
static inline int omp_get_thread_num() { return 0; }
static inline int omp_get_num_threads() { return 1; }
static inline int omp_get_num_procs() { return 1; }
static inline void omp_set_num_threads(int nthread) {}
static inline void omp_set_dynamic(int flag) {}
#endif
#endif //Ball_OMP_H_