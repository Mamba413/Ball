#ifndef Ball_OMP_H_
#define Ball_OMP_H_
#if defined(_OPENMP)
#include <omp.h>
#else
static inline int omp_get_thread_num(void) { return 0; }
static inline int omp_get_num_threads(void) { return 1; }
static inline int omp_get_num_procs(void) { return 1; }
static inline void omp_set_num_threads(int nthread) {}
static inline void omp_set_dynamic(int flag) {}
#endif
#endif //Ball_OMP_H_
