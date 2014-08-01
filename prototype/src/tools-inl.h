#ifndef TOOLS_INL_H
#define TOOLS_INL_H

#include <sys/time.h>
#include <sys/times.h>
#include <sys/resource.h>

// Classic version
inline int calculate_local_span(int dim_proc_, int dim_procs_, int global_span_) {
  int local_span = global_span_ / dim_procs_;
  int remainder = global_span_ - local_span * dim_procs_;
  if (dim_proc_ < remainder) local_span += 1;
  return local_span;
}

// Classic version
inline int calculate_local_offset(int dim_proc_, int dim_procs_, int global_span_) {
  int offset = 0;
  for (int proc = 0; proc < dim_proc_; ++proc) {
    offset += calculate_local_span(proc, dim_procs_, global_span_);
  }
  return offset;
}

// push/pull decompositon version
inline int calculate_local_span(int dim_proc_,
                                int dim_procs_,
                                int global_span_,
                                bool prograde_) {
  int span = calculate_local_span(dim_proc_, dim_procs_, global_span_);
  // As prograde is the initial case, prograde is equivalent to classic
  if (!prograde_) {
    // On a retrograde step, we seek to pull back towards the origin.
    // This means that the origin will have one less row/col, and the 
    // proc at the furthest extent will have one more
    if (dim_proc_ == 0) {
      span -= 1;
    } else if (dim_proc_ == dim_procs_ - 1) {
      span += 1;
    }
  }
  return span;
}

// push/pull decompositon version
inline int calculate_local_offset(int dim_proc_,
                                  int dim_procs_,
                                  int global_span_,
                                  bool prograde_) {
  int offset = 0;
  for (int proc = 0; proc < dim_proc_; ++proc) {
    offset += calculate_local_span(proc, dim_procs_, global_span_, prograde_);
  }
  return offset;
}

// Utility timing function.
inline void timers(double& wall_, double& cpu_)
{
  struct rusage r;
  struct timeval t;
  gettimeofday( &t, (struct timezone *)0 );
  getrusage (RUSAGE_SELF, &r);
  wall_ = t.tv_sec + t.tv_usec*1.0e-6;
  cpu_ = r.ru_utime.tv_sec + r.ru_utime.tv_usec*1.0e-6; 
}

#endif
