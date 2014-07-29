#ifndef TOOLS_INL_H
#define TOOLS_INL_H
inline int calculate_local_span(int dim_proc_, int dim_procs_, int global_span_) {
  int local_span = global_span_ / dim_procs_;
  int remainder = global_span_ - local_span * dim_procs_;
  if (dim_proc_ < remainder) local_span += 1;
  return local_span;
}

inline int calculate_local_offset(int dim_proc_, int dim_procs_, int global_span_) {
  int offset = 0;
  for (int proc = 0; proc < dim_proc_; ++proc) {
    offset += calculate_local_span(proc, dim_procs_, global_span_);
  }
  return offset;
}



#endif
