#ifndef COMPETINGRISK_INTERRUPT_H
#define COMPETINGRISK_INTERRUPT_H
#include <Rcpp.h>
namespace competingrisk {
// Poll outside the numerical inner loops. Rcpp interruption unwinds C++ owners.
inline void interrupt_at(std::size_t iteration) {
  if ((iteration & 63U) == 0U) Rcpp::checkUserInterrupt();
}
}
#endif
