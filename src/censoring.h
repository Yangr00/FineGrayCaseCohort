#ifndef COMPETINGRISK_CENSORING_H
#define COMPETINGRISK_CENSORING_H
#include <Rcpp.h>
#include <vector>
namespace competingrisk {
// Right-continuous censoring KM, shared by all modules within one fit.
class CensoringDistribution {
  std::vector<double> times;
  std::vector<int> cen_ind;
  std::vector<double> sorted_times, survival_prob;
  void computeKaplanMeierEstimate();
  std::vector<size_t> sort_index(const std::vector<double>& v);
public:
  CensoringDistribution(const std::vector<double>& time, const std::vector<int>& censored);
  double getSurvivalProb(double time) const;
};
}
#endif
