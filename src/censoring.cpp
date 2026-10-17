#include "censoring.h"
#include <algorithm>
#include <numeric>
#include <cmath>
#include "interrupt.h"
namespace competingrisk {
using namespace std;
CensoringDistribution::CensoringDistribution(const std::vector<double>& time,
                                           const std::vector<int>& censored)
  : times(time), cen_ind(censored) { computeKaplanMeierEstimate(); }
void CensoringDistribution::computeKaplanMeierEstimate()
  {
    std::vector<size_t> sorted_indices = sort_index(times);
    sorted_times.resize(times.size());
    for (size_t i = 0; i < sorted_indices.size(); ++i)
    {
      sorted_times[i] = times[sorted_indices[i]];
    }

    std::vector<double> unique_times = sorted_times;
    auto last = unique(unique_times.begin(), unique_times.end());
    unique_times.erase(last, unique_times.end());

    std::vector<int> sorted_events(times.size());
    for (size_t i = 0; i < sorted_indices.size(); ++i)
    {
      sorted_events[i] = cen_ind[sorted_indices[i]];
    }

    size_t m = unique_times.size();
    std::vector<double> at_risk(m, 0);
    std::vector<double> num_events(m, 0);

    size_t index = 0;
    for (size_t i = 0; i < m; ++i)
    {
      interrupt_at(i);
      while (index < times.size() && sorted_times[index] <= unique_times[i])
      {
        at_risk[i] += 1;
        if (sorted_times[index] == unique_times[i] && sorted_events[index] == 1)
        {
          num_events[i] += 1;
        }
        index++;
      }
    }

    for (size_t i = m; i > 1; --i)
    {
      at_risk[i-2] += at_risk[i-1];
    }

    // Keep lookup times aligned with the unique-time survival estimates.
    sorted_times = unique_times;
    survival_prob.resize(m, 1.0);
    double survival = 1.0;
    for (size_t j = 0; j < m; ++j)
    {
      interrupt_at(j);
      if (num_events[j] > 0)
      {
        survival *= (1.0 - num_events[j] / at_risk[j]);
      }
      survival_prob[j] = survival;
    }
  }

double CensoringDistribution::getSurvivalProb(double time) const
  {
    if (time < sorted_times[0])
    {
      return 1.0;
    }

    auto it = upper_bound(sorted_times.begin(), sorted_times.end(), time);
    if (it == sorted_times.begin())
    {
      return survival_prob[0];
    }
    else
    {
      return survival_prob[it - sorted_times.begin() - 1];
    }
  }

std::vector<size_t> CensoringDistribution::sort_index(const std::vector<double> &v)
  {
    std::vector<size_t> idx(v.size());
    iota(idx.begin(), idx.end(), 0);
    sort(idx.begin(), idx.end(), [&v](size_t i1, size_t i2)
    { return v[i1] < v[i2]; });
    return idx;
  }
}

// [[Rcpp::export]]
SEXP censoring_create(Rcpp::NumericVector time, Rcpp::IntegerVector censored) {
  if(time.size()==0 || time.size()!=censored.size()) Rcpp::stop("Invalid KM input lengths.");
  for(int i=0;i<time.size();++i)
    if(!std::isfinite(time[i]) || (censored[i]!=0 && censored[i]!=1)) Rcpp::stop("Invalid KM input.");
  return Rcpp::XPtr<competingrisk::CensoringDistribution>(
    new competingrisk::CensoringDistribution(Rcpp::as<std::vector<double>>(time),
                                           Rcpp::as<std::vector<int>>(censored)),true);
}

// Internal reference-test query; fitting shares the pointer created above.
// [[Rcpp::export]]
Rcpp::NumericVector censoring_survival(Rcpp::NumericVector time,
                                     Rcpp::IntegerVector censored,
                                     Rcpp::NumericVector query) {
  Rcpp::XPtr<competingrisk::CensoringDistribution> km(censoring_create(time,censored));
  Rcpp::NumericVector out(query.size());
  for(int i=0;i<query.size();++i) {
    if(!std::isfinite(query[i])) Rcpp::stop("KM query must be finite.");
    out[i]=km->getSurvivalProb(query[i]);
  }
  return out;
}
