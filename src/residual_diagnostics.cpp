#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <numeric>
#include <cmath>
#include "interrupt.h"
#include "residual_correlation.h"
#include "censoring.h"
using namespace Rcpp;
namespace {
using namespace std;
using competingrisk::CensoringDistribution;
// Inclusive event-time risk sets, as in the teacher's residual and score blocks.
class ResidualDiagnostics {
  int n,P;
  vector<int> cc,event,d;
  vector<double> X,beta_estimate,CorrTime,CorrPValue,CorrRankTime,CorrRankPValue;
  vector<vector<double>> Z,schoenfeldLikeResiduals;
  const CensoringDistribution& kme;
  double calc_weight(double Time, double x_i) const
  {
    return kme.getSurvivalProb(Time) / kme.getSurvivalProb(min(Time, x_i));
  }

  int Yi_t(double Xi, double Time, int Eventi) const {
    if ((Xi >= Time) || ((Xi <= Time) && (Eventi > 1))) {
      return 1;
    } else {
      return 0;
    }
  }

  double alphaHat(double Time) const
  {
    std::vector<int> Yit(n, 0);
    for (int i = 0; i < n; ++i)
    {
      Yit[i] = Yi_t(X[i], Time, event[i]);
    }
    double num = 0, denom = 0;
    for (int i = 0; i < n; ++i)
    {
      if (d[i] == 0)
      {
        num += cc[i] * (1 - d[i]) * Yit[i];
        denom += (1 - d[i]) * Yit[i];
      }
    }
    return num / denom;
  }

  double rho(int di, int cci, double alpha_t) const
  {
    if (di == 0 && cci == 0)
    {
      return 0.0;
    }
    else if (di == 1)
    {
      return 1.0;
    }
    else if (di == 0 && cci == 1)
    {
      return 1.0 / alpha_t;
    }
    return 0.0;
  }

  void computeSchoenfeldLikeResiduals() {
    schoenfeldLikeResiduals.assign(n,std::vector<double>(P,NA_REAL));
    for(int i=0;i<n;++i) {
      competingrisk::interrupt_at(i);
      if(!d[i])continue;
      double alpha=alphaHat(X[i]),maximum=-INFINITY;
      std::vector<int> risk;std::vector<double> logits;
      bool valid=true;
      for(int j=0;j<n;++j)if((cc[j] || d[j]) && Yi_t(X[j],X[i],event[j])) {
        double sampling=rho(d[j],cc[j],alpha),w=calc_weight(X[i],X[j]);
        if(sampling==0 || w==0)continue;
        if(!(sampling>0) || !(w>0) || !std::isfinite(sampling) || !std::isfinite(w)) {valid=false;break;}
        double value=std::log(sampling)+std::log(w);
        for(int k=0;k<P;++k)value+=(Z[j][k]-Z[i][k])*beta_estimate[k];
        if(!std::isfinite(value)){valid=false;break;}
        risk.push_back(j);logits.push_back(value);maximum=std::max(maximum,value);
      }
      if(!valid || risk.empty())continue;
      double total=0;std::vector<double> contrast(P);
      for(std::size_t pos=0;pos<risk.size();++pos) {
        double w=std::exp(logits[pos]-maximum);total+=w;
        for(int k=0;k<P;++k)contrast[k]+=w*(Z[i][k]-Z[risk[pos]][k]);
      }
      if(!(total>0) || !std::isfinite(total))continue;
      for(int k=0;k<P;++k)schoenfeldLikeResiduals[i][k]=contrast[k]/total;
    }
  }

  void computeCorrelation(const std::vector<double>& times,
                          const std::vector<std::vector<double>>& residuals,
                          const std::vector<int>& target) {
    NumericMatrix matrix(n,P);
    for(int i=0;i<n;++i)for(int k=0;k<P;++k)matrix(i,k)=residuals[i][k];
    List value=competingrisk::event_correlations(wrap(times),matrix,wrap(target));
    CorrTime=as<vector<double>>(value["cor"]);CorrPValue=as<vector<double>>(value["cor_p"]);
    CorrRankTime=as<vector<double>>(value["rank_cor"]);CorrRankPValue=as<vector<double>>(value["rank_cor_p"]);
  }

public:
  ResidualDiagnostics(NumericVector time,IntegerVector status,NumericMatrix z,
                      IntegerVector subcohort,NumericVector beta,const CensoringDistribution& km)
    : n(time.size()),P(z.ncol()),cc(as<vector<int>>(subcohort)),event(as<vector<int>>(status)),
      d(n),X(as<vector<double>>(time)),beta_estimate(as<vector<double>>(beta)),
      Z(n,vector<double>(P)),kme(km) {
    for(int i=0;i<n;++i){d[i]=event[i]==1;for(int j=0;j<P;++j)Z[i][j]=z(i,j);}
    computeSchoenfeldLikeResiduals();
    computeCorrelation(X,schoenfeldLikeResiduals,d);
  }
  List result() const {
    NumericMatrix residuals(n,P);
    for(int i=0;i<n;++i)for(int j=0;j<P;++j)residuals(i,j)=schoenfeldLikeResiduals[i][j];
    return List::create(_["residuals"]=residuals,_["cor"]=CorrTime,_["cor_p"]=CorrPValue,
                        _["rank_cor"]=CorrRankTime,_["rank_cor_p"]=CorrRankPValue);
  }
};
}

// [[Rcpp::export]]
List residual_diagnostics(NumericVector time,IntegerVector status,NumericMatrix z,
                          IntegerVector cc,NumericVector beta,SEXP censoring) {
  XPtr<competingrisk::CensoringDistribution> km(censoring);
  return ResidualDiagnostics(time,status,z,cc,beta,*km).result();
}

// Internal shared correlation stage; also permits direct fault-injection validation.
// [[Rcpp::export]]
List residual_correlations(NumericVector time,NumericMatrix residuals,IntegerVector target) {
  return competingrisk::event_correlations(time,residuals,target);
}
