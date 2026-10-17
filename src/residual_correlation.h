#ifndef COMPETINGRISK_RESIDUAL_CORRELATION_H
#define COMPETINGRISK_RESIDUAL_CORRELATION_H
#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <numeric>
#include "interrupt.h"
namespace competingrisk {
inline double residual_pearson(const std::vector<double>& x,const std::vector<double>& y) {
  if(x.size()<2 || x.size()!=y.size())return NA_REAL;
  long double mx=0,my=0;
  for(std::size_t i=0;i<x.size();++i) {
    if(!std::isfinite(x[i]) || !std::isfinite(y[i]))return NA_REAL;
    mx+=((long double)x[i]-mx)/(i+1);my+=((long double)y[i]-my)/(i+1);
  }
  long double sx=0,sy=0;
  for(std::size_t i=0;i<x.size();++i) {
    sx=std::max(sx,std::abs((long double)x[i]-mx));
    sy=std::max(sy,std::abs((long double)y[i]-my));
  }
  if(!(sx>0) || !(sy>0))return NA_REAL;
  long double xx=0,yy=0,xy=0;
  for(std::size_t i=0;i<x.size();++i) {
    long double dx=((long double)x[i]-mx)/sx,dy=((long double)y[i]-my)/sy;
    xx+=dx*dx;yy+=dy*dy;xy+=dx*dy;
  }
  double r=static_cast<double>(xy/(std::sqrt(xx)*std::sqrt(yy)));
  if(!std::isfinite(r))return NA_REAL;
  // Clamp only finite roundoff at the mathematical boundary.
  return std::max(-1.0,std::min(1.0,r));
}
inline double residual_p_value(double r,int n) {
  if(n<=2 || !std::isfinite(r))return NA_REAL;
  if(std::abs(r)==1)return 0;
  double t=std::abs(r)*std::sqrt((n-2.0)/((1-r)*(1+r)));
  return 2*R::pt(t,n-2.0,false,false);
}
inline Rcpp::List event_correlations(Rcpp::NumericVector time,Rcpp::NumericMatrix residuals,
                                     Rcpp::IntegerVector target) {
  int n=time.size(),p=residuals.ncol();
  if(residuals.nrow()!=n || target.size()!=n)Rcpp::stop("Incompatible residual correlation dimensions.");
  std::vector<double> x,ranks;std::vector<int> rows;
  bool finite_time=true;
  for(int i=0;i<n;++i)if(target[i]==1) {
    rows.push_back(i);x.push_back(time[i]);finite_time=finite_time && std::isfinite(time[i]);
  }
  ranks.resize(x.size());
  if(finite_time) {
    std::vector<int> order(x.size());std::iota(order.begin(),order.end(),0);
    std::stable_sort(order.begin(),order.end(),[&](int i,int j){return x[i]<x[j];});
    for(std::size_t begin=0;begin<order.size();) {
      std::size_t end=begin+1;while(end<order.size() && x[order[end]]==x[order[begin]])++end;
      for(std::size_t pos=begin;pos<end;++pos)ranks[order[pos]]=(begin+1.0+end)/2;
      begin=end;
    }
  }
  Rcpp::NumericVector cor(p,NA_REAL),cp(p,NA_REAL),rank_cor(p,NA_REAL),rp(p,NA_REAL);
  if(finite_time)for(int k=0;k<p;++k) {
    interrupt_at(k);std::vector<double> y;for(int i:rows)y.push_back(residuals(i,k));
    cor[k]=residual_pearson(x,y);cp[k]=residual_p_value(cor[k],rows.size());
    rank_cor[k]=residual_pearson(ranks,y);rp[k]=residual_p_value(rank_cor[k],rows.size());
  }
  return Rcpp::List::create(Rcpp::_["cor"]=cor,Rcpp::_["cor_p"]=cp,
    Rcpp::_["rank_cor"]=rank_cor,Rcpp::_["rank_cor_p"]=rp);
}
}
#endif
