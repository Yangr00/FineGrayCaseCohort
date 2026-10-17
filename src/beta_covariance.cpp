#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <numeric>
#include <cmath>
#include "interrupt.h"
#include "censoring.h"
using namespace Rcpp;
namespace {
// Algebraically identical teacher covariance: strict risk sets in S0/S1/S2,
// inclusive risk sets in eta and alpha. No pair-by-predictor temporary arrays.
class BetaCovariance {
  int n,p;
  NumericVector time,beta;
  IntegerVector event,cc;
  NumericMatrix z,original;
  std::vector<int> order,events;
  std::vector<double> g,ex,alpha,tilde,s0,mean,eta,psi,mass,scale;
  const competingrisk::CensoringDistribution& km;
  bool observed(int i) const {return cc[i] || event[i]==1;}
  bool inclusive(int i,double t) const {return time[i]>=t || event[i]>1;}
  bool strict(int i,double t) const {return time[i]>t || (time[i]<=t && event[i]>1);}
  double weight(int i,int j) const {return time[j]<=time[i]?1.0:g[j]/g[i];}
  double rho(int i,double a) const {return event[i]==1?1.0:(cc[i]?1.0/a:0.0);}
public:
  BetaCovariance(NumericVector t,IntegerVector e,NumericMatrix x,IntegerVector c,
                 NumericVector b,const competingrisk::CensoringDistribution& k)
    :n(t.size()),p(x.ncol()),time(t),beta(b),event(e),cc(c),z(clone(x)),original(x),
     order(n),g(n),ex(n),alpha(n),tilde(n),s0(n),mean(n*p),eta(n*p),psi(n*p),mass(n),scale(p),km(k) {
    // Centered, dimensionless coordinates prevent offset cancellation and
    // ill-conditioned matrix solves caused solely by measurement units.
    for(int k=0;k<p;++k) {
      long double center=0,ss=0;int m=0;
      for(int i=0;i<n;++i)if(observed(i)){++m;center+=((long double)x(i,k)-center)/m;}
      for(int i=0;i<n;++i)if(observed(i)){long double d=(long double)x(i,k)-center;ss+=d*d;}
      scale[k]=std::sqrt(ss/m);
      if(!(scale[k]>0) || !std::isfinite(scale[k]))stop("Invalid covariate scale in beta covariance.");
      for(int i=0;i<n;++i)z(i,k)=observed(i)?static_cast<double>(((long double)x(i,k)-center)/scale[k]):0;
    }
    double maximum=-INFINITY;
    for(int i=0;i<n;++i)if(observed(i)) {
      for(int k=0;k<p;++k)ex[i]+=z(i,k)*(beta[k]*scale[k]);
      if(!std::isfinite(ex[i]))stop("Non-finite beta covariance linear predictor.");
      maximum=std::max(maximum,ex[i]);
    }
    for(int i=0;i<n;++i)if(observed(i))ex[i]=std::exp(ex[i]-maximum);
    std::iota(order.begin(),order.end(),0);
    std::stable_sort(order.begin(),order.end(),[&](int i,int j){return time[i]<time[j];});
    double alpha0=std::accumulate(cc.begin(),cc.end(),0.0)/n;
    for(int i=0;i<n;++i) {
      competingrisk::interrupt_at(i);
      g[i]=km.getSurvivalProb(time[i]);
      if(event[i]==1)events.push_back(i);
      double num=0,den=0;
      for(int j=0;j<n;++j)if(event[j]!=1 && inclusive(j,time[i])){num+=cc[j];den+=1;}
      alpha[i]=num/den;tilde[i]=alpha[i]==0?alpha0:alpha[i];
    }
    // Cache S0 and S1 once per target event. S2 is accumulated by exchanging
    // finite summations, reducing O(events*n*p^2) work to O(n*p^2).
    for(int j:events) {
      competingrisk::interrupt_at(j);
      for(int i=0;i<n;++i)if(observed(i) && strict(i,time[j])) {
        double w=rho(i,alpha[j])*weight(i,j)*ex[i];s0[j]+=w;
        for(int k=0;k<p;++k)mean[j*p+k]+=w*z(i,k);
      }
      if(!std::isfinite(s0[j]) || s0[j]<0)stop("Invalid weighted beta covariance risk denominator.");
      if(s0[j]==0) {
        for(int i=0;i<n;++i)if(observed(i) && strict(i,time[j]))
          stop("Weighted beta covariance risk mass underflowed.");
      }
      if(s0[j]!=0)for(int k=0;k<p;++k)mean[j*p+k]/=s0[j];
      s0[j]/=n;
    }
    for(int i=0;i<n;++i)if(observed(i)) {
      competingrisk::interrupt_at(i);
      if(event[i]==1)for(int k=0;k<p;++k)
        eta[i*p+k]=s0[i]==0?original(i,k)/scale[k]:z(i,k)-mean[i*p+k];
      for(int j:events)if(inclusive(i,time[j])) {
        double w=(s0[j]==0?0:1/s0[j])*ex[i]*weight(i,j)*rho(i,tilde[j])/n;
        for(int k=0;k<p;++k)eta[i*p+k]-=w*(z(i,k)-mean[j*p+k]);
        if(strict(i,time[j]))mass[i]+=w;
      }
    }
    // In q(X_i), X_l<X_i<=X_j and Y_l(X_j)=1 imply that l is a
    // subcohort competing event. Prefix competing-event moments and suffix
    // target-event moments therefore replace the original triple sum exactly.
    std::vector<double> suffix0(n+1),suffix1((n+1)*p),prefix1(p);
    for(int pos=n-1;pos>=0;--pos) {
      competingrisk::interrupt_at(pos);
      int j=order[pos];double w=event[j]==1 && s0[j]!=0?g[j]/(tilde[j]*s0[j]):0;
      suffix0[pos]=suffix0[pos+1]+w;
      for(int k=0;k<p;++k)suffix1[pos*p+k]=suffix1[(pos+1)*p+k]+w*mean[j*p+k];
    }
    std::vector<double> q(n*p),pi(n),fixedrho(n);
    for(int i=0;i<n;++i)fixedrho[i]=rho(i,tilde[i]);
    double at_risk=0;
    for(int end=n;end>0;) {
      Rcpp::checkUserInterrupt();
      int begin=end-1;while(begin>0 && time[order[begin-1]]==time[order[end-1]])--begin;
      for(int pos=begin;pos<end;++pos)at_risk+=fixedrho[order[pos]];
      for(int pos=begin;pos<end;++pos)pi[order[pos]]=at_risk/n;
      end=begin;
    }
    double prefix0=0;
    for(int begin=0;begin<n;) {
      Rcpp::checkUserInterrupt();
      int end=begin+1;while(end<n && time[order[end]]==time[order[begin]])++end;
      for(int pos=begin;pos<end;++pos) {
        int i=order[pos];if(event[i]==0 && cc[i])for(int k=0;k<p;++k)
          q[i*p+k]=(prefix1[k]*suffix0[begin]-prefix0*suffix1[begin*p+k])/(static_cast<double>(n)*n);
      }
      for(int pos=begin;pos<end;++pos) {
        int i=order[pos];if(event[i]>1 && cc[i]) {
          double w=ex[i]/g[i];prefix0+=w;
          for(int k=0;k<p;++k)prefix1[k]+=w*z(i,k);
        }
      }
      begin=end;
    }
    std::vector<double> cumulative(p);
    for(int begin=0;begin<n;) {
      Rcpp::checkUserInterrupt();
      int end=begin+1;while(end<n && time[order[end]]==time[order[begin]])++end;
      for(int pos=begin;pos<end;++pos) {
        int i=order[pos];if(event[i]==0 && cc[i])for(int k=0;k<p;++k)
          cumulative[k]+=q[i*p+k]/(pi[i]*pi[i]*n);
      }
      for(int pos=begin;pos<end;++pos) {
        int i=order[pos];if(event[i]==0 && cc[i])for(int k=0;k<p;++k)
          psi[i*p+k]=q[i*p+k]/pi[i]-cumulative[k];
      }
      begin=end;
    }
  }
  List result(bool details) const {
    NumericMatrix sigma(p,p),omega(p,p);
    for(int k=0;k<p;++k)for(int l=0;l<p;++l) {
      competingrisk::interrupt_at(k*p+l);
      for(int i=0;i<n;++i) {
        sigma(k,l)+=(eta[i*p+k]+psi[i*p+k])*(eta[i*p+l]+psi[i*p+l])/n;
        omega(k,l)+=mass[i]*z(i,k)*z(i,l)/n;
      }
      for(int j:events)if(s0[j]!=0)omega(k,l)-=mean[j*p+k]*mean[j*p+l]/n;
    }
    Function solve=Environment::base_env()["solve"];
    NumericMatrix inverse=solve(omega),tmp(p,p),variance(p,p),naive(p,p);
    for(int i=0;i<p;++i)for(int j=0;j<p;++j)for(int k=0;k<p;++k)tmp(i,j)+=inverse(i,k)*sigma(k,j);
    for(int i=0;i<p;++i)for(int j=0;j<p;++j) {
      for(int k=0;k<p;++k)variance(i,j)+=tmp(i,k)*inverse(k,j)/n;
      variance(i,j)/=scale[i]*scale[j];
      naive(i,j)=inverse(i,j)/(n*scale[i]*scale[j]);
    }
    List out=List::create(_["var"]=variance,_["naive.var"]=naive);
    if(details) {
      for(int k=0;k<p;++k)for(int l=0;l<p;++l) {
        sigma(k,l)*=scale[k]*scale[l];omega(k,l)*=scale[k]*scale[l];
      }
      out["score_covariance"]=sigma;
      NumericMatrix eta_out(n,p),psi_out(n,p);
      for(int i=0;i<n;++i)for(int k=0;k<p;++k) {
        eta_out(i,k)=eta[i*p+k]*scale[k];psi_out(i,k)=psi[i*p+k]*scale[k];
      }
      out["eta"]=eta_out;out["psi"]=psi_out;out["Omega"]=omega;
    }
    return out;
  }
};
}
// [[Rcpp::export]]
List beta_covariance(NumericVector time,IntegerVector status,NumericMatrix z,
                     IntegerVector cc,NumericVector beta,SEXP censoring,bool details=false) {
  XPtr<competingrisk::CensoringDistribution> km(censoring);
  return BetaCovariance(time,status,z,cc,beta,*km).result(details);
}
