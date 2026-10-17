#include <Rcpp.h>
#include <numeric>
#include <algorithm>
#include "interrupt.h"
#include "censoring.h"
using namespace Rcpp;

// Joint beta/theta influence covariance for beta(t)=beta+diag(q(t))*theta.
// Uses the manuscript's eta, psi (censoring), and mu (subcohort sampling)
// terms. Q multiplies both contributions to eta_theta, as required by
// differentiating U_theta; see the documented appendix interpretation.
// [[Rcpp::export]]
List regression_score_components(NumericVector time,IntegerVector event,
                                 NumericMatrix z,IntegerVector cc,
                                 NumericVector beta,NumericMatrix q,SEXP censoring,bool details=false) {
  int n=time.size(),p=z.ncol(),dim=2*p;
  if(n==0 || z.nrow()!=n || event.size()!=n || cc.size()!=n || beta.size()!=p || q.nrow()!=n || q.ncol()!=p)
    stop("Incompatible regression-test input dimensions.");
  std::vector<int> d(n),cens(n),events,observed,censored;
  std::vector<double> rho0(n),linear(n),ex(n);
  double alpha0=0,maximum=-INFINITY;
  for(int i=0;i<n;++i){
    d[i]=(event[i]==1);cens[i]=(event[i]==0);alpha0+=cc[i];
    if(d[i])events.push_back(i);
    if(cc[i] || d[i]){observed.push_back(i);if(cens[i])censored.push_back(i);
      for(int k=0;k<p;++k)linear[i]+=z(i,k)*beta[k];
      maximum=std::max(maximum,linear[i]);
    }
  }
  alpha0/=n;
  if(!(alpha0>0) || events.size()<3)stop("Regression test requires a nonempty subcohort and at least three target events.");
  for(int i:observed){rho0[i]=d[i]?1:cc[i]/alpha0;ex[i]=std::exp(linear[i]-maximum);}
  XPtr<competingrisk::CensoringDistribution> km(censoring);
  auto risk=[&](int i,double t){return time[i]>=t || event[i]>1;};
  auto alpha=[&](double t){double numerator=0,denominator=0;
    for(int i=0;i<n;++i)if(!d[i] && risk(i,t)){numerator+=cc[i];denominator+=1;}
    return (denominator>0 && numerator>0)?numerator/denominator:alpha0;
  };
  auto weight=[&](int i,double t){return t<=time[i]?1.0:km->getSurvivalProb(t)/km->getSurvivalProb(time[i]);};
  int ne=events.size();
  NumericMatrix mean(ne,p),ib(p,p),itb(p,p),eta(n,dim),mu(n,dim),psi(n,dim),censor_q(n,dim);
  NumericVector score(dim),den(ne),a(ne),pi(n),rho_censor(n);
  std::vector<double> probability(n),information_mass(n),theta_mass(n*p);
  for(int e=0;e<ne;++e){
    std::fill(probability.begin(),probability.end(),0.0);
    int j=events[e];double t=time[j];a[e]=alpha(t);
    for(int i:observed)if(risk(i,t)){
      double r=d[i]?1:cc[i]/a[e];double w=weight(i,t);
      if(!std::isfinite(w))stop("Non-finite censoring weight in regression test.");
      probability[i]=r*w*ex[i];den[e]+=probability[i];
      for(int k=0;k<p;++k)mean(e,k)+=probability[i]*z(i,k);
    }
    if(!(den[e]>0))stop("Empty weighted regression-test risk set.");
    for(int k=0;k<p;++k)mean(e,k)/=den[e];
    for(int i:observed)probability[i]/=den[e];
    for(int k=0;k<p;++k){double r=z(j,k)-mean(e,k);score[k]+=r;score[p+k]+=q(j,k)*r;}
    for(int i:observed){
      information_mass[i]+=probability[i]/n;
      for(int k=0;k<p;++k)theta_mass[i*p+k]+=q(j,k)*probability[i]/n;
    }
    double phi=0,g0=0;std::vector<double> g1(p,0);
    for(int i:observed)if(!d[i] && cc[i] && risk(i,t)){
      phi+=1/alpha0;double v=weight(i,t)*ex[i]/alpha0;g0+=v;
      for(int k=0;k<p;++k)g1[k]+=v*z(i,k);
    }
    for(int i:observed){
      double base=risk(i,t)?weight(i,t)*ex[i]/den[e]:0;
      for(int k=0;k<p;++k){
        double contrast=z(i,k)-mean(e,k);
        double component=(i==j?contrast:0)-base*contrast;
        eta(i,k)+=component;eta(i,p+k)+=q(j,k)*component;
        if(!d[i] && risk(i,t)){
          double adjustment=phi>0?(g1[k]-mean(e,k)*g0)/phi:0;
          double value=(weight(i,t)*ex[i]*contrast-adjustment)/den[e];
          mu(i,k)+=value;mu(i,p+k)+=q(j,k)*value;
        }
      }
    }
    checkUserInterrupt();
  }
  for(int k=0;k<p;++k)for(int l=0;l<p;++l){
    competingrisk::interrupt_at(k*p+l);
    for(int i:observed){
      ib(k,l)+=information_mass[i]*z(i,k)*z(i,l);
      itb(k,l)+=theta_mass[i*p+k]*z(i,k)*z(i,l);
    }
    for(int e=0;e<ne;++e){
      double v=mean(e,k)*mean(e,l)/n;
      ib(k,l)-=v;itb(k,l)-=q(events[e],k)*v;
    }
  }
  // The censoring contrast has Xi<Xl<=Xevent, so only sampled competing
  // events contribute. Sorted moments replace its cubic loop and avoid
  // retaining an events-by-cohort probability matrix.
  std::vector<int> order(n),event_at(n,-1);
  std::iota(order.begin(),order.end(),0);
  std::stable_sort(order.begin(),order.end(),[&](int i,int j){return time[i]<time[j];});
  for(int e=0;e<ne;++e)event_at[events[e]]=e;
  std::vector<double> suffix0(n+1),suffix1((n+1)*p),suffix_q((n+1)*p),suffix_qz((n+1)*p);
  for(int pos=n-1;pos>=0;--pos){
    competingrisk::interrupt_at(pos);
    int j=order[pos],e=event_at[j];
    double v=e>=0?km->getSurvivalProb(time[j])/(a[e]*den[e]):0;
    suffix0[pos]=suffix0[pos+1]+v;
    for(int k=0;k<p;++k){
      double m=e>=0?mean(e,k):0,qt=e>=0?q(j,k):0;
      suffix1[pos*p+k]=suffix1[(pos+1)*p+k]+v*m;
      suffix_q[pos*p+k]=suffix_q[(pos+1)*p+k]+v*qt;
      suffix_qz[pos*p+k]=suffix_qz[(pos+1)*p+k]+v*qt*m;
    }
  }
  double prefix0=0;std::vector<double> prefix1(p);
  for(int begin=0;begin<n;){
    Rcpp::checkUserInterrupt();
    int end=begin+1;while(end<n && time[order[end]]==time[order[begin]])++end;
    for(int pos=begin;pos<end;++pos){
      int l=order[pos];if(event[l]==0 && cc[l])for(int k=0;k<p;++k){
        censor_q(l,k)=(prefix1[k]*suffix0[begin]-prefix0*suffix1[begin*p+k])/n;
        censor_q(l,p+k)=(prefix1[k]*suffix_q[begin*p+k]-prefix0*suffix_qz[begin*p+k])/n;
      }
    }
    for(int pos=begin;pos<end;++pos){
      int i=order[pos];if(event[i]>1 && cc[i]){
        double v=ex[i]/km->getSurvivalProb(time[i]);prefix0+=v;
        for(int k=0;k<p;++k)prefix1[k]+=v*z(i,k);
      }
    }
    begin=end;
  }
  for(int l:censored){
    competingrisk::interrupt_at(l);
    double al=alpha(time[l]);rho_censor[l]=1/al;
    for(int i:observed)if(time[i]>=time[l])pi[l]+=(d[i]?1:cc[i]/al)/n;
    if(!(pi[l]>0))stop("Zero censoring-risk denominator in regression test.");
  }
  for(int i:observed){
    competingrisk::interrupt_at(i);
    for(int k=0;k<dim;++k){
      if(cens[i])psi(i,k)=censor_q(i,k)/pi[i];
      for(int l:censored)if(time[i]>=time[l])
        psi(i,k)-=censor_q(l,k)*rho_censor[l]/(n*pi[l]*pi[l]);
    }
  }
  NumericMatrix sigma(dim,dim);
  for(int i:observed){competingrisk::interrupt_at(i);for(int k=0;k<dim;++k)for(int l=0;l<dim;++l)
    sigma(k,l)+=rho0[i]*((eta(i,k)+psi(i,k))*(eta(i,l)+psi(i,l))+
                             (1-alpha0)/alpha0*mu(i,k)*mu(i,l))/n;
  }
  List out=List::create(_["score"]=score,_["Omega_beta_beta"]=ib,
                        _["Omega_theta_beta"]=itb,_["Sigma"]=sigma);
  if(details) {
    out["eta"]=eta;out["psi"]=psi;out["mu"]=mu;out["rho"]=rho0;out["alpha"]=alpha0;
  }
  return out;
}
