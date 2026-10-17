#include <Rcpp.h>
#include <vector>
#include <numeric>
#include <algorithm>
#include <cmath>
#include <memory>
#include "interrupt.h"
#include "censoring.h"
using namespace Rcpp;
namespace {
using competingrisk::CensoringDistribution;
// Stable evaluation of the canonical unpenalized weighted partial likelihood.
// Standardization is a change of coordinates; estimates are returned in original units.
struct BetaRiskBlock {int event;std::vector<int> index;std::vector<double> w;};
struct BetaLikelihood {
  int p,n;
  std::vector<BetaRiskBlock> blocks;
  std::vector<double> scale,last;
  std::vector<double> design,centered;
  int last_order=-1;
  List cached;
  BetaLikelihood(NumericVector time,IntegerVector status,NumericMatrix z,IntegerVector cc,const CensoringDistribution& km):p(z.ncol()),n(time.size()),scale(p,0),design(n*p),centered(n*p) {
    int n=time.size(),m=0;
    std::vector<double> mean(p,0);
    for(int i=0;i<n;++i) {
      competingrisk::interrupt_at(i);
      if(cc[i] || status[i]==1){++m;for(int k=0;k<p;++k)mean[k]+=z(i,k);}
    }
    for(double& x:mean)x/=m;
    for(int i=0;i<n;++i)if(cc[i] || status[i]==1)for(int k=0;k<p;++k)scale[k]+=std::pow(z(i,k)-mean[k],2);
    for(double& x:scale){x=std::sqrt(x/m);if(!(x>0) || !std::isfinite(x))stop("Constant covariate.");}
    for(int i=0;i<n;++i)for(int k=0;k<p;++k) {
      design[i*p+k]=z(i,k);centered[i*p+k]=(z(i,k)-mean[k])/scale[k];
    }
    for(int i=0;i<n;++i)if(status[i]==1) {
      competingrisk::interrupt_at(i);
      double numerator=0,denominator=0;
      for(int j=0;j<n;++j)if(status[j]!=1 && (time[j]>=time[i] || status[j]>1)) {
        numerator+=cc[j];denominator+=1;
      }
      double alpha=numerator/denominator;
      BetaRiskBlock b;b.event=i;
      for(int j=0;j<n;++j)if((cc[j] || status[j]==1) && (time[j]>=time[i] || (time[j]<=time[i] && status[j]>1))) {
        double sampling=status[j]==1 ? 1 : cc[j]/alpha;
        double w=sampling*km.getSurvivalProb(time[i])/km.getSurvivalProb(std::min(time[i],time[j]));
        if(w==0)continue;
        if(!(w>0) || !std::isfinite(w))stop("Invalid weight.");
        b.index.push_back(j);b.w.push_back(std::log(w));
      }
      if(b.index.empty())stop("An event has an empty weighted risk set.");
      blocks.push_back(std::move(b));
    }
  }
  bool recession(const std::vector<double>& v) {
    bool strict=false;
    std::size_t iteration=0;
    for(const auto& b:blocks){competingrisk::interrupt_at(iteration++);for(int i:b.index){double s=0;for(int k=0;k<p;++k)s+=(design[i*p+k]-design[b.event*p+k])/scale[k]*v[k];if(s>0)return false;if(s<0)strict=true;}}
    return strict;
  }
  bool axis_separation() {
    for(int k=0;k<p;++k)for(double sign:{-1.0,1.0}){std::vector<double> v(p,0);v[k]=sign;if(recession(v))return true;}
    return false;
  }
  List eval(const std::vector<double>& beta_scaled,int order) {
    if(order<=last_order && beta_scaled==last)return cached;
    double loss=0;
    NumericVector gradient(p);NumericMatrix hessian(p,p);
    std::vector<double> linear(n),mass(n);
    std::vector<std::vector<double>> means;
    if(p>8)for(int i=0;i<n;++i)for(int k=0;k<p;++k)linear[i]+=centered[i*p+k]*beta_scaled[k];
    std::size_t iteration=0;
    for(const auto& b:blocks) {
      competingrisk::interrupt_at(iteration++);
      std::vector<double> logits=b.w,pr(logits.size());
      for(size_t j=0;j<logits.size();++j) {
        int i=b.index[j];
        if(p>8)logits[j]+=linear[i]-linear[b.event];
        else for(int k=0;k<p;++k)logits[j]+=(design[i*p+k]-design[b.event*p+k])/scale[k]*beta_scaled[k];
      }
      double maximum=*std::max_element(logits.begin(),logits.end()),total=0;
      if(!std::isfinite(maximum))stop("Non-finite likelihood during beta optimization.");
      for(size_t j=0;j<logits.size();++j)total+=pr[j]=std::exp(logits[j]-maximum);
      loss+=maximum+std::log(total);
      if(order==0)continue;
      std::vector<double> average(p,0);
      for(size_t j=0;j<logits.size();++j){pr[j]/=total;int i=b.index[j];for(int k=0;k<p;++k)average[k]+=pr[j]*(design[i*p+k]-design[b.event*p+k])/scale[k];}
      for(int k=0;k<p;++k)gradient[k]+=average[k];
      if(order<2)continue;
      if(p<=8) {
        for(size_t j=0;j<logits.size();++j)for(int k=0;k<p;++k)for(int l=0;l<p;++l)
          hessian(k,l)+=pr[j]*((design[b.index[j]*p+k]-design[b.event*p+k])/scale[k]-average[k])*
                            ((design[b.index[j]*p+l]-design[b.event*p+l])/scale[l]-average[l]);
      } else {
        for(size_t j=0;j<pr.size();++j)mass[b.index[j]]+=pr[j];
        for(int k=0;k<p;++k)average[k]+=centered[b.event*p+k];
        means.push_back(std::move(average));
      }
    }
    if(p>8 && order>=2)for(int k=0;k<p;++k)for(int l=0;l<p;++l) {
      competingrisk::interrupt_at(k*p+l);
      for(int i=0;i<n;++i)hessian(k,l)+=mass[i]*centered[i*p+k]*centered[i*p+l];
      for(const auto& m:means)hessian(k,l)-=m[k]*m[l];
    }
    last=beta_scaled;last_order=order;
    cached=List::create(_["value"]=loss,_["gradient"]=gradient,_["hessian"]=hessian);
    return cached;
  }
};

}

// Release the large cache on scope exit instead of waiting for R's GC trigger.
// [[Rcpp::export]]
void beta_likelihood_release(SEXP pointer) {
  auto* problem=static_cast<BetaLikelihood*>(R_ExternalPtrAddr(pointer));
  if(problem){delete problem;R_ClearExternalPtr(pointer);}
}

// [[Rcpp::export]]
List beta_likelihood_create(NumericVector time,IntegerVector status,NumericMatrix z,
                            IntegerVector cc,SEXP censoring) {
  XPtr<competingrisk::CensoringDistribution> km(censoring);
  std::unique_ptr<BetaLikelihood> owner(new BetaLikelihood(time,status,z,cc,*km));
  bool separated=owner->axis_separation();
  XPtr<BetaLikelihood> problem(owner.get(),true);
  owner.release();
  return List::create(_["pointer"]=problem,_["separated"]=separated);
}

// Objective, gradient and Hessian share the last-parameter evaluation cache.
// [[Rcpp::export]]
List beta_likelihood_eval(SEXP pointer,NumericVector beta_scaled) {
  XPtr<BetaLikelihood> problem(pointer);
  return problem->eval(as<std::vector<double>>(beta_scaled),2);
}

// [[Rcpp::export]]
List beta_likelihood_check(SEXP pointer,List result) {
  XPtr<BetaLikelihood> problem(pointer);
  int p=problem->p;NumericVector par=result["par"];
  const double score_tolerance=1e-7,step_tolerance=1e-6;
  List diagnostics=List::create(_["optimizer_run"]=true,
    _["optimizer_code"]=result["convergence"],_["optimizer_message"]=result["message"],
    _["iterations"]=result["iterations"],_["evaluations"]=result["evaluations"],
    _["scaled_parameters"]=par,_["n_events"]=static_cast<int>(problem->blocks.size()),
    _["score_tolerance"]=score_tolerance,_["step_tolerance"]=step_tolerance,
    _["verified"]=false);
  auto finish=[&](std::string status,std::string message,bool valid=false) {
    NumericVector beta(p,NA_REAL);
    if(valid)for(int k=0;k<p;++k)beta[k]=par[k]/problem->scale[k];
    diagnostics["verified"]=valid;
    return List::create(_["coefficients"]=beta,_["status"]=status,_["message"]=message,
                        _["solver_diagnostics"]=diagnostics);
  };
  if(par.size()!=p)return finish("nonconverged","nlminb returned incompatible coefficient dimensions.");
  std::vector<double> scaled=as<std::vector<double>>(par);
  for(double v:scaled)if(!std::isfinite(v))return finish("nonconverged","nlminb returned a non-finite coefficient.");
  if(problem->recession(scaled))return finish("no_finite_estimate",
    "No finite unpenalized estimate: a separating risk-set direction was detected.");
  List value=problem->eval(scaled,2);
  NumericVector gradient=value["gradient"];NumericMatrix information=value["hessian"];
  double objective=as<double>(value["value"]),max_score=0;
  diagnostics["objective"]=objective;diagnostics["scaled_score"]=gradient;
  if(!std::isfinite(objective))return finish("nonconverged","The final likelihood is non-finite.");
  for(double v:gradient) {
    if(!std::isfinite(v))return finish("nonconverged","The final score is non-finite.");
    max_score=std::max(max_score,std::abs(v));
  }
  for(double v:information)if(!std::isfinite(v))return finish("nonconverged","The final information matrix is non-finite.");
  Function eigen=Environment::base_env()["eigen"];
  List decomposition=eigen(information,Named("symmetric")=true);
  NumericVector values=decomposition["values"];NumericMatrix vectors=decomposition["vectors"];
  double floor=std::max(1e-12,std::max(0.0,values[0])*1e-12);
  diagnostics["information_min_eigenvalue"]=values[p-1];
  diagnostics["information_max_eigenvalue"]=values[0];
  diagnostics["information_condition"]=values[0]/values[p-1];
  diagnostics["max_score"]=max_score;
  double average_score=max_score/problem->blocks.size();
  diagnostics["average_score"]=average_score;
  if(!std::isfinite(values[p-1]) || values[p-1]<=floor)
    return finish("nonconverged","The final information matrix is singular or not positive definite.");
  NumericVector step(p);std::vector<double> direction(p);
  for(int j=0;j<p;++j) {
    double projection=0;for(int k=0;k<p;++k)projection+=vectors(k,j)*gradient[k];
    for(int k=0;k<p;++k)step[k]+=vectors(k,j)*projection/values[j];
  }
  double relative_step=0;
  for(int k=0;k<p;++k) {
    if(!std::isfinite(step[k]))return finish("nonconverged","The final Newton correction is non-finite.");
    direction[k]=-step[k];relative_step=std::max(relative_step,std::abs(step[k])/(1+std::abs(par[k])));
  }
  diagnostics["scaled_newton_step"]=step;diagnostics["relative_step"]=relative_step;
  if(problem->recession(direction))return finish("no_finite_estimate",
    "No finite unpenalized estimate: a separating risk-set direction was detected.");
  if(average_score>score_tolerance || relative_step>step_tolerance)
    return finish("nonconverged","nlminb failed the independent average-score or parameter-correction check.");
  return finish("converged",std::string("nlminb: ")+as<std::string>(result["message"])+
    "; verified on the same unpenalized likelihood.",true);
}
