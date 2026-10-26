#ifndef STAR_NUCLEAR_CORE_TAIL_FIT_COMMON_H
#define STAR_NUCLEAR_CORE_TAIL_FIT_COMMON_H
// Study-only fixed-profile Gaussian core+tail wrapper. The original engine is unchanged.
// f, delta, rho are NOT fitted. Original isotope yield remains a single fit column.
#include "../../nuclear_pid_tof/NuclearFitCommon.h"
#include <limits>

namespace NuclearCoreTailFit {
typedef NuclearFit::Request Request;
typedef NuclearFit::Result Result;
typedef NuclearFit::Component Component;
typedef NuclearFit::Node Node;
typedef NuclearFit::Settings Settings;
using NuclearFit::Require;
using NuclearFit::Finite;

struct Shape {
  double fraction,offsetSeedSigma,widthRatio;
  Shape(double f,double delta,double rho):fraction(f),offsetSeedSigma(delta),widthRatio(rho){}
};

inline void ValidateShape(const Shape& shape) {
  Require(Finite(shape.fraction)&&shape.fraction>=0&&shape.fraction<.5,
          "Core+tail fraction must be finite in [0,0.5); the core is the majority component");
  Require(Finite(shape.offsetSeedSigma),"Nonfinite tail offset in seed-core-sigma units");
  Require(Finite(shape.widthRatio)&&shape.widthRatio>=1&&
          (shape.fraction==0||shape.widthRatio>1),
          "A nonzero tail must be strictly broader than the core; Gaussian control permits rho=1");
}

inline Shape ReadShape(const Settings& settings,const std::string& prefix="") {
  Shape shape(NuclearFit::Number(settings,prefix+"tailFraction"),
              NuclearFit::Number(settings,prefix+"tailOffsetSeedSigma"),
              NuclearFit::Number(settings,prefix+"tailWidthRatio"));
  ValidateShape(shape);return shape;
}

inline void ValidateNodes(const Component& component) {
  Require(Finite(component.mean)&&Finite(component.sigma)&&component.sigma>0,
          "Nonfinite/nonpositive Gaussian core seed");
  double sum=0;
  for(unsigned q=0;q<component.nodes.size();++q) {
    const Node& node=component.nodes[q];
    Require(Finite(node.weight)&&node.weight>=0&&Finite(node.mean)&&
            Finite(node.sigma)&&node.sigma>0,"Invalid core/R-quadrature node");
    sum+=node.weight;
  }
  Require(component.nodes.empty()||std::fabs(sum-1)<1.e-9,
          "Core R-quadrature weights must sum to one");
}

inline Component BuildMixtureComponent(const Component& core,const Shape& shape) {
  ValidateShape(shape);ValidateNodes(core);
  // Exact Gaussian delegation: preserve the original nodes, arithmetic and fit parameter count.
  if(shape.fraction==0)return core;
  Component mixed=core;mixed.nodes.clear();
  const unsigned count=core.nodes.empty()?1:core.nodes.size();
  for(unsigned q=0;q<count;++q) {
    const Node seed=core.nodes.empty()?Node(1,core.mean,core.sigma):core.nodes[q];
    const double tailMean=seed.mean+shape.offsetSeedSigma*seed.sigma;
    const double tailSigma=shape.widthRatio*seed.sigma;
    Require(Finite(tailMean)&&Finite(tailSigma)&&tailSigma>0,"Tail seed overflows finite support");
    mixed.nodes.push_back(Node(seed.weight*(1-shape.fraction),seed.mean,seed.sigma));
    mixed.nodes.push_back(Node(seed.weight*shape.fraction,tailMean,tailSigma));
  }
  // The scalar mean/sigma are the CORE reference parameters. Existing Objective.Model applies
  // a common mean shift and width scale to every node. Therefore the initial separation
  // delta*SEED_sigma stays fixed in dE/dx during width fitting. It is not delta*FITTED_sigma.
  ValidateNodes(mixed);return mixed;
}

inline Request BuildMixtureRequest(const Request& core,const std::vector<Shape>& shapes) {
  Require(shapes.size()==core.components.size(),"Supply one explicit tail profile per isotope");
  Request mixed=core;
  for(unsigned s=0;s<core.components.size();++s)
    mixed.components[s]=BuildMixtureComponent(core.components[s],shapes[s]);
  return mixed;
}

inline Result Fit(const Request& core,const Settings& settings,const std::vector<Shape>& shapes) {
  return NuclearFit::Fit(BuildMixtureRequest(core,shapes),settings);
}

inline unsigned RNodeCount(const Component& mixed,const Shape& shape) {
  ValidateShape(shape);
  if(shape.fraction==0)return mixed.nodes.empty()?1:mixed.nodes.size();
  Require(!mixed.nodes.empty()&&mixed.nodes.size()%2==0,
          "Core+tail nodes must be interleaved [core(q),tail(q)]");
  return mixed.nodes.size()/2;
}

inline double RNodeWeight(const Component& mixed,const Shape& shape,unsigned q) {
  Require(q<RNodeCount(mixed,shape),"R-node index outside mixture");
  if(shape.fraction==0)return mixed.nodes.empty()?1:mixed.nodes[q].weight;
  const double weight=mixed.nodes[2*q].weight+mixed.nodes[2*q+1].weight;
  Require(Finite(weight)&&weight>=0,"Invalid combined R-node weight");
  if(weight>0)Require(std::fabs(mixed.nodes[2*q+1].weight/weight-shape.fraction)<1.e-12,
                      "Wrong fixed tail profile for this expanded component");
  return weight;
}

inline Node CoreNode(const Component& mixed,const Shape& shape,unsigned q) {
  Require(q<RNodeCount(mixed,shape),"Core-node index outside mixture");
  const Node node=shape.fraction==0?(mixed.nodes.empty()?Node(1,mixed.mean,mixed.sigma):mixed.nodes[q]):mixed.nodes[2*q];
  return Node(RNodeWeight(mixed,shape,q),node.mean,node.sigma);
}

inline double ComponentIntegral(const Component& mixed,const Shape& shape,unsigned q,
                                double lo,double hi) {
  // Conditional on R node q: NO R-quadrature weight here. Caller applies it once.
  Require(!std::isnan(lo)&&!std::isnan(hi)&&hi>=lo,"Invalid mixture integration bounds");
  const double weight=RNodeWeight(mixed,shape,q);
  if(!mixed.active||hi==lo||weight==0)return 0;
  const Node core=CoreNode(mixed,shape,q);
  if(shape.fraction==0)return NuclearFit::GaussianIntegral(lo,hi,core.mean,core.sigma);
  const Node& tail=mixed.nodes[2*q+1];
  const double result=(1-shape.fraction)*NuclearFit::GaussianIntegral(lo,hi,core.mean,core.sigma)+
                      shape.fraction*NuclearFit::GaussianIntegral(lo,hi,tail.mean,tail.sigma);
  Require(Finite(result)&&result>=0&&result<=1+1.e-12,"Invalid normalized mixture-bin integral");
  return result;
}

inline double KernelIntegral(const Component& mixed,double lo,double hi) {
  // Full R-averaged integral. Expanded nodes already contain all quadrature weights.
  return NuclearFit::KernelIntegral(mixed,lo,hi);
}

struct Moments {double mean,sigma;Moments(double m,double s):mean(m),sigma(s){} };
inline Moments MixtureMoments(const Component& mixed) {
  ValidateNodes(mixed);
  if(mixed.nodes.empty())return Moments(mixed.mean,mixed.sigma);
  double weight=0,mean=0;
  for(unsigned q=0;q<mixed.nodes.size();++q) {
    weight+=mixed.nodes[q].weight;mean+=mixed.nodes[q].weight*mixed.nodes[q].mean;
  }
  Require(Finite(weight)&&weight>0,"Cannot calculate zero-weight mixture moments");mean/=weight;
  double variance=0;
  for(unsigned q=0;q<mixed.nodes.size();++q) {
    const Node& node=mixed.nodes[q];const double shift=node.mean-mean;
    variance+=node.weight*(node.sigma*node.sigma+shift*shift)/weight;
  }
  Require(Finite(mean)&&Finite(variance)&&variance>0,"Nonfinite mixture moments");
  return Moments(mean,std::sqrt(variance));
}

inline const char* ResponseMeaning() {
  return "gaussian_core_tail_fixed; isotope density=(1-f)G(core)+fG(tail); f is full-line fraction; "
         "rho=tail/core width ratio; delta is initial tail displacement in each SEED core sigma, "
         "fixed in dE/dx during width profiling, not fitted-sigma units; majority narrow core; "
         "returned mean/sigma and their errors describe the core, NOT PID calibration sigma, "
         "mixture moments, intrinsic detector width, or isotope identity; fixed profiles add no "
         "continuous fitted parameters; grid/model selection is not included in conditionalNdf "
         "and cannot be used as an uncorrected goodness-of-fit p-value; no production mutation";
}
}
#endif
