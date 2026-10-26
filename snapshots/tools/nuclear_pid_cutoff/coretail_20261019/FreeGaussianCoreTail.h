#ifndef STAR_NUCLEAR_FREE_GAUSSIAN_CORE_TAIL_H
#define STAR_NUCLEAR_FREE_GAUSSIAN_CORE_TAIL_H
// Study only. Independent Gaussian columns delegate to the immutable original fitter.
// Fits are never relabelled, clipped, or promoted to calibrated isotope purity here.
#include "NuclearCoreTailFitCommon.h"

namespace NuclearFreeCoreTail {
typedef NuclearFit::Request Request;
typedef NuclearFit::Component Component;
typedef NuclearFit::Settings Settings;
typedef NuclearFit::Vec Vec;
typedef NuclearFit::Matrix Matrix;
typedef NuclearCoreTailFit::Shape Shape;
using NuclearFit::Require;
using NuclearFit::Finite;

struct IsotopeColumns {
  std::string name;unsigned core;int tail;bool active;
  double referenceMean,referenceSigma;
  IsotopeColumns():core(0),tail(-1),active(false),referenceMean(0),referenceSigma(1){}
};
struct Prepared {
  Request request,original;
  std::vector<IsotopeColumns> isotopeMap;
};
struct Result {
  NuclearFit::Result model;
  std::vector<IsotopeColumns> isotopeMap;
  Vec fullLineTailFraction,windowTailFraction;
  std::vector<std::vector<std::string> > isotopeFlags;
  std::vector<std::string> flags;
  bool groupingSupported;
  Result():groupingSupported(false){}
};

inline Component TailSeed(const Component& core,const Shape& shape) {
  Component tail=core;tail.name=core.name+"_tail";
  tail.mean=core.mean+shape.offsetSeedSigma*core.sigma;
  tail.sigma=shape.widthRatio*core.sigma;
  for(unsigned q=0;q<tail.nodes.size();++q) {
    tail.nodes[q].mean=core.nodes[q].mean+shape.offsetSeedSigma*core.nodes[q].sigma;
    tail.nodes[q].sigma=shape.widthRatio*core.nodes[q].sigma;
  }
  NuclearCoreTailFit::ValidateNodes(tail);return tail;
}

inline Prepared Prepare(const Request& original,const std::vector<Shape>& shapes,
                        const Settings& cfg) {
  Require(shapes.size()==original.components.size(),"Explicit tail seed required for every original isotope");
  Require(!original.bins.empty(),"No fit bins for free core+tail preparation");
  const double lo=original.bins.front().lo,hi=original.bins.back().hi;
  const double minimum=NuclearFit::Number(cfg,"minimumGaussianAcceptance");
  Require(minimum>0&&minimum<1,"Invalid declared Gaussian acceptance floor");
  Prepared out;out.original=original;out.request=original;out.request.components.clear();
  for(unsigned s=0;s<original.components.size();++s) {
    const Component& seed=original.components[s];const Shape& shape=shapes[s];
    NuclearCoreTailFit::ValidateShape(shape);NuclearCoreTailFit::ValidateNodes(seed);
    Require(Finite(seed.startWindowYield)&&seed.startWindowYield>=0,"Invalid isotope window-yield seed");
    IsotopeColumns map;map.name=seed.name;map.active=seed.active;
    map.referenceMean=seed.mean;map.referenceSigma=seed.sigma;
    map.core=out.request.components.size();
    // f=0 retains the EXACT original representation, parameter count, and name.
    if(shape.fraction==0||!seed.active) {
      out.request.components.push_back(seed);out.isotopeMap.push_back(map);continue;
    }
    Require(seed.freeMean&&seed.freeSigma,"Nonzero free tail requires independently free core mean and width");
    Component core=seed,tail=TailSeed(seed,shape);core.name=seed.name+"_core";
    const double ac=NuclearFit::KernelIntegral(core,lo,hi);
    const double at=NuclearFit::KernelIntegral(tail,lo,hi);
    Require(Finite(ac)&&Finite(at)&&ac>minimum&&at>minimum,
            "Free core/tail seed acceptance outside explicit support; no clipping/deactivation");
    const double mix=(1-shape.fraction)*ac+shape.fraction*at;
    // The input isotope yield is in the FIT WINDOW; f is on the FULL LINE.
    core.startWindowYield=seed.startWindowYield*(1-shape.fraction)*ac/mix;
    tail.startWindowYield=seed.startWindowYield*shape.fraction*at/mix;
    out.request.components.push_back(core);map.tail=out.request.components.size();
    out.request.components.push_back(tail);out.isotopeMap.push_back(map);
  }
  return out;
}

inline Prepared BuildRequest(const Request& original,const std::vector<Shape>& shapes,
                             const Settings& cfg) {return Prepare(original,shapes,cfg);}

inline void Flag(std::vector<std::string>& flags,const std::string& value) {
  if(std::find(flags.begin(),flags.end(),value)==flags.end())flags.push_back(value);
}
inline bool InteriorYield(const NuclearFit::Result& model,unsigned j) {
  return j<model.columnWindowYields.size()&&Finite(model.columnWindowYields[j])&&
         model.columnWindowYields[j]>0&&j<model.yieldCovarianceColumnActive.size()&&
         model.yieldCovarianceColumnActive[j];
}
inline double FullLineYield(const NuclearFit::Result& model,unsigned j) {
  Require(j<model.columnWindowYields.size()&&j<model.windowAcceptance.size(),"Missing fitted-column acceptance/yield");
  const double y=model.columnWindowYields[j],a=model.windowAcceptance[j];
  Require(Finite(y)&&y>=0&&Finite(a)&&a>0,"Invalid fitted-column full-line normalization");
  return y/a;
}

inline Result Screen(const Prepared& prepared,const NuclearFit::Result& model,const Settings& cfg) {
  Result out;out.model=model;out.isotopeMap=prepared.isotopeMap;
  const unsigned n=out.isotopeMap.size();
  out.fullLineTailFraction.assign(n,-1);out.windowTailFraction.assign(n,-1);
  out.isotopeFlags.resize(n);
  if(!model.finite)Flag(out.flags,"nonfinite_model");
  if(model.status!=0)Flag(out.flags,"minimizer_status_nonzero");
  if(!model.pointAccepted)Flag(out.flags,"fitted_point_not_accepted");
  if(!model.amplitudeConverged)Flag(out.flags,"amplitude_profile_not_converged");
  if(model.atBoundary)Flag(out.flags,"shape_on_mathematical_boundary");
  if(!model.shapeNames.empty()&&!model.shapeCovarianceValid)Flag(out.flags,"shape_covariance_invalid");
  if(!model.yieldCovarianceValid)Flag(out.flags,"conditional_yield_covariance_invalid");
  if(model.components.size()!=prepared.request.components.size()) {
    Flag(out.flags,"model_component_dimension_mismatch");return out;
  }
  const double minimum=NuclearFit::Number(cfg,"minimumGaussianAcceptance");
  for(unsigned s=0;s<n;++s) {
    const IsotopeColumns& map=out.isotopeMap[s];std::vector<std::string>& flags=out.isotopeFlags[s];
    if(!map.active)continue;
    const Component& core=model.components[map.core];
    if(!InteriorYield(model,map.core))Flag(flags,"core_yield_not_interior");
    if(map.core>=model.windowAcceptance.size()||!Finite(model.windowAcceptance[map.core])||
       model.windowAcceptance[map.core]<=minimum)Flag(flags,"core_acceptance_unsupported");
    // Core labels retain their named seed's nearest-mean cell. No post-fit swapping.
    if(!Finite(core.mean)||!Finite(core.sigma)||core.sigma<=0)Flag(flags,"nonfinite_core_shape");
    for(unsigned other=0;other<n;++other)if(other!=s&&out.isotopeMap[other].active) {
      const double reference=out.isotopeMap[other].referenceMean;
      if(std::fabs(core.mean-reference)<std::fabs(core.mean-map.referenceMean))
        Flag(flags,"core_moved_into_neighbor_seed_cell");
    }
    if(map.tail<0) {out.fullLineTailFraction[s]=0;out.windowTailFraction[s]=0;continue;}
    const unsigned ti=unsigned(map.tail);const Component& tail=model.components[ti];
    if(!InteriorYield(model,ti))Flag(flags,"tail_yield_not_interior");
    if(ti>=model.windowAcceptance.size()||!Finite(model.windowAcceptance[ti])||
       model.windowAcceptance[ti]<=minimum)Flag(flags,"tail_acceptance_unsupported");
    if(!Finite(tail.mean)||!Finite(tail.sigma)||tail.sigma<=0)Flag(flags,"nonfinite_tail_shape");
    if(!(tail.sigma>core.sigma))Flag(flags,"tail_not_broader_than_core");
    // Keep descriptive fractions even when identity/width/interior screens fail.
    // Their availability is NOT an approval of the fitted isotope interpretation.
    if(map.core<model.columnWindowYields.size()&&ti<model.columnWindowYields.size()&&
       map.core<model.windowAcceptance.size()&&ti<model.windowAcceptance.size()&&
       Finite(model.columnWindowYields[map.core])&&model.columnWindowYields[map.core]>=0&&
       Finite(model.columnWindowYields[ti])&&model.columnWindowYields[ti]>=0&&
       Finite(model.windowAcceptance[map.core])&&model.windowAcceptance[map.core]>0&&
       Finite(model.windowAcceptance[ti])&&model.windowAcceptance[ti]>0) {
      const double nc=FullLineYield(model,map.core),nt=FullLineYield(model,ti);
      const double yc=model.columnWindowYields[map.core],yt=model.columnWindowYields[ti];
      if(Finite(nc+nt)&&nc+nt>0&&Finite(yc+yt)&&yc+yt>0) {
        out.fullLineTailFraction[s]=nt/(nc+nt);out.windowTailFraction[s]=yt/(yc+yt);
        if(!(out.fullLineTailFraction[s]<.5))Flag(flags,"core_not_full_line_majority");
      } else Flag(flags,"invalid_group_normalization");
    } else Flag(flags,"invalid_group_normalization");
  }
  out.groupingSupported=out.flags.empty();
  for(unsigned s=0;s<n;++s)out.groupingSupported=out.groupingSupported&&out.isotopeFlags[s].empty();
  return out;
}

inline Result Fit(const Prepared& prepared,const Settings& cfg) {
  return Screen(prepared,NuclearFit::Fit(prepared.request,cfg),cfg);
}
inline Result Fit(const Request& original,const Settings& cfg,const std::vector<Shape>& shapes) {
  return Fit(Prepare(original,shapes,cfg),cfg);
}

inline Matrix FoldCovariance(const Matrix& covariance,const std::vector<std::vector<unsigned> >& groups) {
  Matrix out(groups.size(),Vec(groups.size(),0));
  for(unsigned i=0;i<groups.size();++i)for(unsigned j=0;j<groups.size();++j)
    for(unsigned a=0;a<groups[i].size();++a)for(unsigned b=0;b<groups[j].size();++b) {
      const unsigned ia=groups[i][a],ib=groups[j][b];
      Require(ia<covariance.size()&&ib<covariance[ia].size(),"Invalid covariance group index");
      out[i][j]+=covariance[ia][ib];
    }
  return out;
}

inline NuclearFit::Result FoldIsotopes(const Result& result) {
  // DISPLAY adapter only: keeps the original 4-isotope histogram contract, not a new fit.
  // The original dynamic model/covariance MUST also be persisted. Never recalculate NDF.
  const NuclearFit::Result& model=result.model;NuclearFit::Result out=model;
  const unsigned dynamic=model.components.size(),n=result.isotopeMap.size();
  Require(model.columnWindowYields.size()>=dynamic,"Cannot fold incomplete fitted columns");
  std::vector<std::vector<unsigned> > groups;
  out.components.clear();out.columnNames.clear();out.columnWindowYields.clear();
  out.windowAcceptance.clear();out.meanError.clear();out.sigmaError.clear();
  out.effectiveCount.assign(n,0);out.conditionalDropDeltaD.assign(n,-1);
  for(unsigned s=0;s<n;++s) {
    const IsotopeColumns& map=result.isotopeMap[s];Require(map.core<dynamic,"Invalid isotope core index");
    std::vector<unsigned> group(1,map.core);
    if(map.tail>=0) {Require(unsigned(map.tail)<dynamic,"Invalid isotope tail index");group.push_back(unsigned(map.tail));}
    groups.push_back(group);Component combined=model.components[map.core];combined.name=map.name;
    double y=0,full=0;
    for(unsigned j=0;j<group.size();++j) {
      y+=model.columnWindowYields[group[j]];
      if(combined.active)full+=FullLineYield(model,group[j]);
    }
    if(map.tail>=0&&combined.active&&full>0) {
      combined.nodes.clear();
      for(unsigned j=0;j<group.size();++j) {
        const Component& component=model.components[group[j]];
        const double share=FullLineYield(model,group[j])/full;
        if(component.nodes.empty())combined.nodes.push_back(NuclearFit::Node(share,component.mean,component.sigma));
        else for(unsigned q=0;q<component.nodes.size();++q) {
          NuclearFit::Node node=component.nodes[q];node.weight*=share;combined.nodes.push_back(node);
        }
      }
    }
    combined.startWindowYield=y;out.components.push_back(combined);out.columnNames.push_back(map.name);
    out.columnWindowYields.push_back(y);out.windowAcceptance.push_back(full>0?y/full:0);
    out.meanError.push_back(map.core<model.meanError.size()?model.meanError[map.core]:-1);
    out.sigmaError.push_back(map.core<model.sigmaError.size()?model.sigmaError[map.core]:-1);
  }
  for(unsigned j=dynamic;j<model.columnWindowYields.size();++j) {
    groups.push_back(std::vector<unsigned>(1,j));out.columnWindowYields.push_back(model.columnWindowYields[j]);
    out.columnNames.push_back(j<model.columnNames.size()?model.columnNames[j]:"BG_unrecorded");
  }
  out.componentCounts.assign(model.componentCounts.size(),Vec(n,0));
  out.totalTemplates.assign(model.totalTemplates.size(),Vec(groups.size(),0));
  for(unsigned b=0;b<model.totalTemplates.size();++b)for(unsigned j=0;j<groups.size();++j) {
    double count=0;
    for(unsigned k=0;k<groups[j].size();++k) {
      const unsigned col=groups[j][k];Require(col<model.totalTemplates[b].size(),"Incomplete model template columns");
      count+=model.columnWindowYields[col]*model.totalTemplates[b][col];
    }
    out.totalTemplates[b][j]=out.columnWindowYields[j]>0?count/out.columnWindowYields[j]:0;
    if(j<n) {
      out.componentCounts[b][j]=count;
      if(b<model.totalCounts.size()&&model.totalCounts[b]>0)out.effectiveCount[j]+=count*count/model.totalCounts[b];
    }
  }
  out.yieldCovarianceColumnActive.assign(groups.size(),false);
  for(unsigned j=0;j<groups.size();++j) {
    bool active=true;
    for(unsigned k=0;k<groups[j].size();++k)active=active&&InteriorYield(model,groups[j][k]);
    out.yieldCovarianceColumnActive[j]=active;
  }
  if(model.yieldCovarianceValid)out.fixedShapeYieldCovariance=FoldCovariance(model.fixedShapeYieldCovariance,groups);
  else out.fixedShapeYieldCovariance.clear();
  // shapePoint/shapeNames/shapeCovariance/conditionalNdf remain DYNAMIC, not 4-isotope counts.
  return out;
}

inline const char* ResponseMeaning() {
  return "gaussian_coretail_free_components; each named isotope expands to independent Gaussian core/tail means, "
         "widths and nonnegative fit-window yields; f=0 delegates exactly to old Gaussian without tail columns; "
         "seed f and fitted majority refer to FULL-LINE yields Ywindow/acceptance, not local yield ratio; "
         "core-majority, broader-tail, interior/covariance and named-seed-neighbor checks are POST-FIT screens; "
         "no swapping or clipping; folded four-isotope histograms are display only, preserve true dynamic NDF; "
         "conditional covariance is not joint shape/yield covariance; grouped purity is UNAVAILABLE, not zero; "
         "no production PID mutation";
}
}
#endif
