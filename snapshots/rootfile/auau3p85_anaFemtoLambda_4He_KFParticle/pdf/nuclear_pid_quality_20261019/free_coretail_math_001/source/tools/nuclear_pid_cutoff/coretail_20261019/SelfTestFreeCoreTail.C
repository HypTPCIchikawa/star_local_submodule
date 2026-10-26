void SelfTestFreeCoreTail(const char*,const char*);
#ifndef __CINT__
#include "FreeGaussianCoreTail.h"
#include <iostream>
namespace FreeCoreTailTest {
using namespace NuclearFreeCoreTail;
unsigned long checks=0;
void Check(bool ok,const char* message){++checks;Require(ok,message);}
bool Near(double a,double b,double tolerance=1.e-11) {
  return std::fabs(a-b)<=tolerance*std::max(1.,std::max(std::fabs(a),std::fabs(b)));
}
bool Has(const std::vector<std::string>& flags,const char* flag) {
  return std::find(flags.begin(),flags.end(),flag)!=flags.end();
}
Request Seed() {
  Request r;r.nativeAxisMin=0;r.nativeAxisMax=120;r.backgroundWindowYield=20;
  Component a;a.name="he3";a.mean=28;a.sigma=2;a.startWindowYield=1000;
  a.nodes.push_back(NuclearFit::Node(5./18.,27.5,1.8));
  a.nodes.push_back(NuclearFit::Node(8./18.,28,2));
  a.nodes.push_back(NuclearFit::Node(5./18.,28.5,2.2));
  Component b;b.name="he4";b.mean=42;b.sigma=3;b.startWindowYield=500;
  r.components.push_back(a);r.components.push_back(b);
  for(unsigned k=0;k<100;++k)r.bins.push_back(NuclearFit::Bin(23+.25*k,23+.25*(k+1),10));
  return r;
}
NuclearFit::Result FakeModel(const Prepared& p,const Settings& cfg) {
  NuclearFit::Result model;model.components=p.request.components;
  model.finite=true;model.status=0;model.pointAccepted=true;model.amplitudeConverged=true;
  model.shapeCovarianceValid=true;model.yieldCovarianceValid=true;model.conditionalNdf=71;
  model.totalTemplates=NuclearFit::BuildTemplates(p.request,cfg,&model.windowAcceptance);
  const unsigned n=model.components.size(),columns=n+2;
  for(unsigned j=0;j<n;++j) {model.columnNames.push_back(model.components[j].name);model.columnWindowYields.push_back(model.components[j].startWindowYield);}
  model.columnNames.push_back("BG_left");model.columnNames.push_back("BG_right");
  model.columnWindowYields.push_back(10);model.columnWindowYields.push_back(10);
  model.yieldCovarianceColumnActive.assign(columns,true);
  model.fixedShapeYieldCovariance.assign(columns,Vec(columns,0));
  for(unsigned j=0;j<columns;++j)model.fixedShapeYieldCovariance[j][j]=model.columnWindowYields[j];
  model.fixedShapeYieldCovariance[0][1]=model.fixedShapeYieldCovariance[1][0]=.3;
  model.componentCounts.assign(model.totalTemplates.size(),Vec(n,0));
  model.totalCounts.assign(model.totalTemplates.size(),0);model.smoothCounts.assign(model.totalTemplates.size(),0);
  for(unsigned b=0;b<model.totalTemplates.size();++b)for(unsigned j=0;j<columns;++j) {
    const double count=model.columnWindowYields[j]*model.totalTemplates[b][j];
    model.totalCounts[b]+=count;
    if(j<n)model.componentCounts[b][j]=count;else model.smoothCounts[b]+=count;
  }
  model.shapeNames.push_back("test_dynamic_mean");model.shapePoint.push_back(28);
  model.shapeCovariance.assign(1,Vec(1,.2));model.meanError.assign(n,.1);model.sigmaError.assign(n,.2);
  return model;
}
}
#endif

void SelfTestFreeCoreTail(const char* settingsPath,const char* outputPath) {
#ifndef __CINT__
  using namespace FreeCoreTailTest;
  try {
    Require(gSystem->AccessPathName(outputPath),"Existing free core+tail result will not be replaced");
    const Settings cfg=NuclearFit::ReadSettings(settingsPath);
    const Shape shape=NuclearCoreTailFit::ReadShape(cfg),zero(0,0,1);
    const double infinity=std::numeric_limits<double>::infinity();
    const Request original=Seed();std::vector<Shape> zeroShapes(2,zero);
    const Prepared control=Prepare(original,zeroShapes,cfg);
    Check(control.request.components.size()==original.components.size(),"f=0 adds no columns");
    NuclearFit::Objective oldObjective(original,cfg),controlObjective(control.request,cfg);
    Check(oldObjective.Point()==controlObjective.Point()&&oldObjective.names==controlObjective.names,
          "f=0 exact old parameter vector and names");
    const Matrix oldTemplates=NuclearFit::BuildTemplates(original,cfg);
    Check(oldTemplates==NuclearFit::BuildTemplates(control.request,cfg),"f=0 exact old bin templates");
    for(unsigned j=0;j<original.components.size();++j) {
      const Component& a=original.components[j];const Component& b=control.request.components[j];
      Check(a.name==b.name&&a.mean==b.mean&&a.sigma==b.sigma&&a.startWindowYield==b.startWindowYield&&
            a.active==b.active&&a.freeMean==b.freeMean&&a.freeSigma==b.freeSigma&&a.nodes.size()==b.nodes.size(),
            "f=0 component field identity");
      for(unsigned q=0;q<a.nodes.size();++q)Check(a.nodes[q].weight==b.nodes[q].weight&&
          a.nodes[q].mean==b.nodes[q].mean&&a.nodes[q].sigma==b.nodes[q].sigma,"f=0 node identity");
    }
    std::vector<Shape> profiles(2,shape);const Prepared prepared=Prepare(original,profiles,cfg);
    Check(prepared.request.components.size()==4,"Two nonzero-tail isotopes have four independent yield columns");
    Check(prepared.isotopeMap[0].core==0&&prepared.isotopeMap[0].tail==1&&
          prepared.isotopeMap[1].core==2&&prepared.isotopeMap[1].tail==3,"Stable named isotope index map");
    Check(prepared.request.components[0].name=="he3_core"&&prepared.request.components[1].name=="he3_tail"&&
          prepared.request.components[2].name=="he4_core"&&prepared.request.components[3].name=="he4_tail",
          "Independent components retain named isotope ownership");
    std::vector<Shape> mixedProfiles;mixedProfiles.push_back(shape);mixedProfiles.push_back(zero);
    const Prepared mixedControl=Prepare(original,mixedProfiles,cfg);
    Check(mixedControl.request.components.size()==3&&mixedControl.isotopeMap[1].core==2&&
          mixedControl.isotopeMap[1].tail==-1&&mixedControl.request.components[2].name=="he4",
          "Mixed control/free index map has no phantom Gaussian-control tail");
    Request inactive=original;inactive.components[1].active=false;
    const Prepared inactiveControl=Prepare(inactive,profiles,cfg);
    Check(inactiveControl.request.components.size()==3&&!inactiveControl.isotopeMap[1].active&&
          inactiveControl.isotopeMap[1].tail==-1,"Inactive isotope does not add unsupported free-tail columns");
    const double lo=original.bins.front().lo,hi=original.bins.back().hi;
    for(unsigned s=0;s<2;++s) {
      const IsotopeColumns& map=prepared.isotopeMap[s];
      const Component& core=prepared.request.components[map.core];
      const Component& tail=prepared.request.components[unsigned(map.tail)];
      const double ac=NuclearFit::KernelIntegral(core,lo,hi),at=NuclearFit::KernelIntegral(tail,lo,hi);
      const double nc=core.startWindowYield/ac,nt=tail.startWindowYield/at;
      Check(Near(core.startWindowYield+tail.startWindowYield,original.components[s].startWindowYield),
            "Seed preserves total FIT-WINDOW isotope yield");
      Check(Near(nt/(nc+nt),shape.fraction),"Seed fraction is FULL-LINE, not window fraction");
      Check(std::fabs(tail.startWindowYield/(core.startWindowYield+tail.startWindowYield)-shape.fraction)>1.e-5,
            "Finite fit window distinguishes local and full-line fractions");
      Check(Near(NuclearFit::KernelIntegral(core,-infinity,infinity),1)&&
            Near(NuclearFit::KernelIntegral(tail,-infinity,infinity),1),"Each column is a normalized Gaussian quadrature");
      for(unsigned q=0;q<core.nodes.size();++q)
        Check(core.nodes[q].weight==tail.nodes[q].weight,"Core/tail preserve independent normalized R quadrature weights");
      double sum=0;
      for(unsigned b=0;b<original.bins.size();++b)sum+=NuclearFit::KernelIntegral(tail,original.bins[b].lo,original.bins[b].hi);
      Check(Near(sum,at),"Independent tail native-bin closure");
    }
    NuclearFit::Objective freeObjective(prepared.request,cfg);
    Check(freeObjective.names.size()==2*oldObjective.names.size(),"Independent core/tail means and sigmas count in NDF");
    Vec point=freeObjective.Point();const Vec before=point;
    for(unsigned k=0;k<freeObjective.species.size();++k)if(freeObjective.species[k]==1)
      point[k]+=freeObjective.kind[k]==0?1.25:std::log(1.3);
    const Request moved=freeObjective.Model(&point[0]);
    Check(moved.components[0].mean==prepared.request.components[0].mean&&
          moved.components[0].sigma==prepared.request.components[0].sigma,"Tail motion does not tie the core");
    Check(Near(moved.components[1].mean,prepared.request.components[1].mean+1.25)&&
          Near(moved.components[1].sigma,prepared.request.components[1].sigma*1.3),"Tail mean and sigma independently move");
    Check(moved.components[2].mean==prepared.request.components[2].mean&&
          moved.components[3].sigma==prepared.request.components[3].sigma,"Other isotope remains unchanged");
    for(unsigned q=0;q<moved.components[1].nodes.size();++q) {
      Check(Near(moved.components[1].nodes[q].mean,prepared.request.components[1].nodes[q].mean+1.25),
            "Tail R-node displacement follows its OWN fitted mean");
      Check(Near(moved.components[1].nodes[q].sigma,prepared.request.components[1].nodes[q].sigma*1.3),
            "Tail R-node scale follows its OWN fitted width");
    }
    const NuclearFit::Result model=FakeModel(prepared,cfg);const Result screened=Screen(prepared,model,cfg);
    Check(screened.groupingSupported,"Valid interior synthetic model passes declared post-fit screens");
    for(unsigned s=0;s<2;++s) {
      Check(Near(screened.fullLineTailFraction[s],shape.fraction),"Recovered full-line fraction from fitted-window Y/A");
      Check(std::fabs(screened.windowTailFraction[s]-screened.fullLineTailFraction[s])>1.e-5,
            "Reported full-line and fit-window fractions remain distinct");
    }
    const NuclearFit::Result folded=FoldIsotopes(screened);
    Check(folded.components.size()==2&&folded.columnWindowYields.size()==4,"Display folds independent columns plus untouched backgrounds");
    Check(folded.conditionalNdf==model.conditionalNdf&&folded.shapeNames==model.shapeNames&&
          folded.shapePoint==model.shapePoint&&folded.shapeCovariance==model.shapeCovariance,
          "Display fold preserves dynamic NDF and covariance parameterization");
    Check(folded.totalCounts==model.totalCounts&&folded.smoothCounts==model.smoothCounts,
          "Display fold does not change total or background counts");
    for(unsigned b=0;b<model.totalCounts.size();++b) {
      Check(Near(folded.componentCounts[b][0],model.componentCounts[b][0]+model.componentCounts[b][1])&&
            Near(folded.componentCounts[b][1],model.componentCounts[b][2]+model.componentCounts[b][3]),
            "Per-bin isotope display sum closes");
      double total=0;
      for(unsigned j=0;j<folded.columnWindowYields.size();++j)total+=folded.columnWindowYields[j]*folded.totalTemplates[b][j];
      Check(Near(total,model.totalCounts[b]),"Folded display templates and yields reproduce dynamic total");
    }
    Check(Near(folded.fixedShapeYieldCovariance[0][0],model.fixedShapeYieldCovariance[0][0]+
          model.fixedShapeYieldCovariance[1][1]+2*model.fixedShapeYieldCovariance[0][1]),
          "Isotope covariance includes TWO core-tail cross terms");
    Check(Near(NuclearFit::KernelIntegral(folded.components[0],-infinity,infinity),1),
          "Folded display isotope density normalizes on full line");
    NuclearFit::Result bad=model;bad.components[1].sigma=bad.components[0].sigma;
    Result rejected=Screen(prepared,bad,cfg);
    Check(!rejected.groupingSupported&&Has(rejected.isotopeFlags[0],"tail_not_broader_than_core"),
          "Narrow/equal tail flagged, never swapped");
    Check(rejected.fullLineTailFraction[0]>=0&&rejected.windowTailFraction[0]>=0,
          "Rejected width still retains computable descriptive fractions");
    Check(rejected.model.components[1].sigma==bad.components[1].sigma,"Rejected model shape is retained without clipping");
    bad=model;bad.columnWindowYields[1]=100*bad.columnWindowYields[0];
    rejected=Screen(prepared,bad,cfg);
    Check(!rejected.groupingSupported&&Has(rejected.isotopeFlags[0],"core_not_full_line_majority"),
          "Full-line majority checked after fitting");
    bad=model;bad.components[0].mean=prepared.isotopeMap[1].referenceMean;
    rejected=Screen(prepared,bad,cfg);
    Check(!rejected.groupingSupported&&Has(rejected.isotopeFlags[0],"core_moved_into_neighbor_seed_cell"),
          "Named core crossing neighbor is unsupported, never relabelled");
    bad=model;bad.shapeCovarianceValid=false;
    Check(!Screen(prepared,bad,cfg).groupingSupported,"Invalid shape covariance is unsupported");
    bad=model;bad.yieldCovarianceColumnActive[1]=false;
    rejected=Screen(prepared,bad,cfg);
    Check(Has(rejected.isotopeFlags[0],"tail_yield_not_interior"),"Boundary tail yield is unsupported");
    bool caught=false;try {std::vector<Shape> shortShapes;Prepare(original,shortShapes,cfg);}catch(const std::exception&){caught=true;}
    Check(caught,"Missing isotope seed profile rejected");
    caught=false;try {Request fixed=original;fixed.components[0].freeMean=false;Prepare(fixed,profiles,cfg);}catch(const std::exception&){caught=true;}
    Check(caught,"Nonzero free tail on a fixed core rejected rather than silently tied");
    caught=false;try {std::vector<Shape> outside(2,Shape(.1,1000,2));Prepare(original,outside,cfg);}catch(const std::exception&){caught=true;}
    Check(caught,"Tiny tail acceptance rejected without silently clipping shape or yield");
    std::ofstream out(outputPath);Require(bool(out),"Cannot CREATE free core+tail math result");
    out<<"{\"passed\":true,\"checks\":"<<checks<<",\"fitInvoked\":false,\"noRealDataInput\":true,"
       <<"\"productionChanged\":false,\"scope\":\"Independent-column mathematics, control identity, full-line fractions, group display covariance and post-fit screens; NOT real fitting or isotope purity validation\"}\n";
    out.close();Require(out.good(),"Free core+tail math result close failed");
    std::cout<<"NUCLEAR_FREE_CORE_TAIL_SELFTEST_PASS checks="<<checks<<std::endl;
  }catch(const std::exception& error) {
    std::cerr<<"SelfTestFreeCoreTail: "<<error.what()<<std::endl;gSystem->Exit(2);
  }
#endif
}
