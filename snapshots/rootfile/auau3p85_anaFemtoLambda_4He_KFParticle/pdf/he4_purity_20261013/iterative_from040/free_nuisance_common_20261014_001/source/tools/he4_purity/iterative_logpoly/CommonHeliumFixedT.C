// Conditional sensitivity only: T mean/width are empirical FREE-fit values.
// The frozen main engine and the original iterative study remain untouched.
void CommonHeliumFixedT(const char* cacheRoot,const char* output,const char* settings,
                        const char* seeds,const char* targets);
#ifndef __CINT__
#include "IterativeHeliumFit.C"
#endif
void CommonHeliumFixedT(const char* cacheRoot,const char* output,const char* settings,
                        const char* seeds,const char* targets){
#ifndef __CINT__
  using namespace HeliumIterative;
  try{
    Settings cfg=ReadSettings(settings);Target target=ReadTarget(targets);Study study(cfg,target);
    if(study.phase!="common"||Value(cfg,"penaltyStrength")!=0)throw std::runtime_error("fixed empirical T is a COMMON-only zero-penalty sensitivity");
    Protect(TString(output),settings,seeds,targets);
    if(gSystem->Load("libMinuit2")<0)throw std::runtime_error("cannot load Minuit2");
    TFile cache(cacheRoot,"READ");if(cache.IsZombie())throw std::runtime_error("cannot open frozen from040 ROOT");
    LoadData(study,cache,ReadSeeds(seeds));
    for(unsigned j=0;j<study.rows.size();++j)study.rows[j].shapeFree[0]=false;
    CommonResult common=FitCommon(study);
    WriteProducts(study,TString(output),common);
    for(unsigned j=0;j<study.rows.size();++j){delete study.rows[j].native;delete study.rows[j].h;}
    std::cout<<"ITERATIVE_HELIUM_COMPLETE phase=common fixedEmpiricalFreeT=true slices="<<study.rows.size()<<std::endl;
  }catch(const std::exception& e){std::cerr<<"CommonHeliumFixedT: "<<e.what()<<std::endl;gSystem->Exit(2);}
#endif
}
