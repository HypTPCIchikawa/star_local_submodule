// Fast diagnostic: refit cached original projections, not the multi-GB merged file.
void CheckLowRigidity(const char* cache,const char* settings,const char* output,const char* indicesCsv="");
#ifndef __CINT__
#include "FitHe4Purity.C"
#endif
void CheckLowRigidity(const char* cache,const char* settings,const char* output,const char* indicesCsv) {
#ifndef __CINT__
  using namespace He4Purity;
  try {
    const std::map<std::string,std::string> cfg=Settings(settings);
    const TString out(output);if(!gSystem->AccessPathName(out+"/prototype.json"))throw std::runtime_error("prototype output exists");gSystem->mkdir(out,kTRUE);
    TFile file(cache,"READ");if(file.IsZombie())throw std::runtime_error("cannot open cached projections");
    std::ofstream json((out+"/prototype.json").Data());json<<std::setprecision(17)<<"{\"cache\":"<<Quote(cache)<<",\"settings\":"<<Quote(settings)<<",\"slices\":[";
    std::vector<int> indices;
    if(indicesCsv&&*indicesCsv){
      std::istringstream values(indicesCsv);std::string value;
      while(std::getline(values,value,',')){
        std::istringstream token(value);int index=-1;char extra=0;token>>index;
        if(!token||token>>extra||index<0||index>=210||std::find(indices.begin(),indices.end(),index)!=indices.end())throw std::runtime_error("invalid/duplicate cached slice index");
        indices.push_back(index);
      }
      if(indices.empty())throw std::runtime_error("empty cached slice indices");
    }else {const int defaults[]={10,15,20,25,30,35,45,80,120};indices.assign(defaults,defaults+9);}
    TCanvas canvas("lowR_prototype","lowR prototype",1500,1000);gStyle->SetPaperSize(30,20);gStyle->SetOptStat(0);
    const TString book=out+"/prototype_linear.pdf";canvas.Print(book+"[");
    for(unsigned j=0;j<indices.size();++j){
      Slice s;const int index=indices[j];s.lo=.4+.01*index;s.hi=s.lo+.01;s.r=(s.lo+s.hi)/2;
      s.h=dynamic_cast<TH1D*>(file.Get(Form("projection_%03d",index)));if(!s.h)throw std::runtime_error("missing cached projection");s.h=static_cast<TH1D*>(s.h->Clone(Form("proto_%d",index)));s.h->SetDirectory(0);
      s.mean=NuclearIdDeDxVsMom::GetMean(NuclearIdDeDxVsMom::kHe4,s.r);s.sigma=NuclearIdDeDxVsMom::GetSigma(NuclearIdDeDxVsMom::kHe4,s.r);
      if(FeatureEnabled(cfg,"dedxRebinFactor")){
        const int factor=int(Required(cfg,"dedxRebinFactor"));
        const double physicalWidth=s.h->GetBinWidth(1)*s.sigma;
        const double targetWidth=.1*factor;
        if(factor<1)throw std::runtime_error("invalid dEdx rebin factor");
        // Historical cache has 0.1-keV/cm projections; refined caches may already
        // use the requested width. Never rebin an already-rebinned cache again.
        if(std::fabs(physicalWidth-targetWidth)<1.e-6){}
        else if(std::fabs(physicalWidth-.1)<1.e-6&&s.h->GetNbinsX()%factor==0)s.h->Rebin(factor);
        else throw std::runtime_error("cached projection has unexpected physical bin width");
      }
      if(j)json<<',';json<<"{\"r\":"<<s.r<<",\"fits\":[";
      for(int model=0;model<2;++model){s.fits.push_back(FitOne(s.h,s.r,s.mean,s.sigma,model==1,cfg,index));if(model)json<<',';WriteFit(json,s.fits.back());if(s.fits.back().valid&&(s.chosen<0||s.fits.back().aic<s.fits[s.chosen].aic))s.chosen=model;}
      json<<"],\"chosen\":"<<s.chosen<<"}";canvas.Clear();canvas.Divide(2,1);Slice linear=s;linear.h=static_cast<TH1D*>(s.h->Clone(Form("proto_linear_%d",index)));linear.h->SetDirectory(0);canvas.cd(1);DrawSlice(linear,Required(cfg,"gateNSigma"),false,Required(cfg,"fitNSigmaMin"),Required(cfg,"fitNSigmaMax"));canvas.cd(2);DrawSlice(s,Required(cfg,"gateNSigma"),true,Required(cfg,"fitNSigmaMin"),Required(cfg,"fitNSigmaMax"));canvas.Print(book);delete linear.h;
      std::cout<<"prototype R="<<s.r<<" chosen="<<s.chosen<<std::endl;
      for(unsigned k=0;k<s.fits.size();++k){delete s.fits[k].total;delete s.fits[k].signal;delete s.fits[k].background;}delete s.h;
    }
    json<<"]}\n";json.close();canvas.Print(book+"]");std::cout<<"LOW_R_PROTOTYPE_COMPLETE"<<std::endl;
  }catch(const std::exception& e){std::cerr<<"CheckLowRigidity: "<<e.what()<<std::endl;gSystem->Exit(2);}
#endif
}
