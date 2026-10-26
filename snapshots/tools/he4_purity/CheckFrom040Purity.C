// ROOT5 Asimov tests for the first 0.40-0.41 rigidity bin and saved-axis clipping.
void CheckFrom040Purity(const char* settings,const char* output);
#ifndef __CINT__
#include "FitHe4Purity.C"

namespace {
void From040Require(bool ok,const char* message){if(!ok)throw std::runtime_error(message);}
std::vector<double> From040Truth(double r,double mean,double sigma,bool exponential,double shift,double width){
  std::vector<double> p(11);
  // Generate a genuinely separate raw maximum after the saved .5 keV/cm binning
  // and diagnostic smoothing. A weaker component can instead be just a shoulder.
  p[0]=10000.;p[1]=shift;p[2]=width;p[3]=100.;p[4]=exponential?-.2:180.;
  const int species[]={1,2};
  for(int j=0;j<2;++j){
    const NuclearIdDeDxVsMom::SpeciesIndex type=static_cast<NuclearIdDeDxVsMom::SpeciesIndex>(species[j]);
    p[5+3*j]=j?14000.:20000.;
    p[6+3*j]=(NuclearIdDeDxVsMom::GetMean(type,r)-mean)/sigma;
    p[7+3*j]=NuclearIdDeDxVsMom::GetSigma(type,r)/sigma;
  }
  return p;
}
void From040Fill(TH1D& h,const std::vector<double>& p,bool exponential,double low,double high){
  h.SetDirectory(0);h.Sumw2();double count=0.;
  for(int bin=1;bin<=h.GetNbinsX();++bin){
    const double a=h.GetXaxis()->GetBinLowEdge(bin),b=h.GetXaxis()->GetBinUpEdge(bin);
    double n=He4Purity::GaussianIntegral(p,0,a,b)+He4Purity::SmoothIntegral(p,exponential,low,high,a,b);
    for(unsigned j=5;j+2<p.size();j+=3)n+=He4Purity::GaussianIntegral(p,j,a,b);
    n/=h.GetBinWidth(bin);From040Require(n>0.&&TMath::Finite(n),"finite positive axis-clipped bin expectations");
    h.SetBinContent(bin,n);h.SetBinError(bin,std::sqrt(n));count+=n;
  }
  h.SetEntries(count);
}
double From040Numeric(const std::vector<double>& p,bool exponential,double low,double high,double a,double b){
  TF1 signal("from040_signal_integral","[0]*TMath::Gaus(x,[1],[2],0)",-20.,8.);
  signal.SetParameters(p[0],p[1],p[2]);
  const TString linear=Form("[0]+([1]-[0])*(x-(%.17g))/(%.17g)",low,high-low);
  TF1 smooth("from040_bg_integral",exponential?"[0]*exp([1]*x)":linear.Data(),-20.,8.);
  smooth.SetParameters(p[3],p[4]);
  TF1 nucleus("from040_neighbor_integral","[0]*TMath::Gaus(x,[1],[2],0)",-20.,8.);
  const double s=signal.Integral(a,b);double bg=smooth.Integral(a,b);
  for(unsigned j=5;j+2<p.size();j+=3){nucleus.SetParameters(p[j],p[j+1],p[j+2]);bg+=nucleus.Integral(a,b);}
  return s/(s+bg);
}
void From040Free(He4Purity::Fit& fit){delete fit.total;delete fit.signal;delete fit.background;}
}
#endif
void CheckFrom040Purity(const char* settings,const char* output){
#ifndef __CINT__
  try{
    const std::map<std::string,std::string> cfg=He4Purity::Settings(settings);
    From040Require(He4Purity::FeatureEnabled(cfg,"peakSearchEnabled")&&He4Purity::Required(cfg,"peakSearchMinRigidity")<=.4,
                   "explicit He4 search begins at 0.40 rigidity");
    const double gate=He4Purity::Required(cfg,"gateNSigma"),requestedLow=He4Purity::Required(cfg,"fitNSigmaMin"),requestedHigh=He4Purity::Required(cfg,"fitNSigmaMax");
    const TString folder(output),jsonPath=folder+"/from040_closure.json",rootPath=folder+"/from040_closure.root";
    From040Require(gSystem->AccessPathName(jsonPath)&&gSystem->AccessPathName(rootPath),"preserve earlier from040 tests");
    gSystem->mkdir(folder,kTRUE);TFile roots(rootPath,"NEW");From040Require(!roots.IsZombie(),"new from040 ROOT output opens");
    std::ofstream json(jsonPath.Data());From040Require(json.good(),"new from040 JSON opens");json<<std::setprecision(17)<<"{\"settings\":"<<He4Purity::Quote(settings)<<",\"physicalAxis\":[0,120],\"physicalBinWidth\":0.5,\"cases\":[";
    int sequence=0;
    for(int ri=0;ri<2;++ri)for(int model=0;model<2;++model){
      const double r=ri?.425:.405;const bool exponential=model==1;
      const double mean=NuclearIdDeDxVsMom::GetMean(NuclearIdDeDxVsMom::kHe4,r),sigma=NuclearIdDeDxVsMom::GetSigma(NuclearIdDeDxVsMom::kHe4,r);
      const double axisLo=-mean/sigma,axisHi=(120.-mean)/sigma,low=std::max(requestedLow,axisLo),high=std::min(requestedHigh,axisHi);
      From040Require(axisHi<requestedHigh && high>gate,"fit upper edge clipped, but fixed calibration gate fully stored");
      // Exercise substantially shifted first-bin candidates, close to the newly
      // explored region while still separated from the two calibrated BG peaks.
      const std::vector<double> truth=From040Truth(r,mean,sigma,exponential,ri?-2.8:-3.0,.7);
      TH1D h(Form("from040_r%d_%d",ri,model),";original calibrated He4 n#sigma;expected tracks / saved 0.5 keV/cm bin",240,axisLo,axisHi);
      From040Fill(h,truth,exponential,low,high);
      He4Purity::Fit fit=He4Purity::FitOne(&h,r,mean,sigma,exponential,cfg,14000+sequence);
      From040Require(fit.searchApplied&&fit.lowerRigiditySearchApplied&&fit.neighborSpecies.size()==2&&fit.neighborSpecies[0]==1&&fit.neighborSpecies[1]==2,
                     "first bins explicitly fit separate triton and He3 background identities");
      From040Require(fit.rawLocalMaximum,"strong shifted first-bin signals are recognized inside the lower-rigidity diagnostic window");
      From040Require(fit.total&&fit.total->GetXmax()==high,"TF1 uses the supported saved-axis upper limit");
      From040Require(fit.valid&&fit.quality&&fit.status==0&&fit.cov>=2&&fit.signalEvidence,"clipped-axis shifted signal fits successfully");
      From040Require(std::fabs(fit.par[1]-truth[1])<.02&&std::fabs(fit.par[2]-truth[2])<.02,"first-bin signal shift and width recovered");
      for(int neighbor=0;neighbor<2;++neighbor)From040Require(fit.par[1]-fit.par[6+3*neighbor]>=He4Purity::Required(cfg,"minimumSignalNeighborSeparation"),"He4 ordering remains separate from both lower background peaks");
      const double oldTruth=He4Purity::Purity(truth,exponential,low,high,gate),movingTruth=He4Purity::PurityWindow(truth,exponential,low,high,truth[1]-gate*truth[2],truth[1]+gate*truth[2]);
      From040Require(std::fabs(oldTruth-From040Numeric(truth,exponential,low,high,-gate,gate))<1.e-8,"axis-clipped linear/exponential integral normalization closes numerically");
      From040Require(std::fabs(movingTruth-From040Numeric(truth,exponential,low,high,truth[1]-gate*truth[2],truth[1]+gate*truth[2]))<1.e-8,"axis-clipped moving-gate truth closes numerically");
      From040Require(std::fabs(fit.purity-oldTruth)<1.e-3,"original fixed gate purity recovered in first rigidity bins");
      From040Require(fit.fittedGateWithinFitRange&&std::fabs(fit.fittedGatePurity-movingTruth)<1.e-3,"supported shifted gate remains distinct and recovers truth");
      From040Require(std::fabs(fit.purity-From040Numeric(fit.par,exponential,low,high,-gate,gate))<1.e-8,"fitted fixed-gate analytical/numerical closure");
      From040Require(std::fabs(fit.fittedGatePurity-From040Numeric(fit.par,exponential,low,high,fit.par[1]-gate*fit.par[2],fit.par[1]+gate*fit.par[2]))<1.e-8,"fitted moving-gate analytical/numerical closure");
      if(sequence)json<<",";
      json<<"{\"rigidity\":"<<r<<",\"axisUpperNSigma\":"<<axisHi<<",\"effectiveFitUpper\":"<<high<<",\"truthOriginalGatePurity\":"<<oldTruth<<",\"truthMovingGatePurity\":"<<movingTruth<<",\"fit\":";
      He4Purity::WriteFit(json,fit);json<<"}";
      roots.cd();h.Write();fit.total->Write();fit.signal->Write();fit.background->Write();
      std::cout<<"PASS from040 R="<<r<<" "<<fit.model<<" fitHigh="<<high<<" mean="<<fit.par[1]<<" old/moving="<<fit.purity<<"/"<<fit.fittedGatePurity<<std::endl;
      From040Free(fit);++sequence;
    }
    {
      // Exercise the new lower-mean interval itself, beyond the previous -3.2
      // bound, and a raw peak that the old diagnostic window would have missed.
      const double r=.405,mean=NuclearIdDeDxVsMom::GetMean(NuclearIdDeDxVsMom::kHe4,r),sigma=NuclearIdDeDxVsMom::GetSigma(NuclearIdDeDxVsMom::kHe4,r);
      const double axisLo=-mean/sigma,axisHi=(120.-mean)/sigma,low=std::max(requestedLow,axisLo),high=std::min(requestedHigh,axisHi);
      std::vector<double> truth=From040Truth(r,mean,sigma,true,-3.5,.4);truth[0]=20000.;
      TH1D h("from040_beyond_previous_lower_bound",";calibrated He4 n#sigma;expected counts",240,axisLo,axisHi);From040Fill(h,truth,true,low,high);
      He4Purity::Fit fit=He4Purity::FitOne(&h,r,mean,sigma,true,cfg,14000+sequence);
      From040Require(fit.lowerRigiditySearchApplied&&fit.valid&&fit.quality&&fit.signalEvidence,"new lower-bound interval supports a well-determined signal");
      From040Require(fit.par[1]<-3.2&&std::fabs(fit.par[1]-truth[1])<.02&&std::fabs(fit.par[2]-truth[2])<.02,"He4 peak beyond the previous mean bound is recovered");
      From040Require(fit.rawLocalMaximum&&fit.rawPeakMean<-3.2,"full lower-rigidity diagnostic window locates the out-of-old-range raw peak");
      From040Require(fit.neighborSpecies.size()==2&&fit.neighborSpecies[0]==1&&fit.neighborSpecies[1]==2,"extended signal is still separated from triton and He3 identities");
      const double oldTruth=He4Purity::Purity(truth,true,low,high,gate),movingTruth=He4Purity::PurityWindow(truth,true,low,high,truth[1]-gate*truth[2],truth[1]+gate*truth[2]);
      From040Require(std::fabs(fit.purity-oldTruth)<1.e-3&&fit.fittedGateWithinFitRange&&std::fabs(fit.fittedGatePurity-movingTruth)<1.e-3,"fixed and moving purities recover truth beyond old mean bounds");
      From040Require(std::fabs(oldTruth-From040Numeric(truth,true,low,high,-gate,gate))<1.e-8&&
                     std::fabs(movingTruth-From040Numeric(truth,true,low,high,truth[1]-gate*truth[2],truth[1]+gate*truth[2]))<1.e-8,"independent fixed/moving integrals close in the additional lower-bound case");
      json<<",{\"control\":\"signal_beyond_previous_minus3p2_limit\",\"truthOriginalGatePurity\":"<<oldTruth<<",\"truthMovingGatePurity\":"<<movingTruth<<",\"fit\":";He4Purity::WriteFit(json,fit);json<<"}";
      roots.cd();h.Write();fit.total->Write();fit.signal->Write();fit.background->Write();
      std::cout<<"PASS extended lower mean="<<fit.par[1]<<" rawPeak="<<fit.rawPeakMean<<" old/moving="<<fit.purity<<"/"<<fit.fittedGatePurity<<std::endl;
      From040Free(fit);++sequence;
    }
    {
      const double r=.405,mean=NuclearIdDeDxVsMom::GetMean(NuclearIdDeDxVsMom::kHe4,r),sigma=NuclearIdDeDxVsMom::GetSigma(NuclearIdDeDxVsMom::kHe4,r);
      const double axisLo=-mean/sigma,axisHi=(120.-mean)/sigma,low=std::max(requestedLow,axisLo),high=std::min(requestedHigh,axisHi);
      const std::vector<double> truth=From040Truth(r,mean,sigma,true,.8,1.);
      TH1D h("from040_unsupported_moving_gate",";calibrated He4 n#sigma;expected counts",240,axisLo,axisHi);From040Fill(h,truth,true,low,high);
      He4Purity::Fit fit=He4Purity::FitOne(&h,r,mean,sigma,true,cfg,14000+sequence);
      From040Require(fit.valid&&fit.status==0&&fit.cov>=2,"upper-shifted peak still determines a valid original-gate fit");
      From040Require(!fit.fittedGateWithinFitRange&&fit.fittedGatePurity==-1.&&fit.fittedGateError==0.,"unsupported moving gate is flagged and never extrapolated as measured purity");
      From040Require(std::fabs(fit.purity-He4Purity::Purity(truth,true,low,high,gate))<1.e-3,"fixed old gate remains supported when moving gate is clipped");
      json<<",{\"control\":\"moving_gate_beyond_saved_axis\",\"fit\":";He4Purity::WriteFit(json,fit);json<<"}";
      roots.cd();h.Write();fit.total->Write();From040Free(fit);++sequence;
    }
    {
      const double r=.4,mean=NuclearIdDeDxVsMom::GetMean(NuclearIdDeDxVsMom::kHe4,r),sigma=NuclearIdDeDxVsMom::GetSigma(NuclearIdDeDxVsMom::kHe4,r),axisLo=-mean/sigma,axisHi=(120.-mean)/sigma;
      From040Require(axisHi<gate,"exact R=.400 endpoint lacks part of fixed +/-2 gate");
      TH1D h("from040_exact_endpoint",";calibrated He4 n#sigma;expected counts",240,axisLo,axisHi);
      const std::vector<double> truth=From040Truth(r,mean,sigma,true,-1.6,.7);From040Fill(h,truth,true,requestedLow,axisHi);
      He4Purity::Fit fit=He4Purity::FitOne(&h,r,mean,sigma,true,cfg,14000+sequence);
      From040Require(!fit.valid&&!fit.quality&&fit.total==0&&fit.purity==-1.,"exact endpoint with unsupported fixed gate is rejected");
      json<<",{\"control\":\"exact_0400_endpoint_not_first_bin_center\",\"axisUpperNSigma\":"<<axisHi<<",\"valid\":false}";
      roots.cd();h.Write();++sequence;
    }
    json<<"],\"status\":\"passed\"}\n";json.close();roots.Close();
    std::cout<<"PASS all "<<sequence<<" from040 shifted-signal and saved-axis controls"<<std::endl;
  }catch(const std::exception& error){std::cerr<<"FAIL from040 tests: "<<error.what()<<std::endl;gSystem->Exit(1);}
#endif
}
