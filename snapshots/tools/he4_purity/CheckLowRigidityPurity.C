// ROOT5 ACLiC tests for shifted He4 peaks, outside-window triton tails and controls.
void CheckLowRigidityPurity(const char* settings,const char* output);
#ifndef __CINT__
#include "FitHe4Purity.C"

namespace {
void LowRequire(bool condition,const char* message) {
  if(!condition)throw std::runtime_error(message);
}
std::vector<double> LowTruth(double rigidity,bool exponential,int control,
                             double calibrationMean,double calibrationSigma) {
  // Synthetic truth; these values are fixture data, not analysis cuts.
  std::vector<double> p(11);
  p[0]=control==1?0.:(control==2?100.:5000.);
  p[1]=-1.2;p[2]=.7;
  p[3]=control==2?20000.:(control==1?1000.:100.);
  p[4]=exponential?(control==2?-.4:-.2):180.;
  const int species[]={1,2};
  for(int neighbor=0;neighbor<2;++neighbor) {
    const NuclearIdDeDxVsMom::SpeciesIndex type=static_cast<NuclearIdDeDxVsMom::SpeciesIndex>(species[neighbor]);
    const int at=5+3*neighbor;
    p[at]=control==1?(neighbor?4000.:5000.):(neighbor?12000.:18000.);
    p[at+1]=(NuclearIdDeDxVsMom::GetMean(type,rigidity)-calibrationMean)/calibrationSigma;
    p[at+2]=NuclearIdDeDxVsMom::GetSigma(type,rigidity)/calibrationSigma;
  }
  return p;
}
double NumericalWindow(const std::vector<double>& p,bool exponential,double fitLo,double fitHi,
                       double gateLo,double gateHi) {
  TF1 signal("low_numeric_signal","[0]*TMath::Gaus(x,[1],[2],0)",-10.,6.);
  signal.SetParameters(p[0],p[1],p[2]);
  TF1 neighbor("low_numeric_neighbor","[0]*TMath::Gaus(x,[1],[2],0)",-10.,6.);
  const TString linear=Form("[0]+([1]-[0])*(x-(%.17g))/(%.17g)",fitLo,fitHi-fitLo);
  TF1 smooth("low_numeric_smooth",exponential?"[0]*exp([1]*x)":linear.Data(),-10.,6.);
  smooth.SetParameters(p[3],p[4]);
  const double s=signal.Integral(gateLo,gateHi);
  double b=smooth.Integral(gateLo,gateHi);
  for(unsigned at=5;at+2<p.size();at+=3) {
    neighbor.SetParameters(p[at],p[at+1],p[at+2]);b+=neighbor.Integral(gateLo,gateHi);
  }
  return s/(s+b);
}
void FillLowToy(TH1D& h,const std::vector<double>& p,bool exponential,double low,double high) {
  h.SetDirectory(0);h.Sumw2();double total=0.;
  for(int bin=1;bin<=h.GetNbinsX();++bin) {
    const double a=h.GetXaxis()->GetBinLowEdge(bin),b=h.GetXaxis()->GetBinUpEdge(bin);
    double value=He4Purity::GaussianIntegral(p,0,a,b)+He4Purity::SmoothIntegral(p,exponential,low,high,a,b);
    for(unsigned at=5;at+2<p.size();at+=3)value+=He4Purity::GaussianIntegral(p,at,a,b);
    value/=h.GetBinWidth(bin);
    LowRequire(value>0. && TMath::Finite(value),"low-rigidity toy expectations are finite and positive");
    h.SetBinContent(bin,value);h.SetBinError(bin,std::sqrt(value));total+=value;
  }
  h.SetEntries(total);
}
void FreeLowFit(He4Purity::Fit& fit) { delete fit.total;delete fit.signal;delete fit.background; }
}
#endif

void CheckLowRigidityPurity(const char* settings,const char* output) {
#ifndef __CINT__
  try {
    const std::map<std::string,std::string> cfg=He4Purity::Settings(settings);
    LowRequire(He4Purity::FeatureEnabled(cfg,"peakSearchEnabled"),"low-rigidity search enabled explicitly");
    const double low=He4Purity::Required(cfg,"fitNSigmaMin"),high=He4Purity::Required(cfg,"fitNSigmaMax");
    const double gate=He4Purity::Required(cfg,"gateNSigma");
    const TString directory(output),jsonPath=directory+"/lowrigidity_closure.json",rootPath=directory+"/lowrigidity_closure.root";
    LowRequire(gSystem->AccessPathName(jsonPath)&&gSystem->AccessPathName(rootPath),"refusing to overwrite earlier low-rigidity tests");
    gSystem->mkdir(directory,kTRUE);
    TFile roots(rootPath,"NEW");LowRequire(!roots.IsZombie(),"low-rigidity toy ROOT output opens");
    std::ofstream json(jsonPath.Data());LowRequire(json.good(),"low-rigidity toy JSON output opens");
    json<<std::setprecision(17)<<"{\"settings\":"<<He4Purity::Quote(settings)<<",\"cases\":[";
    const double rigidities[]={.55,.70};int caseIndex=0;
    for(int rindex=0;rindex<2;++rindex)for(int model=0;model<2;++model) {
      const double r=rigidities[rindex];const bool exponential=model==1;
      const double mean=NuclearIdDeDxVsMom::GetMean(NuclearIdDeDxVsMom::kHe4,r);
      const double sigma=NuclearIdDeDxVsMom::GetSigma(NuclearIdDeDxVsMom::kHe4,r);
      const std::vector<double> truth=LowTruth(r,exponential,0,mean,sigma);
      LowRequire(truth[1]<-He4Purity::Required(cfg,"signalMeanShiftMax"),"injected signal mean lies outside the old allowed bound");
      if(rindex==1)LowRequire(truth[6]<low,"triton center lies outside the fit window in the R=.70 fixture");
      const double expectedOld=He4Purity::Purity(truth,exponential,low,high,gate);
      const double expectedMoving=He4Purity::PurityWindow(truth,exponential,low,high,truth[1]-gate*truth[2],truth[1]+gate*truth[2]);
      LowRequire(std::fabs(expectedOld-NumericalWindow(truth,exponential,low,high,-gate,gate))<1.e-8,"fixed calibration gate analytical/numerical integral closure");
      LowRequire(std::fabs(expectedMoving-NumericalWindow(truth,exponential,low,high,truth[1]-gate*truth[2],truth[1]+gate*truth[2]))<1.e-8,"moving fitted gate analytical/numerical integral closure");
      LowRequire(std::fabs(expectedOld-expectedMoving)>1.e-4,"fixture distinguishes the original and recentered gates");
      TH1D h(Form("low_r%d_%s",rindex,exponential?"exp":"linear"),";original calibrated He4 n#sigma;expected tracks / bin",320,-10.,6.);
      FillLowToy(h,truth,exponential,low,high);
      He4Purity::Fit fit=He4Purity::FitOne(&h,r,mean,sigma,exponential,cfg,12000+caseIndex);
      LowRequire(fit.searchApplied,"search applies to both low-rigidity samples");
      LowRequire(fit.neighborSpecies.size()==2 && fit.neighborSpecies[0]==1 && fit.neighborSpecies[1]==2,
                 "triton and He3 identities retained, including outside-window triton tails");
      LowRequire(fit.valid && fit.quality && fit.status==0 && fit.cov>=2,"shifted-signal toy fit converges with good covariance and residuals");
      LowRequire(fit.signalEvidence && fit.nullValid && fit.nullDeltaAic>He4Purity::Required(cfg,"minimumSignalDeltaAic"),
                 "strong injected He4 is supported over a no-signal model");
      LowRequire(std::fabs(fit.par[1]-truth[1])<.02 && std::fabs(fit.par[2]-truth[2])<.02,"shifted He4 peak location and width recovered");
      LowRequire(std::fabs(fit.rawSignalMean-(mean+sigma*truth[1]))<.02*sigma &&
                 std::fabs(fit.rawSignalSigma-sigma*truth[2])<.02*sigma,"physical dEdx peak mean and sigma recovered");
      LowRequire(std::fabs(fit.purity-expectedOld)<1.e-3,"original calibration-gate purity recovered");
      LowRequire(std::fabs(fit.fittedGatePurity-expectedMoving)<1.e-3,"separate moving-gate purity recovered");
      LowRequire(std::fabs(fit.purity-NumericalWindow(fit.par,exponential,low,high,-gate,gate))<1.e-8,"fitted original-gate analytical/numerical closure");
      LowRequire(std::fabs(fit.fittedGatePurity-NumericalWindow(fit.par,exponential,low,high,fit.par[1]-gate*fit.par[2],fit.par[1]+gate*fit.par[2]))<1.e-8,"fitted moving-gate analytical/numerical closure");
      LowRequire(fit.error>0. && fit.fittedGateError>0. && TMath::Finite(fit.error) && TMath::Finite(fit.fittedGateError),"both gate purities carry finite positive covariance uncertainties");
      for(unsigned j=0;j<fit.neighborSpecies.size();++j)
        LowRequire(fit.par[1]-fit.par[6+3*j]>=He4Purity::Required(cfg,"minimumSignalNeighborSeparation"),"He4 cannot exchange identity with a neighboring nucleus peak");
      if(caseIndex)json<<",";
      json<<"{\"rigidity\":"<<r<<",\"truthOldGatePurity\":"<<expectedOld<<",\"truthMovingGatePurity\":"<<expectedMoving<<",\"fit\":";
      He4Purity::WriteFit(json,fit);json<<"}";
      roots.cd();h.Write();fit.total->Write();fit.signal->Write();fit.background->Write();
      std::cout<<"PASS low-rigidity R="<<r<<" "<<fit.model<<" mean="<<fit.par[1]<<" sigma="<<fit.par[2]
               <<" old/moving purity="<<fit.purity<<"/"<<fit.fittedGatePurity<<" nullDeltaAic="<<fit.nullDeltaAic<<std::endl;
      FreeLowFit(fit);++caseIndex;
    }
    for(int control=1;control<=2;++control) {
      const double r=.70,mean=NuclearIdDeDxVsMom::GetMean(NuclearIdDeDxVsMom::kHe4,r),sigma=NuclearIdDeDxVsMom::GetSigma(NuclearIdDeDxVsMom::kHe4,r);
      const std::vector<double> truth=LowTruth(r,true,control,mean,sigma);
      TH1D h(Form("low_control_%d",control),";original calibrated He4 n#sigma;expected tracks / bin",320,-10.,6.);
      FillLowToy(h,truth,true,low,high);
      He4Purity::Fit fit=He4Purity::FitOne(&h,r,mean,sigma,true,cfg,12000+caseIndex);
      if(control==1)LowRequire(!fit.signalEvidence,"signal-free control must not report He4 evidence");
      else {
        for(int bin=1;bin<h.GetNbinsX();++bin) {
          const double u=h.GetBinCenter(bin);
          if(u>=He4Purity::Required(cfg,"signalSearchMeanMin") && u<=He4Purity::Required(cfg,"signalSearchMeanMax"))
            LowRequire(h.GetBinContent(bin)>h.GetBinContent(bin+1),"weak-shoulder control is strictly falling across the search region");
        }
        LowRequire(!fit.rawLocalMaximum,"a weak Gaussian shoulder on monotonic background is not labeled a raw local maximum");
      }
      json<<",{\"control\":"<<He4Purity::Quote(control==1?"no_signal":"monotonic_weak_shoulder")<<",\"fit\":";
      He4Purity::WriteFit(json,fit);json<<"}";
      roots.cd();h.Write();fit.total->Write();fit.signal->Write();fit.background->Write();
      std::cout<<"PASS control="<<control<<" rawLocalMaximum="<<fit.rawLocalMaximum<<" signalEvidence="<<fit.signalEvidence
               <<" nullDeltaAic="<<fit.nullDeltaAic<<std::endl;
      FreeLowFit(fit);++caseIndex;
    }
    json<<"],\"status\":\"passed\"}\n";json.close();roots.Close();
    std::cout<<"PASS all 6 low-rigidity shifted-peak and control cases"<<std::endl;
  } catch(const std::exception& error) {
    std::cerr<<"FAIL low-rigidity He4 purity tests: "<<error.what()<<std::endl;gSystem->Exit(1);
  }
#endif
}
