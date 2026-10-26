// Deterministic closure test, run with ROOT5 ACLiC (.C+).
void CheckToyPurity(const char* settings,const char* output);
#ifndef __CINT__
#include "FitHe4Purity.C"
#include <limits>

namespace {
void ToyRequire(bool condition,const char* explanation) {
  if(!condition)throw std::runtime_error(explanation);
}
double NumericalTruth(const std::vector<double>& parameters,bool exponential,
                      double fitLo,double fitHi,double gate) {
  // Use independently integrated ROOT formulas, not the analytical helpers.
  TF1 signal("toy_integral_signal","[0]*TMath::Gaus(x,[1],[2],0)",-10.,6.);
  signal.SetParameters(parameters[0],parameters[1],parameters[2]);
  TF1 neighbor("toy_integral_neighbor","[0]*TMath::Gaus(x,[1],[2],0)",-10.,6.);
  neighbor.SetParameters(parameters[5],parameters[6],parameters[7]);
  const TString linear=Form("[0]+([1]-[0])*(x-(%.17g))/(%.17g)",fitLo,fitHi-fitLo);
  TF1 background("toy_integral_background",exponential?"[0]*exp([1]*x)":linear.Data(),-10.,6.);
  background.SetParameters(parameters[3],parameters[4]);
  const double s=signal.Integral(-gate,gate);
  const double b=neighbor.Integral(-gate,gate)+background.Integral(-gate,gate);
  return s/(s+b);
}
}
#endif

void CheckToyPurity(const char* settings,const char* output) {
#ifndef __CINT__
  try {
    const std::map<std::string,std::string> cfg=He4Purity::Settings(settings);
    const double low=He4Purity::Required(cfg,"fitNSigmaMin");
    const double high=He4Purity::Required(cfg,"fitNSigmaMax");
    const double gate=He4Purity::Required(cfg,"gateNSigma");
    ToyRequire(low==-7. && high==4. && gate==2.,"fixture expects published [-7,4] fit range and fixed +/-2 calibration gate");
    const TString directory(output),jsonPath=directory+"/toy_closure.json",rootPath=directory+"/toy_closure.root";
    ToyRequire(gSystem->AccessPathName(jsonPath) && gSystem->AccessPathName(rootPath),"refusing to overwrite previous toy results");
    gSystem->mkdir(directory,kTRUE);
    const double rigidity=1.5;
    const double calMean=NuclearIdDeDxVsMom::GetMean(NuclearIdDeDxVsMom::kHe4,rigidity);
    const double calSigma=NuclearIdDeDxVsMom::GetSigma(NuclearIdDeDxVsMom::kHe4,rigidity);
    const double neighborMean=(NuclearIdDeDxVsMom::GetMean(NuclearIdDeDxVsMom::kHe3,rigidity)-calMean)/calSigma;
    const double neighborSigma=NuclearIdDeDxVsMom::GetSigma(NuclearIdDeDxVsMom::kHe3,rigidity)/calSigma;
    ToyRequire(neighborMean>low && neighborMean<high && neighborSigma>0.,"He3 nuisance peak lies inside the toy fit range");
    for(int species=0;species<2;++species) {
      const double mu=(NuclearIdDeDxVsMom::GetMean(static_cast<NuclearIdDeDxVsMom::SpeciesIndex>(species),rigidity)-calMean)/calSigma;
      ToyRequire(mu<low || mu>high,"only He3 contributes a neighboring nucleus peak in this fixture");
    }
    TFile resultFile(rootPath,"NEW");ToyRequire(!resultFile.IsZombie(),"toy ROOT output opens");
    std::ofstream report(jsonPath.Data());ToyRequire(report.good(),"toy JSON output opens");
    report<<std::setprecision(17)<<"{\"rigidity\":"<<rigidity<<",\"settings\":"<<He4Purity::Quote(settings)
          <<",\"binWidth\":0.05,\"fitRange\":["<<low<<","<<high<<"],\"gateNSigma\":"<<gate<<",\"cases\":[";
    for(int model=0;model<2;++model) {
      const bool exponential=model==1;
      std::vector<double> truth(8);
      truth[0]=5000.;truth[1]=.05;truth[2]=.8;
      truth[3]=10.;truth[4]=exponential?-.2:30.;
      truth[5]=12000.;truth[6]=neighborMean;truth[7]=neighborSigma;
      const double expected=He4Purity::Purity(truth,exponential,low,high,gate);
      const double independent=NumericalTruth(truth,exponential,low,high,gate);
      ToyRequire(std::fabs(expected-independent)<1.e-8,"analytic Gaussian/background gate integrals match independent ROOT numerical integration");
      TH1D histogram(Form("toy_%s",exponential?"exponential":"linear"),";calibrated n#sigma_{^{4}He};expected tracks / bin",320,-10.,6.);
      histogram.SetDirectory(0);histogram.Sumw2();
      double count=0.;
      for(int bin=1;bin<=histogram.GetNbinsX();++bin) {
        const double a=histogram.GetXaxis()->GetBinLowEdge(bin),b=histogram.GetXaxis()->GetBinUpEdge(bin);
        const double expectation=(He4Purity::GaussianIntegral(truth,0,a,b)+
          He4Purity::GaussianIntegral(truth,5,a,b)+He4Purity::SmoothIntegral(truth,exponential,low,high,a,b))/histogram.GetBinWidth(bin);
        ToyRequire(expectation>0. && TMath::Finite(expectation),"toy bin expectation is finite and positive");
        histogram.SetBinContent(bin,expectation);histogram.SetBinError(bin,std::sqrt(expectation));count+=expectation;
      }
      histogram.SetEntries(count);
      He4Purity::Fit fit=He4Purity::FitOne(&histogram,rigidity,calMean,calSigma,exponential,cfg,9000+model);
      const double numericalFitted=NumericalTruth(fit.par,exponential,low,high,gate);
      ToyRequire(std::fabs(fit.purity-numericalFitted)<1.e-8,"fitted analytical purity agrees with independent numerical integration");
      ToyRequire(fit.status==0 && fit.cov>=2 && fit.valid && fit.quality,"toy fit has valid minimizer, covariance, model residuals and precision");
      ToyRequire(TMath::Finite(fit.error) && fit.error>0.,"toy purity covariance gives a positive finite uncertainty");
      ToyRequire(std::fabs(fit.purity-expected)<1.e-3,"recovered purity agrees with known truth to better than 0.001 absolute");
      ToyRequire(std::fabs(fit.par[1]-.05)<.01 && std::fabs(fit.par[2]-.8)<.01,"toy recovers the injected He4 mean and width");
      resultFile.cd();histogram.Write();fit.total->Write();fit.signal->Write();fit.background->Write();
      if(model)report<<",";
      report<<"{\"truthPurity\":"<<expected<<",\"numericalTruthPurity\":"<<independent<<",\"absoluteBias\":"<<std::fabs(fit.purity-expected)<<",\"fit\":";
      He4Purity::WriteFit(report,fit);report<<"}";
      std::cout<<"PASS toy "<<fit.model<<" truth="<<expected<<" fit="<<fit.purity<<" +/- "<<fit.error
               <<" status="<<fit.status<<" covariance="<<fit.cov<<std::endl;
      delete fit.total;delete fit.signal;delete fit.background;
    }
    report<<"],\"status\":\"passed\"}\n";report.close();resultFile.Close();
    std::cout<<"PASS both He4 purity toy closure cases"<<std::endl;
  } catch(const std::exception& error) {
    std::cerr<<"FAIL He4 purity toy closure: "<<error.what()<<std::endl;gSystem->Exit(1);
  }
#endif
}
