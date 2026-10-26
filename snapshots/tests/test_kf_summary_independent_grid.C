// ROOT5 toy-only check: unequal CF/purity bin boundaries and display-independent fits.
#include "../tools/kfparticle_summary/KFSummarySupport.h"
#include "TH2D.h"
#include "TCanvas.h"
#include <stdexcept>

void test_kf_summary_independent_grid(const char* auditPath) {
  try {
    const double cfg[9]={1.08,1.15,1.11,1.12,.001,.008,.003,10,2};
    KFReport::Init(auditPath,1.11596,.00296,cfg);
    TH2D source("toy_original","toy",100,0,.5,200,1.05,1.25);source.Sumw2();
    double entries=0;
    for (int x=1;x<=100;++x) for (int y=1;y<=200;++y) {
      const int cell=(x-1)/10;
      const double mass=source.GetYaxis()->GetBinCenter(y);
      const double value=1000*TMath::Exp(-.5*TMath::Power((mass-1.11596)/.00296,2))+20+10*cell;
      source.SetBinContent(x,y,value);source.SetBinError(x,y,TMath::Sqrt(value));entries+=value;
    }
    source.SetEntries(entries);
    TH1D signal("toy_signal","toy",25,0,.5),purity("toy_purity","toy",25,0,.5);
    const double minM=1.11596-3*.00296,maxM=1.11596+3*.00296;
    KFReport::Method5IndependentGrid(&source,&signal,&purity,4,10,2,minM,maxM);
    const int lo=source.GetYaxis()->FindBin(minM),hi=source.GetYaxis()->FindBin(maxM);
    double e0=0,e1=0;
    const double n0=source.IntegralAndError(9,10,lo,hi,e0);
    const double n1=source.IntegralAndError(11,12,lo,hi,e1);
    const double p0=purity.GetBinContent(1),p1=purity.GetBinContent(4);
    const double ep0=purity.GetBinError(1),ep1=purity.GetBinError(4);
    if (!(p0>0 && p1>0 && ep0>0 && ep1>0))throw std::runtime_error("toy purity fits unavailable");
    const double expected=n0*p0+n1*p1;
    const double variance=p0*p0*e0*e0+p1*p1*e1*e1+n0*n0*ep0*ep0+n1*n1*ep1*ep1;
    if (TMath::Abs(signal.GetBinContent(3)-expected)>1e-10*expected)
      throw std::runtime_error("straddling output bin did not retain both purity cells");
    if (TMath::Abs(TMath::Power(signal.GetBinError(3),2)-variance)>1e-10*variance)
      throw std::runtime_error("shared-purity error within an output bin is wrong");
    if (TMath::Abs(signal.GetBinContent(3)-(n0+n1)*p1)<1e-6*expected)
      throw std::runtime_error("toy did not distinguish center assignment from fine-bin weighting");

    TH1D range("toy_range","toy",20,0,.4);
    for (int b=1;b<=20;++b) {range.SetBinContent(b,b<=15?1:5);range.SetBinError(b,.1);}
    range.GetXaxis()->SetRangeUser(0,.3);
    const int first=range.GetXaxis()->GetFirst(),last=range.GetXaxis()->GetLast();
    TCanvas canvas("toy_canvas","toy",400,300);range.Draw();
    TF1 constant("toy_range_fit","pol0",0,.4);
    KFReport::FitLL(&range,&constant,"R");
    if (TMath::Abs(constant.GetParameter(0)-2)>1e-6)
      throw std::runtime_error("display range truncated the independent fit interval");
    if (range.GetXaxis()->GetFirst()!=first || range.GetXaxis()->GetLast()!=last)
      throw std::runtime_error("display range was not restored after fitting");
    std::cout<<"INDEPENDENT_GRID_AND_FIT_RANGE_PASS"<<std::endl;
  } catch (const std::exception& e) {
    std::cerr<<"INDEPENDENT_GRID_TEST_FAILED: "<<e.what()<<std::endl;gSystem->Exit(1);
  }
}
