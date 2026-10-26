// INPUT ONLY: READ fixed source once; preserve native stored counts, no fit/PID.
#include "TFile.h"
#include "TList.h"
#include "TH2.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TGaxis.h"
#include "TLatex.h"
#include "TSystem.h"
#include "TPad.h"
#include "TMath.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>

#ifndef __CINT__
namespace NuclearTofInputs {
typedef std::map<std::string,double> Settings;
const char* inputKeys[]={"true/hDedxP","true/hDedxP_d_m2","true/hDedxP_t_m2","true/hDedxP_3He_m2"};
const char* names[]={"ALL","TAG_D","TAG_T","TAG_HE3","REMAINING"};
const char* labels[]={"ALL pre-PID","TOF D / He4-compatible (not pure isotope)","TOF T window (not pure isotope)","TOF 3He window (not pure isotope)","ALL minus three disjoint TOF tags"};
void Require(bool ok,const std::string& message){if(!ok)throw std::runtime_error(message);}
double Need(const Settings& c,const std::string& key){Settings::const_iterator i=c.find(key);Require(i!=c.end()&&TMath::Finite(i->second),"Missing/nonfinite explicit field: "+key);return i->second;}
int Integer(const Settings& c,const std::string& key){double x=Need(c,key);Require(x==std::floor(x),"Noninteger field: "+key);return int(x);}
bool Near(double a,double b,double tolerance){return std::fabs(a-b)<=tolerance*std::max(1.,std::max(std::fabs(a),std::fabs(b)));}
Settings Read(const char* path){std::ifstream f(path);Require(bool(f),"Cannot read explicit settings");Settings c;std::string key;double value;while(f>>key>>value)Require(c.insert(std::make_pair(key,value)).second,"Duplicate setting");Require(f.eof(),"Malformed settings");return c;}
void Axis(TAxis* actual,TAxis* other){Require(actual->GetNbins()==other->GetNbins(),"Native axis dimensions differ");for(int b=1;b<=actual->GetNbins()+1;++b)Require(actual->GetBinLowEdge(b)==other->GetBinLowEdge(b),"Native axes are not exactly equal");}
void Validate(const Settings& c){
  Require(Integer(c,"schemaVersion")==1,"Unsupported INPUT schema");
  Require(Integer(c,"expectedXBins")==600&&Integer(c,"expectedYBins")==1200&&Need(c,"expectedXMin")==0&&Need(c,"expectedXMax")==6&&Need(c,"expectedYMin")==0&&Need(c,"expectedYMax")==120,"Unsupported authorized native dimensions");
  Require(Need(c,"expectedXWidth")==.01&&Need(c,"expectedYWidth")==.1,"Unsupported native widths");
  Require(Integer(c,"displayDedxRebinFactor")>0&&Integer(c,"expectedYBins")%Integer(c,"displayDedxRebinFactor")==0,"Invalid display-only rebin");
  Require(Need(c,"axisRelativeTolerance")>0&&Need(c,"axisRelativeTolerance")<=1.e-6&&Need(c,"varianceRelativeTolerance")>0&&Need(c,"varianceRelativeTolerance")<=1.e-6&&Integer(c,"axisMaxDigits")>0&&Integer(c,"modeSmoothingHalfBins")>=0&&Need(c,"modeMomentHalfWindow")>0,"Invalid audit/diagnostic settings");
  Require(Need(c,"fullDedxMin")==0&&Need(c,"fullDedxMax")==120&&Need(c,"supplementalDedxMin")>=0&&Need(c,"supplementalDedxMin")<Need(c,"supplementalDedxMax")&&Need(c,"supplementalDedxMax")<=120,"Invalid display ranges");
  Require(Need(c,"overviewRigidityMin")==0&&Need(c,"overviewRigidityMax")==6,"Overview must preserve entire native R domain");
  Require(Need(c,"atlasRigidityMin")>=0&&Need(c,"atlasRigidityMax")<=6&&Need(c,"atlasRigidityMax")>Need(c,"atlasRigidityMin")&&Need(c,"coarseRigidityStep")>=Need(c,"expectedXWidth"),"Invalid non-overlapping R partition");
  const char* aligned[]={"atlasRigidityMin","atlasRigidityMax","coarseRigidityStep","diagnosticRigidityMin","diagnosticRigidityMax"};
  for(unsigned i=0;i<5;++i){double q=Need(c,aligned[i])/Need(c,"expectedXWidth");Require(Near(q,std::floor(q+.5),Need(c,"axisRelativeTolerance")),"R edge/step does not align with native bins");}
  double count=(Need(c,"atlasRigidityMax")-Need(c,"atlasRigidityMin"))/Need(c,"coarseRigidityStep");Require(Near(count,std::floor(count+.5),Need(c,"axisRelativeTolerance")),"Coarse partition has incomplete last group");
  Require(Need(c,"diagnosticRigidityMin")>=Need(c,"atlasRigidityMin")&&Need(c,"diagnosticRigidityMax")<=Need(c,"atlasRigidityMax")&&Need(c,"diagnosticRigidityMax")>Need(c,"diagnosticRigidityMin"),"Diagnostic interval outside atlas");
  Require(Need(c,"logMinimum")>0&&Need(c,"logMaximumFactor")>1&&Need(c,"linearMaximumFactor")>1,"Invalid y display settings");
  for(int kind=0;kind<2;++kind)for(int side=0;side<4;++side){const char* sides[]={"Left","Right","Bottom","Top"};std::string key=std::string(kind?"projectionMargin":"mapMargin")+sides[side];Require(Need(c,key)>0&&Need(c,key)<.5,"Invalid plot margin");}
}
struct Audit{double count[5],variance[5],entries[5],regular[5],minimumRemaining,maximumVarianceRelative;unsigned cells;Audit():minimumRemaining(1.e300),maximumVarianceRelative(0),cells(0){for(int s=0;s<5;++s)count[s]=variance[s]=entries[s]=regular[s]=0;}};
void Native(TH2* source[4],TH2D* cache[5],const Settings& c,Audit& a){
  for(int s=0;s<4;++s){Require(std::string(source[s]->ClassName())=="TH2F"&&source[s]->GetSumw2N()==0,"Expected original TH2F with absent Sumw2");Require(source[s]->GetNbinsX()==Integer(c,"expectedXBins")&&source[s]->GetNbinsY()==Integer(c,"expectedYBins"),"Unexpected native shape");Axis(source[0]->GetXaxis(),source[s]->GetXaxis());Axis(source[0]->GetYaxis(),source[s]->GetYaxis());a.entries[s]=source[s]->GetEntries();}
  Require(source[0]->GetXaxis()->GetXmin()==Need(c,"expectedXMin")&&source[0]->GetXaxis()->GetXmax()==Need(c,"expectedXMax")&&source[0]->GetYaxis()->GetXmin()==Need(c,"expectedYMin")&&source[0]->GetYaxis()->GetXmax()==Need(c,"expectedYMax"),"Native endpoints differ");
  for(int b=1;b<=source[0]->GetNbinsX();++b)Require(Near(source[0]->GetXaxis()->GetBinWidth(b),Need(c,"expectedXWidth"),Need(c,"axisRelativeTolerance")),"Native R width differs");
  for(int b=1;b<=source[0]->GetNbinsY();++b)Require(Near(source[0]->GetYaxis()->GetBinWidth(b),Need(c,"expectedYWidth"),Need(c,"axisRelativeTolerance")),"Native dE/dx width differs");
  for(int s=0;s<5;++s){cache[s]=new TH2D(names[s],";R = |p/q| (GeV/c);TPC dE/dx (keV/cm)",source[0]->GetNbinsX(),source[0]->GetXaxis()->GetXmin(),source[0]->GetXaxis()->GetXmax(),source[0]->GetNbinsY(),source[0]->GetYaxis()->GetXmin(),source[0]->GetYaxis()->GetXmax());cache[s]->SetDirectory(0);cache[s]->Sumw2();Axis(source[0]->GetXaxis(),cache[s]->GetXaxis());Axis(source[0]->GetYaxis(),cache[s]->GetYaxis());}
  for(int x=0;x<=source[0]->GetNbinsX()+1;++x)for(int y=0;y<=source[0]->GetNbinsY()+1;++y){
    double n[5];for(int s=0;s<4;++s){n[s]=source[s]->GetBinContent(x,y);double e=source[s]->GetBinError(x,y),v=e*e;Require(TMath::Finite(n[s])&&n[s]>=0&&TMath::Finite(e)&&Near(v,n[s],Need(c,"varianceRelativeTolerance")),"Source count/Poisson-type error invalid");a.maximumVarianceRelative=std::max(a.maximumVarianceRelative,std::fabs(v-n[s])/std::max(1.,n[s]));}
    n[4]=n[0]-n[1]-n[2]-n[3];Require(n[4]>=0,"Disjoint TOF tag sum exceeds ALL bin");a.minimumRemaining=std::min(a.minimumRemaining,n[4]);
    for(int s=0;s<5;++s){cache[s]->SetBinContent(x,y,n[s]);cache[s]->SetBinError(x,y,std::sqrt(n[s]));Require(cache[s]->GetBinContent(x,y)==n[s],"TH2D copy changed stored count");a.count[s]+=n[s];a.variance[s]+=n[s];if(x>0&&x<=source[0]->GetNbinsX()&&y>0&&y<=source[0]->GetNbinsY())a.regular[s]+=n[s];}++a.cells;
  }
  a.entries[4]=a.count[4];for(int s=0;s<5;++s)cache[s]->SetEntries(a.entries[s]);
}
void WriteAudit(const Audit& a,const char* output,const Settings& c){
  std::ofstream f((std::string(output)+"/native_audit.json").c_str());Require(bool(f),"Cannot create native audit");f<<std::setprecision(17);
  f<<"{\"schemaVersion\":1,\"cellsIncludingFlows\":"<<a.cells<<",\"axesExactlyEqual\":true,\"tagSumSubsetOfAll\":true,\"minimumRemainingCount\":"<<a.minimumRemaining<<",\"maximumSourceVarianceRelativeDifference\":"<<a.maximumVarianceRelative<<",\"nativeAxis\":{\"rBins\":600,\"rRange\":[0,6],\"rBinWidth\":0.01,\"dedxBins\":1200,\"dedxRange\":[0,120],\"dedxBinWidth\":0.1},\"sourceSumw2Sizes\":[0,0,0,0],\"cacheClass\":\"TH2D\",\"cacheErrorPolicy\":\"sqrt(stored count), explicit Sumw2; REMAINING is correlated subset subtraction, variance=remaining stored N\",\"cacheCountErrorFlowPreserved\":true,\"histograms\":[";
  for(int s=0;s<5;++s){if(s)f<<',';f<<"{\"cacheKey\":\""<<names[s]<<"\",\"sourceKey\":\""<<(s<4?inputKeys[s]:"derived ALL-TAG_D-TAG_T-TAG_HE3")<<"\",\"storedCount\":"<<a.count[s]<<",\"varianceSum\":"<<a.variance[s]<<",\"normalCount\":"<<a.regular[s]<<",\"flowCount\":"<<a.count[s]-a.regular[s]<<",\"entries\":"<<a.entries[s]<<",\"storedMinusEntries\":"<<a.count[s]-a.entries[s]<<",\"entriesMeaning\":\""<<(s<4?"preserved source GetEntries":"derived stored sum, not original track fill count")<<"\"}";}
  f<<"],\"writtenCacheReadbackExact\":true,\"noFit\":true,\"productionEnabled\":false,\"tagDIsPureHe4\":false,\"coarseMomentsAreGaussianFit\":false}\n";f.close();Require(f.good(),"Cannot finish native audit");
}
void VerifyWrittenCache(TH2D* expected[5],const char* output){
  TFile readback(TString(output)+"/cache.root","READ");Require(!readback.IsZombie(),"Cannot read back native cache");
  Require(readback.GetListOfKeys()->GetSize()==5,"Written cache must have exactly five flat histograms");
  for(int s=0;s<5;++s){TH2D* actual=dynamic_cast<TH2D*>(readback.Get(names[s]));Require(actual!=0,"Written cache flat key missing or wrong class");Axis(expected[s]->GetXaxis(),actual->GetXaxis());Axis(expected[s]->GetYaxis(),actual->GetYaxis());Require(actual->GetEntries()==expected[s]->GetEntries(),"Written cache Entries changed");
    for(int x=0;x<=actual->GetNbinsX()+1;++x)for(int y=0;y<=actual->GetNbinsY()+1;++y)Require(actual->GetBinContent(x,y)==expected[s]->GetBinContent(x,y)&&actual->GetBinError(x,y)==expected[s]->GetBinError(x,y),"Written native count/error/flow differs from prepared cache");
  }readback.Close();
}
void Map(TH2D* input,int sample,const Settings& c,int mode,TCanvas& page,const TString& pdf,const char* output){
  page.Clear();page.cd();gPad->SetLeftMargin(Need(c,"mapMarginLeft"));gPad->SetRightMargin(Need(c,"mapMarginRight"));gPad->SetBottomMargin(Need(c,"mapMarginBottom"));gPad->SetTopMargin(Need(c,"mapMarginTop"));gPad->SetLogz();
  TH2D* h=static_cast<TH2D*>(input->Clone(Form("map_%d_%d",sample,mode)));h->SetDirectory(0);h->RebinY(Integer(c,"displayDedxRebinFactor"));
  const double xlo=Need(c,mode?"atlasRigidityMin":"overviewRigidityMin"),xhi=Need(c,mode?"atlasRigidityMax":"overviewRigidityMax"),ylo=Need(c,mode==2?"supplementalDedxMin":"fullDedxMin"),yhi=Need(c,mode==2?"supplementalDedxMax":"fullDedxMax");
  h->GetXaxis()->SetRangeUser(xlo,xhi);h->GetYaxis()->SetRangeUser(ylo,yhi);h->SetMinimum(1);h->GetXaxis()->SetTitleSize(Need(c,"mapAxisTitleSize"));h->GetYaxis()->SetTitleSize(Need(c,"mapAxisTitleSize"));h->GetXaxis()->SetLabelSize(Need(c,"mapAxisLabelSize"));h->GetYaxis()->SetLabelSize(Need(c,"mapAxisLabelSize"));h->GetXaxis()->SetTitleOffset(Need(c,"xTitleOffset"));h->GetYaxis()->SetTitleOffset(Need(c,"yTitleOffset"));h->Draw("COLZ");
  TLatex text;text.SetNDC();text.SetTextFont(42);text.SetTextSize(Need(c,"headerTextSize"));text.DrawLatex(Need(c,"mapMarginLeft"),.94,labels[sample]);text.SetTextSize(Need(c,"captionTextSize"));text.DrawLatex(Need(c,"mapMarginLeft"),.075,Form("DISPLAY ONLY R %.2f-%.2f; dE/dx %.1f-%.1f; native cache untouched",xlo,xhi,ylo,yhi));text.DrawLatex(Need(c,"mapMarginLeft"),.035,"Charge-inclusive pre-TPC QA; TOF subsets are not independent of ALL or isotope truth");
  page.Print(pdf);page.Print(TString(output)+Form("/colz_%s_%d.png",names[sample],mode));delete h;
}
struct Moment{double mode,count,mean,rms;Moment():mode(-999),count(0),mean(-999),rms(-999){}};
Moment ModeMoment(TH1D* h,double lo,double hi,const Settings& c){
  Moment m;double maximum=-1;int half=Integer(c,"modeSmoothingHalfBins");for(int b=1;b<=h->GetNbinsX();++b)if(h->GetBinCenter(b)>=lo&&h->GetBinCenter(b)<hi){double n=0;int bins=0;for(int q=std::max(1,b-half);q<=std::min(h->GetNbinsX(),b+half);++q){n+=h->GetBinContent(q);++bins;}n/=bins;if(n>maximum){maximum=n;m.mode=h->GetBinCenter(b);}}
  if(maximum<=0){m.mode=-999;return m;}double first=0,second=0;for(int b=1;b<=h->GetNbinsX();++b){double x=h->GetBinCenter(b),n=h->GetBinContent(b);if(x>=lo&&x<hi&&std::fabs(x-m.mode)<=Need(c,"modeMomentHalfWindow")){m.count+=n;first+=n*x;second+=n*x*x;}}if(m.count>0){m.mean=first/m.count;m.rms=std::sqrt(std::max(0.,second/m.count-m.mean*m.mean));}return m;
}
void Projection(TH1D* source,int sample,double lo,double hi,const Settings& c,bool log){
  gPad->SetLeftMargin(Need(c,"projectionMarginLeft"));gPad->SetRightMargin(Need(c,"projectionMarginRight"));gPad->SetBottomMargin(Need(c,"projectionMarginBottom"));gPad->SetTopMargin(Need(c,"projectionMarginTop"));gPad->SetLogy(log);
  TH1D* h=static_cast<TH1D*>(source->Clone(Form("draw_%s_%g_%d",names[sample],lo,int(log))));h->SetDirectory(0);h->SetBit(TObject::kCanDelete);h->GetXaxis()->SetRange(1,h->GetNbinsX());double maximum=1;for(int b=1;b<=h->GetNbinsX();++b)maximum=std::max(maximum,h->GetBinContent(b)+h->GetBinError(b));
  h->SetMinimum(log?Need(c,"logMinimum"):0);h->SetMaximum(maximum*Need(c,log?"logMaximumFactor":"linearMaximumFactor"));h->SetTitle(Form(";TPC dE/dx (keV/cm);Tracks / %.2f keV/cm",h->GetBinWidth(1)));h->SetMarkerStyle(20);h->SetMarkerSize(Need(c,"markerSize"));h->GetXaxis()->SetTitleSize(Need(c,"projectionAxisTitleSize"));h->GetYaxis()->SetTitleSize(Need(c,"projectionAxisTitleSize"));h->GetXaxis()->SetLabelSize(Need(c,"projectionAxisLabelSize"));h->GetYaxis()->SetLabelSize(Need(c,"projectionAxisLabelSize"));h->GetXaxis()->SetTitleOffset(Need(c,"xTitleOffset"));h->GetYaxis()->SetTitleOffset(Need(c,"yTitleOffset"));h->Draw("E1");
  TLatex t;t.SetNDC();t.SetTextFont(42);t.SetTextSize(Need(c,"headerTextSize"));t.DrawLatex(Need(c,"projectionMarginLeft"),.95,Form("%s | R %.2f-%.2f | %s",names[sample],lo,hi,log?"log":"linear"));t.SetTextSize(Need(c,"captionTextSize"));t.DrawLatex(Need(c,"projectionMarginLeft"),.095,labels[sample]);t.DrawLatex(Need(c,"projectionMarginLeft"),.055,Form("Stored tracks %.6g; y-flow %.0f; full physical axis 0-120",h->Integral(1,h->GetNbinsX()),h->GetBinContent(0)+h->GetBinContent(h->GetNbinsX()+1)));t.DrawLatex(Need(c,"projectionMarginLeft"),.022,"No fit; strongest mode / window moments are NOT isotope identification");
}
void Coarse(TH2D* cache[5],const Settings& c,const char* output){
  TFile p((std::string(output)+"/coarse_projections.root").c_str(),"NEW");Require(!p.IsZombie(),"Cannot create coarse projection ROOT");std::ofstream tsv((std::string(output)+"/coarse_statistics.tsv").c_str());Require(bool(tsv),"Cannot create coarse statistics");tsv<<std::setprecision(17)<<"index\tlo\thi\tfirstNativeR\tlastNativeR\tdiagnosticInterval\tsample\tnormalCount\tunderflow\toverflow\tvarianceSum\tsupplementCount\tglobalMode\tglobalModeWindowCount\tglobalModeWindowMean\tglobalModeWindowRms\tsupplementMode\tsupplementModeWindowCount\tsupplementModeWindowMean\tsupplementModeWindowRms\n";
  TCanvas atlas("nuclear_input_atlas","Four original QA samples / no fit",Integer(c,"atlasCanvasWidth"),Integer(c,"atlasCanvasHeight"));TString pdf=TString(output)+"/coarse_full_dedx_atlas.pdf";atlas.Print(pdf+"[");int previousLast=-1,groups=int(std::floor((Need(c,"atlasRigidityMax")-Need(c,"atlasRigidityMin"))/Need(c,"coarseRigidityStep")+.5));
  for(int group=0;group<groups;++group){double lo=Need(c,"atlasRigidityMin")+group*Need(c,"coarseRigidityStep"),hi=lo+Need(c,"coarseRigidityStep"),tol=Need(c,"axisRelativeTolerance");int first=cache[0]->GetXaxis()->FindFixBin(lo+tol),last=cache[0]->GetXaxis()->FindFixBin(hi-tol);Require(first>=1&&last<=cache[0]->GetNbinsX()&&last>=first&&(group==0||first==previousLast+1),"Coarse R groups overlap/gap or include flow");Require(Near(lo,cache[0]->GetXaxis()->GetBinLowEdge(first),tol)&&Near(hi,cache[0]->GetXaxis()->GetBinUpEdge(last),tol),"Coarse R boundaries differ from native");previousLast=last;atlas.Clear();atlas.Divide(2,4);TH1D* display[5];
    for(int sample=0;sample<5;++sample){TH1D* native=cache[sample]->ProjectionY(Form("coarse_native_%s_%03d",names[sample],group),first,last,"e");native->SetDirectory(0);p.cd();native->Write();display[sample]=static_cast<TH1D*>(native->Clone(Form("coarse_display_%s_%03d",names[sample],group)));display[sample]->SetDirectory(0);int factor=Integer(c,"displayDedxRebinFactor");display[sample]->Rebin(factor);
      for(int b=0;b<=display[sample]->GetNbinsX()+1;++b){double n=0,v=0;if(b==0||b==display[sample]->GetNbinsX()+1){int q=b==0?0:native->GetNbinsX()+1;n=native->GetBinContent(q);v=std::pow(native->GetBinError(q),2);}else for(int q=1+(b-1)*factor;q<=b*factor;++q){n+=native->GetBinContent(q);v+=std::pow(native->GetBinError(q),2);}Require(display[sample]->GetBinContent(b)==n&&Near(std::pow(display[sample]->GetBinError(b),2),v,Need(c,"varianceRelativeTolerance")),"Coarse native/display count/variance/flow changed");}display[sample]->Write();
      Moment full=ModeMoment(display[sample],Need(c,"fullDedxMin"),Need(c,"fullDedxMax"),c),supp=ModeMoment(display[sample],Need(c,"supplementalDedxMin"),Need(c,"supplementalDedxMax"),c);double suppCount=0,variance=0;for(int b=0;b<=display[sample]->GetNbinsX()+1;++b){variance+=std::pow(display[sample]->GetBinError(b),2);if(b>0&&b<=display[sample]->GetNbinsX()&&display[sample]->GetBinCenter(b)>=Need(c,"supplementalDedxMin")&&display[sample]->GetBinCenter(b)<Need(c,"supplementalDedxMax"))suppCount+=display[sample]->GetBinContent(b);}
      bool diagnostic=lo>=Need(c,"diagnosticRigidityMin")-tol&&hi<=Need(c,"diagnosticRigidityMax")+tol;tsv<<group<<'\t'<<lo<<'\t'<<hi<<'\t'<<first<<'\t'<<last<<'\t'<<int(diagnostic)<<'\t'<<names[sample]<<'\t'<<display[sample]->Integral(1,display[sample]->GetNbinsX())<<'\t'<<display[sample]->GetBinContent(0)<<'\t'<<display[sample]->GetBinContent(display[sample]->GetNbinsX()+1)<<'\t'<<variance<<'\t'<<suppCount<<'\t'<<full.mode<<'\t'<<full.count<<'\t'<<full.mean<<'\t'<<full.rms<<'\t'<<supp.mode<<'\t'<<supp.count<<'\t'<<supp.mean<<'\t'<<supp.rms<<'\n';delete native;
      if(sample<4){atlas.cd(2*sample+1);Projection(display[sample],sample,lo,hi,c,true);atlas.cd(2*sample+2);Projection(display[sample],sample,lo,hi,c,false);}
    }atlas.Print(pdf);if(group==0||Near(lo,.4,Need(c,"axisRelativeTolerance"))||Near(lo,.8,Need(c,"axisRelativeTolerance")))atlas.Print(TString(output)+Form("/coarse_atlas_R_%03d.png",int(lo*1000+.5)));for(int s=0;s<5;++s)delete display[s];
  }atlas.Print(pdf+"]");p.Close();tsv.close();Require(tsv.good(),"Cannot finish coarse statistics");
}
}
#endif

void ExtractTofInputs(const char* input,const char* output,const char* settings){
#ifndef __CINT__
  using namespace NuclearTofInputs;
  try{Settings c=Read(settings);Validate(c);Require(gSystem->AccessPathName(TString(output)+"/cache.root"),"Refusing overwrite of cache");TFile source(input,"READ");Require(!source.IsZombie(),"Cannot READ original ROOT");TH2* hist[4];for(int s=0;s<4;++s){hist[s]=dynamic_cast<TH2*>(source.Get(inputKeys[s]));Require(hist[s]!=0,std::string("Missing native source: ")+inputKeys[s]);}TH2D* cache[5]={0,0,0,0,0};Audit audit;Native(hist,cache,c,audit);
    TFile product(TString(output)+"/cache.root","NEW");Require(!product.IsZombie(),"Cannot create cache ROOT");for(int s=0;s<5;++s){product.cd();Require(cache[s]->Write()>0,"Native cache histogram write failed");}product.Close();VerifyWrittenCache(cache,output);WriteAudit(audit,output,c);
    gStyle->SetOptStat(0);gStyle->SetOptFit(0);gStyle->SetPalette(1);TGaxis::SetMaxDigits(Integer(c,"axisMaxDigits"));TCanvas map("nuclear_input_maps","Native nuclear TOF diagnostics",Integer(c,"mapCanvasWidth"),Integer(c,"mapCanvasHeight"));TString maps=TString(output)+"/source_colz_atlas.pdf";map.Print(maps+"[");for(int s=0;s<4;++s)for(int mode=0;mode<3;++mode)Map(cache[s],s,c,mode,map,maps,output);map.Print(maps+"]");Coarse(cache,c,output);source.Close();for(int s=0;s<5;++s)delete cache[s];
    std::cout<<"NUCLEAR_TOF_INPUTS_COMPLETE cells="<<audit.cells<<"; fixed source READ once, native stored counts preserved; no fit"<<std::endl;
  }catch(const std::exception& e){std::cerr<<"NUCLEAR_TOF_INPUTS_ERROR "<<e.what()<<std::endl;gSystem->Exit(1);}
#endif
}
