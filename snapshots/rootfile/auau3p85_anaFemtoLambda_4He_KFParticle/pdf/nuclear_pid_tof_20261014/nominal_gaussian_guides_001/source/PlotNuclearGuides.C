void PlotNuclearGuides(const char* =0,const char* =0,const char* =0,const char* =0,const char* =0);
#ifndef __CINT__
#include "TFile.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TGaxis.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TMath.h"
#include "TPad.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
namespace NuclearGuidePlots {
typedef std::map<std::string,std::string> Row;
const char* species[]={"d","t","he3","he4"};
const int colors[]={kOrange+7,kBlue+1,kCyan+2,kGreen+2};
void Need(bool b,const std::string& m){if(!b)throw std::runtime_error(m);}
std::vector<std::string> Split(const std::string& s){std::vector<std::string> v;std::string x;std::istringstream in(s);while(std::getline(in,x,'\t'))v.push_back(x);return v;}
std::string Txt(const Row& r,const std::string& k){Row::const_iterator i=r.find(k);Need(i!=r.end(),"Missing "+k);return i->second;}
double Num(const Row& r,const std::string& k){std::istringstream in(Txt(r,k));double x=0;in>>x;in>>std::ws;Need(in.eof()&&TMath::Finite(x),"Bad number "+k);return x;}
int Int(const Row& r,const std::string& k){double x=Num(r,k);Need(x==std::floor(x),"Bad integer "+k);return int(x);}
bool Near(double a,double b,double t){return std::fabs(a-b)<=t*std::max(1.,std::max(std::fabs(a),std::fabs(b)));}
std::vector<Row> Rows(const char* path){std::ifstream f(path);Need(bool(f),"Missing manifest");std::string s;Need(bool(std::getline(f,s)),"Empty manifest");std::vector<std::string> h=Split(s);Need(std::set<std::string>(h.begin(),h.end()).size()==h.size(),"Duplicate columns");std::vector<Row> out;while(std::getline(f,s)){std::vector<std::string> v=Split(s);Need(v.size()==h.size(),"Bad column count");Row r;for(unsigned j=0;j<h.size();++j)r[h[j]]=v[j];out.push_back(r);}Need(f.eof(),"Unreadable manifest");return out;}
Row Settings(const char* path){std::ifstream f(path);Need(bool(f),"Missing settings");Row out;std::string s;while(std::getline(f,s)){std::vector<std::string> v=Split(s);Need(v.size()==2&&out.insert(std::make_pair(v[0],v[1])).second,"Malformed setting");}Need(f.eof(),"Unreadable settings");return out;}
int Species(const std::string& s){for(int j=0;j<4;++j)if(s==species[j])return j;throw std::runtime_error("Unknown species "+s);}
bool Show(const Row& r,int selected){return selected<0||Species(Txt(r,"species"))==selected;}
void Curves(const std::vector<Row>& rows,int selected,int family,const Row& c,std::vector<TGraph*>& keep){
  std::map<std::string,std::vector<Row> > groups;
  const char* names[]={"legacy","working","candidate"};
  for(unsigned j=0;j<rows.size();++j)if(Show(rows[j],selected)&&Txt(rows[j],"family")==names[family])groups[Txt(rows[j],"species")+":"+Txt(rows[j],"group")].push_back(rows[j]);
  for(std::map<std::string,std::vector<Row> >::const_iterator g=groups.begin();g!=groups.end();++g){const std::vector<Row>& rr=g->second;Need(rr.size()>=2,"Curve group cannot be approximated by isolated point");int color=family==0?kGray+2:family==1?kMagenta+1:colors[Species(Txt(rr[0],"species"))];for(int band=-1;band<=1;++band){TGraph* line=new TGraph();line->SetBit(TObject::kCanDelete);line->SetLineColor(color);line->SetLineStyle(family==0?3:family==1?7:band?2:1);line->SetLineWidth(family==2&&band==0?3:1);for(unsigned j=0;j<rr.size();++j)line->SetPoint(j,Num(rr[j],"r"),Num(rr[j],"mean")+band*Num(c,"guideNSigma")*Num(rr[j],"sigma"));line->Draw("L same");keep.push_back(line);}}
}
void Panel(TH2D* data,const std::vector<Row>& points,const std::vector<Row>& curves,const Row& c,int selected,bool low,const char* caption,const char* warning){
  gPad->SetLeftMargin(Num(c,"marginLeft"));gPad->SetRightMargin(Num(c,"marginRight"));gPad->SetBottomMargin(Num(c,"marginBottom"));gPad->SetTopMargin(Num(c,"marginTop"));gPad->SetLogz();
  TH2D* h=static_cast<TH2D*>(data->Clone(Form("display_%s_%d_%d",data->GetName(),selected,int(low))));h->SetDirectory(0);h->SetBit(TObject::kCanDelete);h->GetXaxis()->SetRangeUser(Num(c,low?"lowRMin":"fullRMin"),Num(c,low?"lowRMax":"fullRMax"));h->GetYaxis()->SetRangeUser(Num(c,"visibleDedxMin"),Num(c,"visibleDedxMax"));h->SetMinimum(Num(c,"logZMinimum"));h->SetTitle(";Rigidity |p/q| (GeV/c);TPC dE/dx (keV/cm);Tracks / display cell");
  TAxis* axes[]={h->GetXaxis(),h->GetYaxis(),h->GetZaxis()};for(int j=0;j<3;++j){axes[j]->SetTitleFont(Int(c,"axisFont"));axes[j]->SetLabelFont(Int(c,"axisFont"));axes[j]->SetTitleSize(Num(c,"axisTitleSize"));axes[j]->SetLabelSize(Num(c,"axisLabelSize"));}axes[0]->SetTitleOffset(Num(c,"xTitleOffset"));axes[1]->SetTitleOffset(Num(c,"yTitleOffset"));axes[2]->SetTitleOffset(Num(c,"zTitleOffset"));h->Draw("colz");
  std::vector<TGraph*> refs;for(int f=0;f<3;++f)Curves(curves,selected,f,c,refs);
  TGraphErrors* good[4];for(int j=0;j<4;++j){good[j]=new TGraphErrors();good[j]->SetBit(TObject::kCanDelete);good[j]->SetMarkerStyle(20);good[j]->SetMarkerColor(colors[j]);good[j]->SetLineColor(colors[j]);good[j]->SetMarkerSize(Num(c,"markerSize"));}
  TGraphErrors* failed=new TGraphErrors();failed->SetBit(TObject::kCanDelete);failed->SetMarkerStyle(5);failed->SetMarkerColor(kRed+1);failed->SetLineColor(kRed+1);failed->SetMarkerSize(Num(c,"markerSize"));int nGood=0,nBad=0;
  for(unsigned i=0;i<points.size();++i)if(Show(points[i],selected)){const Row& p=points[i];TGraphErrors* g=Int(p,"acceptedNumericalShape")?good[Species(Txt(p,"species"))]:failed;int n=g->GetN();g->SetPoint(n,Num(p,"r"),Num(p,"mean"));g->SetPointError(n,0,Num(c,"guideNSigma")*Num(p,"sigma"));if(Int(p,"acceptedNumericalShape"))++nGood;else ++nBad;}
  // P/E1 only: failed centers and width bars are never connected into a curve.
  if(failed->GetN())failed->Draw("P E1 same");for(int j=0;j<4;++j)if(good[j]->GetN())good[j]->Draw("P E1 same");
  // Legend keys must outlive this function until the canvas is printed.
  TGraph* legacyKey=new TGraph();TGraph* workingKey=new TGraph();TGraph* candidateKey=new TGraph();legacyKey->SetLineColor(kGray+2);legacyKey->SetLineStyle(3);workingKey->SetLineColor(kMagenta+1);workingKey->SetLineStyle(7);candidateKey->SetLineColor(kBlack);candidateKey->SetLineWidth(3);
  TLegend legend(.12,.79,.90,.91);legend.SetFillStyle(0);legend.SetBorderSize(0);legend.SetNColumns(2);legend.SetTextSize(Num(c,"legendTextSize"));legend.AddEntry(good[selected<0?0:selected],"Free centers / numerical gates only","pe");legend.AddEntry(failed,"Failed / unsupported: unlinked","pe");legend.AddEntry(legacyKey,"Legacy center and +/-2 widths","l");legend.AddEntry(workingKey,"Working knot provider (reference)","l");legend.AddEntry(candidateKey,"PCHIP only supported native intervals","l");legend.DrawClone();
  TLatex text;text.SetNDC();text.SetTextFont(Int(c,"axisFont"));text.SetTextSize(Num(c,"headerTextSize"));text.DrawLatex(.12,.98,caption);text.SetTextSize(Num(c,"captionTextSize"));text.SetTextColor(kRed+1);text.DrawLatex(.12,.945,warning);text.SetTextColor(kBlack);text.DrawLatex(.12,.115,Form("%s R range; saved counts, display-only 5-bin dE/dx sum",low?"Additional low":"Main"));text.DrawLatex(.12,.086,"Bars = +/-2 projected Gaussian widths, NOT uncertainty / purity / isotope truth");text.DrawLatex(.12,.057,"PCHIP widths are R-unfolded SEEDS; unsupported gaps are never bridged");text.DrawLatex(.12,.028,selected<0?"Colors: d orange, t blue, He3 cyan, He4 green; He4: NO supported curve":"TAG acceptance / efficiency biases remain; not an independent ALL sample");
}
}
#endif
void PlotNuclearGuides(const char* cache,const char* points,const char* curves,const char* settings,const char* output){
#ifndef __CINT__
  using namespace NuclearGuidePlots;
  if(!cache&&!points&&!curves&&!settings&&!output){std::cout<<"NUCLEAR_GUIDES_COMPILE_OK"<<std::endl;return;}
  try{Need(cache&&points&&curves&&settings&&output,"All arguments explicit");Row c=Settings(settings);Need(Int(c,"schemaVersion")==1&&Int(c,"expectedRBins")==600&&Int(c,"expectedDedxBins")==1200&&Int(c,"dedxDisplayRebin")==5,"Unsupported native/display schema");Need(Num(c,"guideNSigma")==2&&Num(c,"nativeDedxMin")==0&&Num(c,"nativeDedxMax")==120&&Num(c,"visibleDedxMin")==0&&Num(c,"visibleDedxMax")==120,"Full physical 0-120 axis and two-width guide required");std::vector<Row> p=Rows(points),q=Rows(curves);Need(!p.empty(),"Missing free-fit points");TFile in(cache,"READ");Need(!in.IsZombie(),"Unreadable native cache");const char* keys[]={"ALL","TAG_D","TAG_T","TAG_HE3","REMAINING"};TH2D* maps[5];std::ofstream audit((std::string(output)+"/map_audit.tsv").c_str());Need(bool(audit),"Cannot create map audit");audit<<std::setprecision(17)<<"sample\tnativeCells\tnativeCountIncludingFlows\tdisplayCountIncludingFlows\tmaximumDisplayVarianceRelativeDifference\n";
    for(int k=0;k<5;++k){TH2D* src=dynamic_cast<TH2D*>(in.Get(keys[k]));Need(src&&src->GetNbinsX()==600&&src->GetNbinsY()==1200&&src->GetXaxis()->GetXmin()==Num(c,"nativeRMin")&&src->GetXaxis()->GetXmax()==Num(c,"nativeRMax")&&src->GetYaxis()->GetXmin()==0&&src->GetYaxis()->GetXmax()==120,"Native cache dimensions differ");maps[k]=static_cast<TH2D*>(src->Clone(Form("map_%s",keys[k])));maps[k]->SetDirectory(0);maps[k]->RebinY(5);double nativeSum=0,displaySum=0,maxDiff=0;for(int x=0;x<=601;++x){for(int y=0;y<=1201;++y){double n=src->GetBinContent(x,y);Need(TMath::Finite(n)&&n>=0&&Near(std::pow(src->GetBinError(x,y),2),n,Num(c,"axisRelativeTolerance")),"Invalid native count/error");nativeSum+=n;}for(int y=0;y<=241;++y){int first=y==0?0:y==241?1201:1+5*(y-1),last=(y==0||y==241)?first:first+4;double n=0,v=0;for(int b=first;b<=last;++b){n+=src->GetBinContent(x,b);v+=std::pow(src->GetBinError(x,b),2);}Need(maps[k]->GetBinContent(x,y)==n,"Display count/flow changed");double delta=std::fabs(std::pow(maps[k]->GetBinError(x,y),2)-v)/std::max(1.,v);maxDiff=std::max(maxDiff,delta);Need(delta<=Num(c,"axisRelativeTolerance"),"Display variance changed");displaySum+=maps[k]->GetBinContent(x,y);}}Need(displaySum==nativeSum,"Display all-cell count changed");audit<<keys[k]<<'\t'<<src->GetNcells()<<'\t'<<nativeSum<<'\t'<<displaySum<<'\t'<<maxDiff<<'\n';}
    audit.close();Need(audit.good(),"Cannot finish audit");gStyle->SetOptStat(0);gStyle->SetOptFit(0);gStyle->SetPaperSize(Num(c,"paperWidthCm"),Num(c,"paperHeightCm"));TGaxis::SetMaxDigits(Int(c,"axisMaxDigits"));TCanvas page("nuclear_guides","Nominal Gaussian guide comparison",Int(c,"canvasWidth"),Int(c,"canvasHeight"));TString pdf=TString(output)+"/nominal_gaussian_colz_guides.pdf";Need(gSystem->AccessPathName(pdf),"Refusing overwrite");page.Print(pdf+"[");int mapIndex[]={0,1,2,3,1},selection[]={-1,0,1,2,3};const char* titles[]={"ALL / four nominal species","TAG_D / nominal d branch","TAG_T / nominal t branch","TAG_HE3 / nominal He3 branch","SAME TAG_D / nominal He4-compatible branch"};const char* warnings[]={"REFERENCE ONLY; no complete supported 4-species calibration","D / He4-compatible TOF gate; NOT pure d","TOF-conditional sample; isotope identity not validated","TOF-conditional sample; isotope identity not validated","NOT independent He4 tag; NO supported He4 curve"};for(int i=0;i<5;++i){page.Clear();page.Divide(2,1);for(int low=0;low<2;++low){page.cd(low+1);Panel(maps[mapIndex[i]],p,q,c,selection[i],low,titles[i],warnings[i]);}page.Print(pdf);}page.Print(pdf+"]");for(int k=0;k<5;++k)delete maps[k];in.Close();std::cout<<"NUCLEAR_GUIDES_COMPLETE pages=5; no fit/adoption"<<std::endl;
  }catch(const std::exception& e){std::cerr<<"PlotNuclearGuides: "<<e.what()<<std::endl;gSystem->Exit(2);}
#endif
}
