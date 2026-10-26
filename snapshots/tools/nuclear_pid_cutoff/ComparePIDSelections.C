// Read-only configuration comparison; no fit and no production mutation.
void ComparePIDSelections(const char*,const char*,const char*,const char*,const char*,const char*,const char*);
#ifndef __CINT__
#include "TFile.h"
#include "TH2.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TPad.h"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
namespace PidCompare {
typedef std::map<std::string,std::string> Row;
typedef std::vector<Row> Rows;
void Require(bool ok,const std::string& reason){if(!ok)throw std::runtime_error(reason);}
bool Finite(double x){return std::isfinite(x);}
double N(const Row& r,const std::string& key){Row::const_iterator it=r.find(key);Require(it!=r.end(),"Missing column "+key);char* end=0;double value=std::strtod(it->second.c_str(),&end);Require(end&&*end=='\0',"Malformed numeric column "+key);return value;}
std::vector<std::string> Split(const std::string& line){std::vector<std::string> out;std::istringstream input(line);std::string value;while(std::getline(input,value,'\t'))out.push_back(value);return out;}
Rows Read(const char* path){std::ifstream input(path);Require(bool(input),"Cannot read TSV");std::string line;Require(bool(std::getline(input,line)),"Missing TSV header");std::vector<std::string> columns=Split(line);Rows rows;while(std::getline(input,line)){std::vector<std::string> fields=Split(line);Require(fields.size()==columns.size(),"TSV field count differs");Row row;for(unsigned i=0;i<fields.size();++i)Require(row.insert(std::make_pair(columns[i],fields[i])).second,"Duplicate TSV column");rows.push_back(row);}Require(input.eof()&&!rows.empty(),"Empty TSV");return rows;}
Row Settings(const char* path){Rows rows=Read(path);Row result;for(unsigned i=0;i<rows.size();++i)Require(result.insert(std::make_pair(rows[i]["key"],rows[i]["value"])).second,"Duplicate plot setting");return result;}
const char* Species[]={"d","t","he3","he4"};
const char* Labels[]={"d","t","^{3}He","^{4}He"};
const char* Profiles[]={"oldCanonical","oldLegacy","new"};
const char* ProfileLabels[]={"Old calibrated (canonical He4)","Old legacy (dedicated d/t/He3)","New common calibrated"};
const int Colors[]={kBlack,kBlue+1,kRed+1};
const int Styles[]={3,2,1};
void Text(double x,double y,const char* value,double size=.03){TLatex t;t.SetNDC();t.SetTextSize(size);t.SetTextFont(42);t.DrawLatex(x,y,value);}
void Pad(){gPad->SetLeftMargin(.13);gPad->SetRightMargin(.05);gPad->SetBottomMargin(.13);gPad->SetTopMargin(.12);gPad->SetGridx();gPad->SetGridy();}
double Value(const Row& row,int profile,int kind,double ns){double mean=N(row,std::string(Profiles[profile])+"Mean"),sigma=N(row,std::string(Profiles[profile])+"Sigma");if(!Finite(mean)||!Finite(sigma)||sigma<=0)return std::numeric_limits<double>::quiet_NaN();if(kind==0)return mean;if(kind==1)return sigma;return mean+(kind==2?-1:1)*ns*sigma;}
TGraph* Curve(const Rows& rows,const char* species,int profile,int kind,double xmin,double xmax,double ns){TGraph* g=new TGraph();g->SetBit(TObject::kCanDelete);for(unsigned i=0;i<rows.size();++i)if(rows[i].find("species")->second==species){double x=N(rows[i],"R"),y=Value(rows[i],profile,kind,ns);if(x>=xmin&&x<=xmax&&Finite(y))g->SetPoint(g->GetN(),x,y);}g->SetLineColor(Colors[profile]);g->SetLineStyle(Styles[profile]);g->SetLineWidth(kind>1?1:2);return g;}
void ResponsePage(TCanvas& canvas,const Rows& rows,int species,double xmin,double xmax,const Row& cfg,const TString& pdf,const TString& output,const char* tag){
  canvas.Clear();canvas.Divide(2,2);double ns=N(cfg,"gateNSigma");
  for(int panel=0;panel<4;++panel){canvas.cd(panel+1);Pad();double low=1.e300,high=-1.e300;
    for(unsigned i=0;i<rows.size();++i)if(rows[i].find("species")->second==Species[species]&&N(rows[i],"R")>=xmin&&N(rows[i],"R")<=xmax){const Row& row=rows[i];double nm=N(row,"newMean"),nsig=N(row,"newSigma");for(int p=0;p<3;++p){double v=panel<2?Value(row,p,panel,ns):(panel==2?(N(row,std::string(Profiles[p])+"Mean")-nm)/nsig:N(row,std::string(Profiles[p])+"Sigma")/nsig);if(Finite(v)){low=std::min(low,v);high=std::max(high,v);}if(panel==0){double a=Value(row,p,2,ns),b=Value(row,p,3,ns);if(Finite(a)&&Finite(b)){low=std::min(low,a);high=std::max(high,b);}}}}
    Require(Finite(low)&&Finite(high)&&high>=low,"No response points");double span=high-low;if(span<1.e-10)span=panel==3?.1:std::max(.01,std::fabs(high)*.1);low-=.12*span;high+=.22*span;if(panel<2)low=std::max(0.,low);
    const char* titles[]={"Configured dE/dx center and #pm2#sigma edges","Configured dE/dx width","Center difference / new #sigma","Width ratio to new"};const char* axes[]={"dE/dx (keV/cm)","#sigma_{dE/dx} (keV/cm)","(#mu_{old}-#mu_{new})/#sigma_{new}","#sigma_{old}/#sigma_{new}"};
    TH1F* frame=gPad->DrawFrame(xmin,low,xmax,high);frame->SetTitle(Form(";R=|p/q| (GeV/c);%s",axes[panel]));
    TLegend* legend=new TLegend(.16,.69,.82,.86);legend->SetBit(TObject::kCanDelete);legend->SetBorderSize(0);legend->SetFillStyle(0);legend->SetTextSize(.024);
    for(int p=0;p<(panel<2?3:2);++p){TGraph* graph=0;if(panel<2){graph=Curve(rows,Species[species],p,panel,xmin,xmax,ns);if(panel==0)for(int edge=2;edge<=3;++edge){TGraph* band=Curve(rows,Species[species],p,edge,xmin,xmax,ns);band->Draw("L SAME");}}
      else{graph=new TGraph();graph->SetBit(TObject::kCanDelete);graph->SetLineColor(Colors[p]);graph->SetLineStyle(Styles[p]);graph->SetLineWidth(2);for(unsigned i=0;i<rows.size();++i)if(rows[i].find("species")->second==Species[species]){double x=N(rows[i],"R"),nm=N(rows[i],"newMean"),nsig=N(rows[i],"newSigma"),old=N(rows[i],std::string(Profiles[p])+(panel==2?"Mean":"Sigma"));double y=panel==2?(old-nm)/nsig:old/nsig;if(x>=xmin&&x<=xmax&&Finite(y))graph->SetPoint(graph->GetN(),x,y);}}
      graph->Draw("L SAME");legend->AddEntry(graph,ProfileLabels[p],"l");}
    legend->Draw();Text(.13,.94,Form("%s: %s",Labels[species],titles[panel]),.032);Text(.14,.025,"Response functions only; domain/veto and TOF are separate selections",.023);
  }
  canvas.cd();canvas.Print(pdf);canvas.Print(output+Form("/response_%s_%s.png",Species[species],tag));
}
void Line(double x1,double y1,double x2,double y2,int color,int style){TLine* line=new TLine(x1,y1,x2,y2);line->SetBit(TObject::kCanDelete);line->SetLineColor(color);line->SetLineStyle(style);line->SetLineWidth(2);line->Draw();}
void Maps(TCanvas& canvas,TH2& all,TH2& m2,const Rows& rows,const Rows& windows,const Row& cfg,const TString& pdf,const TString& output){
  canvas.Clear();canvas.Divide(2,1);const int speciesColors[]={kBlue+1,kMagenta+1,kGreen+2,kRed+1};
  for(int panel=0;panel<2;++panel){canvas.cd(panel+1);Pad();gPad->SetRightMargin(.14);gPad->SetLogz();TH2* h=dynamic_cast<TH2*>(all.Clone(Form("pid_compare_all_%d",panel)));h->SetDirectory(0);h->SetBit(TObject::kCanDelete);h->GetXaxis()->SetRangeUser(N(cfg,"rigidityMin"),N(cfg,"rigidityMax"));h->GetYaxis()->SetRangeUser(0,N(cfg,"dedxDisplayMax"));h->SetMinimum(1);h->SetTitle(";R=|p/q| (GeV/c);TPC dE/dx (keV/cm)");h->Draw("COLZ");for(int s=0;s<4;++s)for(int kind=0;kind<4;++kind)if(kind!=1){TGraph* g=Curve(rows,Species[s],panel?2:1,kind,N(cfg,"rigidityMin"),N(cfg,"rigidityMax"),N(cfg,"gateNSigma"));g->SetLineColor(speciesColors[s]);g->SetLineStyle(kind?2:1);g->Draw("L SAME");}Text(.13,.94,panel?"New common calibrated TPC curves":"Old dedicated legacy TPC curves");Text(.14,.04,"ALL pre-PID QA; center and #pm2#sigma, not fitted confidence bands",.021);}
  canvas.Print(pdf);canvas.Print(output+"/tpc_maps.png");
  canvas.Clear();canvas.Divide(2,1);
  for(int panel=0;panel<2;++panel){canvas.cd(panel+1);Pad();gPad->SetRightMargin(.14);gPad->SetLogz();TH2* h=dynamic_cast<TH2*>(m2.Clone(Form("pid_compare_m2_%d",panel)));h->SetDirectory(0);h->SetBit(TObject::kCanDelete);double xmin=panel?N(cfg,"tofThreshold"):0,xmax=panel?N(cfg,"rigidityMax"):m2.GetXaxis()->GetXmax();h->GetXaxis()->SetRangeUser(xmin,xmax);h->GetYaxis()->SetRangeUser(N(cfg,panel?"m2DetailMin":"m2DisplayMin"),N(cfg,panel?"m2DetailMax":"m2DisplayMax"));h->SetMinimum(1);h->SetTitle(";R=|p/q| (GeV/c);TOF (m/q)^{2} (GeV^{2}/c^{4})");h->Draw("COLZ");
    for(unsigned s=0;s<windows.size();++s)for(int edge=0;edge<2;++edge){double y=N(windows[s],edge?"studiedMax":"studiedMin");Line(N(cfg,"tofThreshold"),y,xmax,y,speciesColors[s],s==3?2:1);}for(int edge=0;edge<2;++edge){double y=N(windows[3],edge?"currentConfiguredMax":"currentConfiguredMin");Line(N(cfg,"tofThreshold"),y,xmax,y,kBlack,3);}
    Text(.13,.94,panel?"TOF m^{2}: He3 / D / He4 gate detail":"TOF m^{2}: configured study windows");Text(.14,.08,"Red dashed He4 study = D window; black dotted He4 production = [3,4]",.021);Text(.14,.04,"Inclusive TOF QA only: no dE/dx condition, no new TOF fit, no purity claim",.021);}
  canvas.Print(pdf);canvas.Print(output+"/tof_m2_maps.png");
}
}
#endif
void ComparePIDSelections(const char* original,const char* cache,const char* outputPath,const char* curvePath,const char* windowPath,const char* configPath,const char* m2Key){
#ifndef __CINT__
  using namespace PidCompare;
  try{Rows curves=Read(curvePath),windows=Read(windowPath);Row cfg=Settings(configPath);Require(windows.size()==4,"Exactly 4 species needed");TString output(outputPath),pdf=output+"/pid_parameter_comparison.pdf";Require(gSystem->AccessPathName(pdf)&&gSystem->AccessPathName(output+"/tof_m2_native.root"),"Refusing output overwrite");
    TFile source(original,"READ");Require(!source.IsZombie()&&!source.TestBit(TFile::kRecovered),"Cannot open original ROOT");TH2* originalM2=dynamic_cast<TH2*>(source.Get(m2Key));Require(originalM2!=0,"Original single-key TOF map absent");TH2* m2=dynamic_cast<TH2*>(originalM2->Clone("TOF_M2_PRETPC"));m2->SetDirectory(0);source.Close();
    TFile product(output+"/tof_m2_native.root","CREATE");Require(!product.IsZombie(),"Cannot create native m2 cache");product.cd();Require(m2->Write()>0,"Cannot save native m2 cache");product.Close();TFile check(output+"/tof_m2_native.root","READ");TH2* readback=dynamic_cast<TH2*>(check.Get("TOF_M2_PRETPC"));Require(readback&&readback->GetNbinsX()==m2->GetNbinsX()&&readback->GetNbinsY()==m2->GetNbinsY()&&readback->GetEntries()==m2->GetEntries(),"Native M2 readback differs");double normal=0,flow=0;
    for(int x=0;x<=m2->GetNbinsX()+1;++x)for(int y=0;y<=m2->GetNbinsY()+1;++y){Require(readback->GetBinContent(x,y)==m2->GetBinContent(x,y)&&readback->GetBinError(x,y)==m2->GetBinError(x,y),"Native M2 count/error changed");if(x>0&&x<=m2->GetNbinsX()&&y>0&&y<=m2->GetNbinsY())normal+=m2->GetBinContent(x,y);else flow+=m2->GetBinContent(x,y);}check.Close();
    std::ofstream audit((output+"/tof_m2_audit.json").Data());audit<<std::setprecision(17)<<"{\"singleOriginalKeyRead\":true,\"histogram\":\""<<m2Key<<"\",\"entries\":"<<m2->GetEntries()<<",\"regularCount\":"<<normal<<",\"flowCount\":"<<flow<<",\"xBins\":"<<m2->GetNbinsX()<<",\"xMin\":"<<m2->GetXaxis()->GetXmin()<<",\"xMax\":"<<m2->GetXaxis()->GetXmax()<<",\"yBins\":"<<m2->GetNbinsY()<<",\"yMin\":"<<m2->GetYaxis()->GetXmin()<<",\"yMax\":"<<m2->GetYaxis()->GetXmax()<<",\"countErrorFlowReadbackExact\":true,\"noFit\":true,\"noDedxM2Joint\":true}\n";audit.close();Require(audit.good(),"Cannot finish native m2 audit");
    TFile native(cache,"READ");TH2* all=dynamic_cast<TH2*>(native.Get("ALL"));Require(!native.IsZombie()&&all,"Native ALL map absent");gStyle->SetOptStat(0);gStyle->SetPalette(1);TCanvas canvas("pid_parameter_compare","PID selection comparison",int(N(cfg,"canvasWidth")),int(N(cfg,"canvasHeight")));canvas.Print(pdf+"[");canvas.Clear();
    Text(.08,.94,"Configured PID comparison: 2021 P24iy, Au+Au FXT 3.85 GeV",.035);Text(.08,.87,"Old canonical He4 calibrated response = new common calibrated response",.028);Text(.08,.81,"Old dedicated d/t/He3 entries used legacy tables, including the competing species",.028);Text(.08,.75,"These are selection curves and gate edges, not a new response fit or confidence band",.026);Text(.08,.66,"TOF: same saved d/t/He3 mean +/- 1.5 sigma; no new TOF mean or width fit",.026);
    Text(.08,.58,"Species       saved / old reference       current configured         studied m2 window",.025);
    for(unsigned s=0;s<windows.size();++s)Text(.08,.51-.07*s,Form("%-5s       %.6f - %.6f       %.6f - %.6f       %.6f - %.6f",Species[s],N(windows[s],"oldTofExtendedMin"),N(windows[s],"oldTofExtendedMax"),N(windows[s],"currentConfiguredMin"),N(windows[s],"currentConfiguredMax"),N(windows[s],"studiedMin"),N(windows[s],"studiedMax")),.024);
    Text(.08,.18,"Old canonical/legacy CF selection: TPC only; old TOF windows were QA tags",.026);Text(.08,.12,"He4 old tofExtended/configured [3,4] is NOT the narrower D-window study",.026);Text(.08,.06,"Only input TOF QA map is added; original ROOT, PID coefficients and old products unchanged",.023);canvas.Print(pdf);canvas.Print(output+"/configuration_table.png");
    for(int s=0;s<4;++s)ResponsePage(canvas,curves,s,N(cfg,"rigidityMin"),N(cfg,"rigidityMax"),cfg,pdf,output,"full");for(int s=2;s<4;++s)ResponsePage(canvas,curves,s,N(cfg,"rigidityMin"),N(cfg,"lowRigidityMax"),cfg,pdf,output,"lowR");Maps(canvas,*all,*m2,curves,windows,cfg,pdf,output);canvas.Print(pdf+"]");native.Close();delete m2;std::cout<<"PID_SELECTION_COMPARISON_PASS"<<std::endl;
  }catch(const std::exception& e){std::cerr<<"ComparePIDSelections: "<<e.what()<<std::endl;gSystem->Exit(2);}
#endif
}
