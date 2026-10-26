// Saved-prediction comparison only: no fit, PID changes, or isotope assertion.
void CompareNuclearResponses(const char*,const char*,const char*);
#ifndef __CINT__
#include "TFile.h"
#include "TH1.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TMath.h"
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

namespace NuclearResponseCompare {
typedef std::map<std::string,std::string> Fields;
void Need(bool ok,const std::string& message){if(!ok)throw std::runtime_error(message);}
bool Finite(double v){return TMath::Finite(v);}
std::vector<std::string> Split(const std::string& line,char separator='\t'){
  std::vector<std::string> values;std::istringstream input(line);std::string value;
  while(std::getline(input,value,separator))values.push_back(value);return values;
}
double Num(const Fields& fields,const std::string& name){
  Fields::const_iterator found=fields.find(name);Need(found!=fields.end(),"Missing comparison field "+name);
  std::istringstream input(found->second);double value;std::string extra;
  Need(bool(input>>value)&&!(input>>extra)&&Finite(value),"Nonfinite/nonnumeric field "+name);return value;
}
std::string Text(const Fields& fields,const std::string& name){
  Fields::const_iterator found=fields.find(name);Need(found!=fields.end(),"Missing comparison field "+name);return found->second;
}
Fields Settings(const char* path){std::ifstream input(path);Need(bool(input),"Missing explicit settings");Fields cfg;std::string line;
  while(std::getline(input,line)){std::vector<std::string> item=Split(line);Need(item.size()==2&&cfg.insert(std::make_pair(item[0],item[1])).second,"Duplicate/bad settings row");}
  Need(input.eof(),"Settings read error");return cfg;
}
std::vector<Fields> Rows(const char* path){std::ifstream input(path);Need(bool(input),"Missing comparison rows");std::string line;
  Need(bool(std::getline(input,line)),"Missing rows header");std::vector<std::string> names=Split(line);std::vector<Fields> rows;
  while(std::getline(input,line)){std::vector<std::string> values=Split(line);Need(values.size()==names.size(),"Malformed comparison row");Fields row;
    for(unsigned i=0;i<names.size();++i)Need(row.insert(std::make_pair(names[i],values[i])).second,"Duplicate TSV column");rows.push_back(row);}
  Need(input.eof()&&!rows.empty(),"Empty or unreadable comparison rows");return rows;
}
std::vector<double> Numbers(const std::string& text){std::vector<std::string> tokens=Split(text,',');std::vector<double> result;
  for(unsigned i=0;i<tokens.size();++i){Fields f;f["value"]=tokens[i];result.push_back(Num(f,"value"));}return result;
}
TH1* Histogram(TFile& file,const char* name){TH1* h=dynamic_cast<TH1*>(file.Get(name));Need(h&&h->GetDimension()==1,std::string("Missing 1D histogram ")+name);return h;}
bool Close(double a,double b,double tolerance){return std::fabs(a-b)<=tolerance*std::max(1.,std::max(std::fabs(a),std::fabs(b)));}
void CheckAxis(const TH1& a,const TH1& b,double tolerance){
  Need(a.GetNbinsX()==b.GetNbinsX(),"Different native bin counts");
  for(int bin=1;bin<=a.GetNbinsX()+1;++bin)Need(Close(a.GetBinLowEdge(bin),b.GetBinLowEdge(bin),tolerance),"Different native axis edges");
}
void Audit(TFile& gaussian,TFile& lognormal,const Fields& row,double tolerance){
  int gi=int(Num(row,"gaussianIndex")),li=int(Num(row,"lognormalIndex"));
  TH1* gr=Histogram(gaussian,Form("raw_native_%03d",gi));TH1* lr=Histogram(lognormal,Form("raw_native_%03d",li));
  TH1* gm=Histogram(gaussian,Form("prediction_fit_mask_%03d",gi));TH1* lm=Histogram(lognormal,Form("prediction_fit_mask_%03d",li));
  TH1* gp=Histogram(gaussian,Form("prediction_total_%03d",gi));TH1* lp=Histogram(lognormal,Form("prediction_total_%03d",li));
  CheckAxis(*gr,*lr,tolerance);CheckAxis(*gr,*gm,tolerance);CheckAxis(*gr,*lm,tolerance);CheckAxis(*gr,*gp,tolerance);CheckAxis(*gr,*lp,tolerance);
  double observed=0;std::vector<std::string> species=Text(row,"sample")=="TAG_D"?Split("d,t,he3,he4",','):Split(Text(row,"sample")=="TAG_T"?"t":"he3",',');
  std::vector<TH1*> gc,lc;TH1* gs=Histogram(gaussian,Form("prediction_smooth_%03d",gi));TH1* ls=Histogram(lognormal,Form("prediction_smooth_%03d",li));
  CheckAxis(*gr,*gs,tolerance);CheckAxis(*gr,*ls,tolerance);
  for(unsigned j=0;j<species.size();++j){gc.push_back(Histogram(gaussian,Form("prediction_%s_%03d",species[j].c_str(),gi)));lc.push_back(Histogram(lognormal,Form("prediction_%s_%03d",species[j].c_str(),li)));CheckAxis(*gr,*gc.back(),tolerance);CheckAxis(*gr,*lc.back(),tolerance);}
  for(int bin=0;bin<=gr->GetNbinsX()+1;++bin){
    Need(gr->GetBinContent(bin)==lr->GetBinContent(bin)&&gr->GetBinError(bin)==lr->GetBinError(bin),"Native raw counts/errors/flows differ between responses");
    Need(gm->GetBinContent(bin)==lm->GetBinContent(bin),"Different saved native fit masks");
    Need(Finite(gr->GetBinContent(bin))&&gr->GetBinContent(bin)>=0,"Invalid observed count");
    Need(Finite(gp->GetBinContent(bin))&&gp->GetBinContent(bin)>=0&&Finite(lp->GetBinContent(bin))&&lp->GetBinContent(bin)>=0,"Invalid saved prediction");
    if(gm->GetBinContent(bin)>0)observed+=gr->GetBinContent(bin);
    double ga=gs->GetBinContent(bin),la=ls->GetBinContent(bin);for(unsigned j=0;j<gc.size();++j){ga+=gc[j]->GetBinContent(bin);la+=lc[j]->GetBinContent(bin);}
    Need(Close(ga,gp->GetBinContent(bin),tolerance)&&Close(la,lp->GetBinContent(bin),tolerance),"Saved prediction component closure failed");
  }
  Need(Close(observed,Num(row,"fitCount"),tolerance),"Native mask count differs from JSON fitCount");
}
void Note(double y,const std::string& message,int color=kBlack,double size=.027){TLatex text;text.SetNDC();text.SetTextColor(color);text.SetTextSize(size);text.DrawLatex(.12,y,message.c_str());}
int Marker(int response,int background){return response?(background?25:21):(background?24:20);}
int Color(int response){return response?kOrange+7:kBlue+1;}
void PlotMetric(TPad& pad,const std::vector<Fields>& rows,const std::string& species,const std::string& metric,const Fields& cfg){
  pad.cd();double rlo=Num(cfg,"rigidityMin"),rhi=Num(cfg,"rigidityMax");
  double ymin=0,ymax=Num(cfg,metric=="center"?"centerDisplayMax":"widthDisplayMax");
  if(metric=="deviance"){ymin=1.e-6;ymax=10.;pad.SetLogy();}
  std::string title=metric=="center"?"Gaussian mean / log-normal mode (keV/cm)":(metric=="width"?"Gaussian std / log-normal curvature width (keV/cm)":"Saved fit deviance / fit-window count");
  TH1F* frame=pad.DrawFrame(rlo,ymin,rhi,ymax);frame->SetTitle((species+";Rigidity p/q (GeV/c);"+title).c_str());frame->GetYaxis()->SetTitleOffset(1.35);
  TLegend legend(.57,.70,.94,.91);legend.SetBorderSize(0);legend.SetFillStyle(0);legend.SetTextSize(.025);
  for(int response=0;response<2;++response)for(int background=0;background<2;++background){
    std::string prefix=response?"lognormal":"gaussian",model=background?"exponential":"linear";
    TGraphErrors points,bad;points.SetMarkerStyle(Marker(response,background));points.SetMarkerColor(Color(response));points.SetLineColor(Color(response));points.SetMarkerSize(.62);
    bad.SetMarkerStyle(Marker(response,background));bad.SetMarkerColor(kRed+1);bad.SetLineColor(kRed+1);bad.SetMarkerSize(.62);
    for(unsigned i=0;i<rows.size();++i){const Fields& row=rows[i];if(Text(row,"species")!=species||Text(row,"background")!=model)continue;
      double value=Num(row,prefix+(metric=="center"?"Center":(metric=="width"?"Width":"DPerCount")));
      double error=metric=="deviance"?0:std::max(0.,Num(row,prefix+(metric=="center"?"CenterError":"WidthError")));
      bool unsafe=Num(row,prefix+"NumericalValid")==0||Num(row,prefix+"ZeroYield")>0||Num(row,prefix+"OutsideDetector")>0||Num(row,prefix+"Boundary")>0||Num(row,prefix+"IdentityValidated")==0;
      if(metric=="deviance"&&value<=0)continue;
      bool clipped=value<ymin||value>ymax;double displayed=std::max(ymin*1.03,std::min(ymax*.985,value));TGraphErrors& graph=unsafe||clipped?bad:points;
      int n=graph.GetN();graph.SetPoint(n,Num(row,"r"),displayed);graph.SetPointError(n,0,clipped?0:std::min(error,(ymax-ymin)/10));
    }
    std::string label=(response?"Log-normal":"Gaussian")+std::string(background?" / exp BG":" / pol1 BG");
    TObject* goodDraw=points.DrawClone("P SAME");TObject* badDraw=bad.DrawClone("P SAME");legend.AddEntry(points.GetN()?goodDraw:badDraw,label.c_str(),"p");
  }
  legend.DrawClone();Note(.65,"Red: unvalidated identity and/or numerical/zero/axis/boundary warning",kRed+1,.020);
  Note(.61,"Points only; unsupported rows are never joined. Clipped values remain in TSV.",kBlack,.020);
}
TH1* Display(const TH1& source,int factor,const char* name){TH1* h=static_cast<TH1*>(source.Clone(name));h->SetDirectory(0);double before=h->Integral(0,h->GetNbinsX()+1);h->Rebin(factor);Need(Close(before,h->Integral(0,h->GetNbinsX()+1),1.e-13),"Display regrouping changed count/flow sum");return h;}
void RawPanel(TPad& pad,TH1& raw,TH1& gaussian,TH1& lognormal,double lo,double hi,bool log,const Fields& row,const std::string& subtitle){
  pad.cd();pad.SetLogy(log);raw.GetXaxis()->SetRangeUser(lo,hi);gaussian.GetXaxis()->SetRangeUser(lo,hi);lognormal.GetXaxis()->SetRangeUser(lo,hi);
  raw.SetTitle((Text(row,"species")+"-compatible TAG data, R="+Text(row,"r")+" / "+Text(row,"background")+" BG / "+subtitle+";dE/dx (keV/cm);Stored track counts / display bin").c_str());
  raw.SetLineColor(kBlack);raw.SetMarkerColor(kBlack);raw.SetMarkerStyle(20);raw.SetMarkerSize(.38);
  gaussian.SetLineColor(kBlue+1);gaussian.SetLineWidth(2);gaussian.SetLineStyle(1);lognormal.SetLineColor(kOrange+7);lognormal.SetLineWidth(2);lognormal.SetLineStyle(2);
  double maximum=std::max(raw.GetMaximum(),std::max(gaussian.GetMaximum(),lognormal.GetMaximum()));raw.SetMaximum(maximum*(log?15:1.3));raw.SetMinimum(log?.35:0);
  raw.DrawCopy("E");gaussian.DrawCopy("HIST SAME");lognormal.DrawCopy("HIST SAME");raw.DrawCopy("E SAME");
  TLegend legend(.58,.71,.93,.90);legend.SetBorderSize(0);legend.SetFillStyle(0);legend.SetTextSize(.025);legend.AddEntry(&raw,"Same saved native data","lep");legend.AddEntry(&gaussian,"Gaussian saved expectation","l");legend.AddEntry(&lognormal,"Log-normal saved expectation","l");legend.DrawClone();
  Note(.66,"Predictions are saved fit-window values, not extrapolations.",kBlack,.022);
  Note(.62,"Nominal component identity unvalidated; TAG_D is d/4He compatible.",kRed+1,.021);
}
void RawPage(TCanvas& canvas,TFile& gaussian,TFile& lognormal,TFile& displayFile,const Fields& row,const Fields& cfg,const std::string& pdf){
  int gi=int(Num(row,"gaussianIndex")),li=int(Num(row,"lognormalIndex")),factor=int(Num(cfg,"displayRebin"));
  TH1* original=Histogram(gaussian,Form("raw_native_%03d",gi));TH1* gp=Histogram(gaussian,Form("prediction_total_%03d",gi));TH1* lp=Histogram(lognormal,Form("prediction_total_%03d",li));
  Need(original->GetNbinsX()%factor==0,"Display factor does not partition native axis");
  TH1* raw=Display(*original,factor,"compare_raw");TH1* g=Display(*gp,factor,"compare_gauss");TH1* l=Display(*lp,factor,"compare_log");
  double lo=Num(cfg,"axisMin"),hi=Num(cfg,"axisMax"),roiLo=hi,roiHi=lo,nsigma=Num(cfg,"roiHalfWidth");
  const char* prefixes[]={"gaussian","lognormal"};for(int m=0;m<2;++m){std::string prefix=prefixes[m];double center=Num(row,prefix+"Center"),width=Num(row,prefix+"Width");
    if(Num(row,prefix+"ZeroYield")==0&&center>lo&&center<hi&&width>0){roiLo=std::min(roiLo,center-nsigma*width);roiHi=std::max(roiHi,center+nsigma*width);}}
  if(!(roiHi>roiLo)){roiLo=Num(row,"fitLo");roiHi=Num(row,"fitHi");}roiLo=std::max(lo,roiLo);roiHi=std::min(hi,roiHi);
  canvas.Clear();canvas.Divide(2,2);RawPanel(*static_cast<TPad*>(canvas.cd(1)),*raw,*g,*l,lo,hi,true,row,"full axis, log");
  RawPanel(*static_cast<TPad*>(canvas.cd(2)),*raw,*g,*l,lo,hi,false,row,"full axis, linear");
  RawPanel(*static_cast<TPad*>(canvas.cd(3)),*raw,*g,*l,roiLo,roiHi,false,row,"primary union ROI, linear");
  RawPanel(*static_cast<TPad*>(canvas.cd(4)),*raw,*g,*l,roiLo,roiHi,true,row,"primary union ROI, log");
  canvas.Print(pdf.c_str());displayFile.cd();std::string label=Text(row,"background")+"_"+Text(row,"sample")+"_"+Text(row,"species")+"_"+Text(row,"rBin");
  original->Write((label+"_observed_native_unchanged").c_str());raw->Write((label+"_observed_display_grouped").c_str());g->Write((label+"_gaussian_display_grouped").c_str());l->Write((label+"_lognormal_display_grouped").c_str());
  delete raw;delete g;delete l;
}
}
#endif
void CompareNuclearResponses(const char* output,const char* settings,const char* rowsTsv){
#ifndef __CINT__
  using namespace NuclearResponseCompare;
  try{
    Fields cfg=Settings(settings);std::vector<Fields> rows=Rows(rowsTsv);std::string folder(output),summary=folder+"/response_parameters.pdf",overlay=folder+"/saved_spectra_overlays.pdf",root=folder+"/comparison_display_inputs.root";
    Need(!gSystem->AccessPathName((folder+"/context.json").c_str()),"New comparison context required");
    Need(gSystem->AccessPathName(summary.c_str())&&gSystem->AccessPathName(overlay.c_str())&&gSystem->AccessPathName(root.c_str())&&gSystem->AccessPathName((folder+"/native_comparison_audit.json").c_str()),"Existing products will not be overwritten");
    TFile gl(Text(cfg,"gaussianRoot_linear").c_str(),"READ"),ge(Text(cfg,"gaussianRoot_exponential").c_str(),"READ"),ll(Text(cfg,"lognormalRoot_linear").c_str(),"READ"),le(Text(cfg,"lognormalRoot_exponential").c_str(),"READ");
    TFile* sources[]={&gl,&ge,&ll,&le};for(int i=0;i<4;++i)Need(!sources[i]->IsZombie()&&!sources[i]->TestBit(TFile::kRecovered),"Source ROOT unreadable/recovered");
    double tolerance=Num(cfg,"auditRelativeTolerance");std::set<std::string> checked;
    for(unsigned i=0;i<rows.size();++i){Fields row=rows[i];std::string key=Text(row,"background")+":"+Text(row,"sample")+":"+Text(row,"rBin");if(checked.insert(key).second)Audit(Text(row,"background")=="linear"?gl:ge,Text(row,"background")=="linear"?ll:le,row,tolerance);}
    gStyle->SetOptStat(0);gStyle->SetTitleFontSize(.034);gStyle->SetPaperSize(26,20);gStyle->SetPadLeftMargin(.14);gStyle->SetPadRightMargin(.05);gStyle->SetPadBottomMargin(.12);TH1::AddDirectory(kFALSE);
    TCanvas canvas("response_comparison","Saved Gaussian/log-normal responses",1500,1100);canvas.Print((summary+"[").c_str());const char* species[]={"d","t","he3","he4"};
    for(int s=0;s<4;++s){canvas.Clear();canvas.Divide(2,2);PlotMetric(*static_cast<TPad*>(canvas.cd(1)),rows,species[s],"center",cfg);PlotMetric(*static_cast<TPad*>(canvas.cd(2)),rows,species[s],"width",cfg);PlotMetric(*static_cast<TPad*>(canvas.cd(3)),rows,species[s],"deviance",cfg);
      canvas.cd(4);Note(.89,std::string("Nominal primary: ")+species[s],kBlack,.045);Note(.80,"Gaussian: center = arithmetic mean; width = std deviation.",kBlue+1,.030);Note(.73,"Log-normal: center = mode; width = local curvature scale.",kOrange+7,.030);Note(.66,"Widths have different meanings: this is NOT variance-unfolding.",kRed+1,.027);Note(.59,"Native-R projected responses; neither width is intrinsic validation.",kBlack,.027);Note(.52,"D/N uses the whole SAME fitted TAG window, not primary purity.",kBlack,.027);Note(.45,"For d and 4He, D/N repeats the SAME TAG_D fit (not independent).",kRed+1,.026);Note(.38,"Gaussian/log-normal are non-nested; delta D is not significance.",kBlack,.027);Note(.31,"No accepted isotope identities, production curves, or purity.",kRed+1,.028);Note(.24,"Four markers: filled/open circle = Gaussian pol1/exp;",kBlack,.027);Note(.18,"filled/open square = log-normal pol1/exp. See TSV for flags.",kBlack,.027);canvas.Print(summary.c_str());}
    canvas.Print((summary+"]").c_str());TFile display(root.c_str(),"CREATE");Need(!display.IsZombie(),"Cannot CREATE comparison-only ROOT");
    std::vector<double> representatives=Numbers(Text(cfg,"representativeRigidities"));canvas.Print((overlay+"[").c_str());unsigned pages=0;
    for(unsigned j=0;j<representatives.size();++j)for(int b=0;b<2;++b)for(int s=0;s<4;++s){const Fields* chosen=0;double distance=1.e100;
      for(unsigned i=0;i<rows.size();++i)if(Text(rows[i],"species")==species[s]&&Text(rows[i],"background")==std::string(b?"exponential":"linear")){double delta=std::fabs(Num(rows[i],"r")-representatives[j]);if(delta<distance){distance=delta;chosen=&rows[i];}}
      Need(chosen&&distance<Num(cfg,"representativeTolerance"),"Representative native-R point not covered: no silent nearest-bin substitution");RawPage(canvas,b?ge:gl,b?le:ll,display,*chosen,cfg,overlay);++pages;}
    canvas.Print((overlay+"]").c_str());display.Close();gl.Close();ge.Close();ll.Close();le.Close();
    TFile verify(root.c_str(),"READ");Need(!verify.IsZombie()&&!verify.TestBit(TFile::kRecovered)&&verify.GetListOfKeys()->GetSize()==int(pages*4),"Closed comparison ROOT readback/key count failed");verify.Close();
    std::ofstream audit((folder+"/native_comparison_audit.json").c_str());Need(bool(audit),"Cannot create audit JSON");audit<<"{\"noFit\":true,\"productionEnabled\":false,\"nativeCountsErrorsFlowsIdentical\":true,\"nativeAxesFitMasksIdentical\":true,\"predictionComponentClosure\":true,\"closedDisplayRootReadable\":true,\"checkedUniqueFitPairs\":"<<checked.size()<<",\"summaryPages\":4,\"overlayPages\":"<<pages<<",\"displayGroupingOnly\":true,\"nonNestedDeltaDIsNotSignificance\":true,\"identitiesValidated\":false}\n";audit.close();std::cout<<"NUCLEAR_RESPONSE_COMPARISON_COMPLETE rows="<<rows.size()<<" fitPairs="<<checked.size()<<" pages="<<pages<<std::endl;
  }catch(const std::exception& e){std::cerr<<"CompareNuclearResponses: "<<e.what()<<std::endl;gSystem->Exit(2);}
#endif
}
