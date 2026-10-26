#include <TFile.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TStyle.h>
#include <TColor.h>
#include <TEnv.h>
#include <TLatex.h>
#include <TNamed.h>
#include <TSystem.h>
#include <TMath.h>
#include <TString.h>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <stdexcept>
#include <cmath>

std::string BeRequired(TEnv& e,const char* key){const std::string value=e.GetValue(key,"");if(value.empty())throw std::runtime_error(std::string("Missing plot config: ")+key);return value;}
struct BeSample{std::string key,title,selection;TH2D* h;BeSample():h(0){}};
bool BeSameAxis(const TAxis* a,const TAxis* b){if(a->GetNbins()!=b->GetNbins())return false;for(int i=1;i<=a->GetNbins()+1;++i)if(a->GetBinLowEdge(i)!=b->GetBinLowEdge(i))return false;return true;}
void BeText(double x,double y,const char* text,double size=0.033){TLatex latex;latex.SetNDC();latex.SetTextFont(42);latex.SetTextSize(size);latex.DrawLatex(x,y,text);}
void BeDraw(BeSample& s,double etaMin,double etaMax,double zmax,const std::string& source,const std::string& beta){
 gPad->SetLeftMargin(.15);gPad->SetRightMargin(.20);gPad->SetBottomMargin(.15);gPad->SetTopMargin(.20);gPad->SetTicks(1,1);gPad->SetLogz();
 s.h->SetTitle("");s.h->GetXaxis()->SetRangeUser(etaMin,etaMax);s.h->GetYaxis()->SetRange(1,s.h->GetNbinsY());s.h->SetMinimum(.5);s.h->SetMaximum(zmax);
 s.h->GetXaxis()->SetTitle("#eta_{lab}");s.h->GetYaxis()->SetTitle("#beta = |p|/E");s.h->GetZaxis()->SetTitle("Candidates / bin");
 s.h->GetXaxis()->SetTitleSize(.052);s.h->GetYaxis()->SetTitleSize(.052);s.h->GetZaxis()->SetTitleSize(.043);
 s.h->GetXaxis()->SetLabelSize(.043);s.h->GetYaxis()->SetLabelSize(.043);s.h->GetZaxis()->SetLabelSize(.038);
 s.h->GetXaxis()->SetTitleOffset(1.1);s.h->GetYaxis()->SetTitleOffset(1.25);s.h->GetZaxis()->SetTitleOffset(1.5);
 s.h->Draw("COLZ");
 BeText(.15,.95,s.title.c_str(),.038);BeText(.15,.90,s.selection.c_str(),.027);
 const double flow=s.h->GetEntries()-s.h->Integral(1,s.h->GetNbinsX(),1,s.h->GetNbinsY());
 BeText(.15,.845,(flow>0?TString::Format("N = %.0f (flow: %.0f)",s.h->GetEntries(),flow):TString::Format("N = %.0f",s.h->GetEntries())).Data(),.030);
 BeText(.15,.037,source.c_str(),.023);BeText(.15,.014,beta.c_str(),.022);gPad->Modified();gPad->Update();
}
void BePage(std::vector<BeSample>& samples,const std::vector<int>& ids,int width,int height,double etaMin,double etaMax,double zmax,const std::string& pdf,const std::string& png,const std::string& source,const std::string& beta){
 TCanvas c("beta_eta_page","beta eta",width,height);if(ids.size()>1)c.Divide(ids.size(),1,0,0);
 for(unsigned i=0;i<ids.size();++i){if(ids.size()>1)c.cd(i+1);else c.cd();BeDraw(samples[ids[i]],etaMin,etaMax,zmax,source,beta);}
 c.Print(pdf.c_str());if(!png.empty())c.Print(png.c_str());
}
void PlotBetaEta(const char* configPath){
 try{
  TEnv e;if(e.ReadFile(configPath,kEnvLocal)!=0)throw std::runtime_error("Cannot read explicit plot config");
  const std::string inputList=BeRequired(e,"Inputs"),out=BeRequired(e,"Output"),source=BeRequired(e,"SourceCaption"),beta=BeRequired(e,"BetaCaption");
  const int n=e.GetValue("Samples",-1),margin=e.GetValue("EtaMarginBins",-1);if(n!=6||margin<0)throw std::runtime_error("Invalid explicit plotting settings");
  std::vector<BeSample> samples(n);
  for(int i=0;i<n;++i){samples[i].key=BeRequired(e,TString::Format("Sample.%d.Key",i));samples[i].title=BeRequired(e,TString::Format("Sample.%d.Title",i));samples[i].selection=BeRequired(e,TString::Format("Sample.%d.Selection",i));}
  TH1::AddDirectory(kFALSE);std::ifstream list(inputList.c_str());std::string path;std::vector<std::string> sources;
  double inputEntries=0,recoEntries=0;
  while(std::getline(list,path)){
   if(path.empty())continue;TFile f(path.c_str(),"READ");if(f.IsZombie()||f.TestBit(TFile::kRecovered))throw std::runtime_error("Bad input ROOT "+path);
   std::cout<<"SOURCE "<<path<<std::endl;
   TH1* input=dynamic_cast<TH1*>(f.Get("hRefMult"));TH1* reco=dynamic_cast<TH1*>(f.Get("hRefMultCorr"));if(!input||!reco)throw std::runtime_error("Missing event counter histograms");inputEntries+=input->GetEntries();recoEntries+=reco->GetEntries();delete input;delete reco;
   for(int i=0;i<n;++i){
    TH2D* h=dynamic_cast<TH2D*>(f.Get(samples[i].key.c_str()));if(!h)throw std::runtime_error("Missing TH2D "+samples[i].key);
    if(!samples[i].h){samples[i].h=(TH2D*)h->Clone(("sum_"+samples[i].key).c_str());samples[i].h->SetDirectory(0);if(!samples[i].h->GetSumw2N())samples[i].h->Sumw2();}
    else{if(!BeSameAxis(samples[i].h->GetXaxis(),h->GetXaxis())||!BeSameAxis(samples[i].h->GetYaxis(),h->GetYaxis()))throw std::runtime_error("Input axes differ "+samples[i].key);if(!samples[i].h->Add(h))throw std::runtime_error("Cannot add TH2D");}
    delete h;
   }
   sources.push_back(path);f.Close();
  }
  if(sources.size()!=3)throw std::runtime_error("Exactly three reviewed intermediate ROOTs required");
  double zmax=0;int occupiedLo=samples[0].h->GetNbinsX()+1,occupiedHi=0;
  std::ofstream counts((out+"/counts.tsv").c_str());counts<<std::setprecision(17)<<"key\tentries\tregular_bin_sum\tflow_bin_sum\tnonzero_regular_bins\tmaximum_bin\n";
  for(int i=0;i<n;++i){TH2D* h=samples[i].h;double all=0,regular=0;int occupied=0;
   if(h->GetEntries()<=0)throw std::runtime_error("Empty requested histogram "+samples[i].key);
   h->SetTitle(samples[i].title.c_str());h->GetXaxis()->SetTitle(BeRequired(e,"XTitle").c_str());h->GetYaxis()->SetTitle(BeRequired(e,"YTitle").c_str());
   for(int x=0;x<=h->GetNbinsX()+1;++x)for(int y=0;y<=h->GetNbinsY()+1;++y){double v=h->GetBinContent(x,y);if(!TMath::Finite(v)||v<0)throw std::runtime_error("Nonfinite/negative candidate count");all+=v;if(x>=1&&x<=h->GetNbinsX()&&y>=1&&y<=h->GetNbinsY()){regular+=v;if(v>0){++occupied;if(x<occupiedLo)occupiedLo=x;if(x>occupiedHi)occupiedHi=x;if(v>zmax)zmax=v;}}}
   if(all!=h->GetEntries())throw std::runtime_error("Unit-count entries vs all-bin sum mismatch");
   counts<<samples[i].key<<"\t"<<h->GetEntries()<<"\t"<<regular<<"\t"<<all-regular<<"\t"<<occupied<<"\t"<<h->GetMaximum()<<"\n";
  }
  int closure=0,subsets=0;
  for(int x=0;x<=samples[1].h->GetNbinsX()+1;++x)for(int y=0;y<=samples[1].h->GetNbinsY()+1;++y){
   const double low=samples[2].h->GetBinContent(x,y),high=samples[5].h->GetBinContent(x,y),match=samples[3].h->GetBinContent(x,y),m2=samples[4].h->GetBinContent(x,y);
   if(low+high!=samples[1].h->GetBinContent(x,y))++closure;if(m2>match||match>low)++subsets;
  }
  if(closure||subsets)throw std::runtime_error("He4 partition/subset closure failed");
  if(samples[2].h->GetEntries()+samples[5].h->GetEntries()!=samples[1].h->GetEntries())throw std::runtime_error("He4 Entries closure failed");
  counts.close();
  int lo=TMath::Max(1,occupiedLo-margin),hi=TMath::Min(samples[0].h->GetNbinsX(),occupiedHi+margin);
  const double etaMin=samples[0].h->GetXaxis()->GetBinLowEdge(lo),etaMax=samples[0].h->GetXaxis()->GetBinUpEdge(hi);zmax*=1.05;
  {TFile cache((out+"/beta_eta_summed.root").c_str(),"RECREATE");if(cache.IsZombie())throw std::runtime_error("Cannot create summed ROOT");for(int i=0;i<n;++i)if(samples[i].h->Write(samples[i].key.c_str())<=0)throw std::runtime_error("Cannot save merged TH2D");
   std::ostringstream consumed;std::ifstream envfile(configPath);consumed<<envfile.rdbuf();TNamed cfg("BetaEtaPlotConfiguration",consumed.str().c_str());cfg.Write();TNamed definition("BetaEtaDefinition","Kinematic beta=|p|/E, lab eta; raw candidate counts, no pair multiplicity weight; lowR TOF samples overlap");definition.Write();cache.Close();}
  gStyle->SetOptStat(0);gStyle->SetCanvasColor(0);gStyle->SetPadColor(0);gStyle->SetFrameFillColor(0);gStyle->SetTitleFont(42,"XYZ");gStyle->SetLabelFont(42,"XYZ");
  Double_t stops[5]={0,.25,.5,.75,1},red[5]={.267,.230,.128,.369,.993},green[5]={.005,.322,.567,.789,.906},blue[5]={.329,.546,.551,.383,.144};
  Int_t first=TColor::CreateGradientColorTable(5,stops,red,green,blue,255),palette[255];for(int i=0;i<255;++i)palette[i]=first+i;gStyle->SetPalette(255,palette);gStyle->SetNumberContours(255);
  const std::string pdf=out+"/beta_vs_eta_lambda_he4.pdf";
  {TCanvas opener;opener.Print((pdf+"[").c_str());}
  std::vector<int> main;main.push_back(0);main.push_back(2);main.push_back(5);
  BePage(samples,main,1800,700,etaMin,etaMax,zmax,pdf,out+"/overview.png",source,beta);
  std::vector<int> lambda(1,0);BePage(samples,lambda,1050,850,etaMin,etaMax,zmax,pdf,"",source,beta);
  std::vector<int> split;split.push_back(2);split.push_back(5);BePage(samples,split,1500,750,etaMin,etaMax,zmax,pdf,out+"/he4_low_high.png",source,beta);
  std::vector<int> low;low.push_back(2);low.push_back(3);low.push_back(4);BePage(samples,low,1800,700,etaMin,etaMax,zmax,pdf,out+"/he4_low_tof_comparison.png",source,beta);
  BePage(samples,main,1800,700,samples[0].h->GetXaxis()->GetXmin(),samples[0].h->GetXaxis()->GetXmax(),zmax,pdf,out+"/full_eta_range.png",source,beta);
  {TCanvas closer;closer.Print((pdf+"]").c_str());}
  for(int i=0;i<n;++i){TCanvas c("single_beta_eta","single",1050,850);c.cd();BeDraw(samples[i],etaMin,etaMax,zmax,source,beta);c.Print((out+"/"+samples[i].key+".png").c_str());}
  std::ofstream validation((out+"/validation.tsv").c_str());validation<<std::setprecision(17)<<"inputs\t"<<sources.size()<<"\ninputEvents\t"<<inputEntries<<"\nreconstructedEvents\t"<<recoEntries<<"\nhe4_partition_mismatch_bins\t"<<closure<<"\nhe4_subset_mismatch_bins\t"<<subsets<<"\neta_display_min\t"<<etaMin<<"\neta_display_max\t"<<etaMax<<"\nshared_zmax\t"<<zmax<<"\nDONE\t1\n";validation.close();
  for(int i=0;i<n;++i)delete samples[i].h;
  std::cout<<"BETA_ETA_DONE "<<out<<std::endl;gSystem->Exit(0);
 }catch(const std::exception& error){std::cerr<<"ERROR: "<<error.what()<<std::endl;gSystem->Exit(1);}
}
