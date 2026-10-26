// Validation plots of raw CF inputs, never computes a contamination subtraction.
#include "TFile.h"
#include "TH1D.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TText.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TPad.h"
#include "TString.h"
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <stdexcept>
#include <iostream>
#include <algorithm>
#include <cmath>

#ifndef __CINT__
namespace HeliumInputQa {
std::string Trim(const std::string& value){
  const size_t first=value.find_first_not_of(" \t\r");if(first==std::string::npos)return "";
  return value.substr(first,value.find_last_not_of(" \t\r")-first+1);
}
std::map<std::string,std::string> Flat(const std::string& path){
  std::ifstream file(path.c_str());if(!file)throw std::runtime_error("cannot read explicit QA configuration: "+path);
  std::map<std::string,std::string> result;std::string line;
  while(std::getline(file,line)){
    line=line.substr(0,line.find('#'));const size_t colon=line.find(':');if(colon==std::string::npos)continue;
    const std::string key=Trim(line.substr(0,colon));std::string value=Trim(line.substr(colon+1));
    if(value.size()>1&&((value[0]=='"'&&value[value.size()-1]=='"')||(value[0]=='\''&&value[value.size()-1]=='\'')))value=value.substr(1,value.size()-2);
    if(!result.insert(std::make_pair(key,value)).second)throw std::runtime_error("duplicate QA key: "+key);
  }return result;
}
std::string Required(const std::map<std::string,std::string>& values,const char* key){
  std::map<std::string,std::string>::const_iterator i=values.find(key);if(i==values.end()||i->second.empty())throw std::runtime_error(std::string("missing QA key: ")+key);return i->second;
}
int Integer(const std::map<std::string,std::string>& values,const char* key){
  std::istringstream number(Required(values,key));double value;number>>value;
  if(!number||!std::isfinite(value)||value!=std::floor(value))throw std::runtime_error("invalid QA integer");return int(value);
}
TH2* Need2(TFile& file,const TString& name){
  TH2* hist=dynamic_cast<TH2*>(file.Get(name));if(!hist)throw std::runtime_error(std::string("missing 2D CF input: ")+name.Data());return hist;
}
void Label(const char* title,int first,int last){
  const std::string stableTitle(title);
  TText* header=new TText(.15,.94,stableTitle.c_str());header->SetNDC();header->SetTextFont(42);header->SetTextSize(.035);header->Draw();
  std::ostringstream selected;selected<<"cent9 "<<first<<"-"<<last<<"; raw pair counts";
  TText* selection=new TText(.15,.9,selected.str().c_str());selection->SetNDC();selection->SetTextFont(42);selection->SetTextSize(.024);selection->Draw();
}
}
#endif

void checkHistFemtoLambdaHeliumControl(const char* input,const char* output,const char* mainconf,int rebinFactor){
#ifndef __CINT__
  using namespace HeliumInputQa;
  try {
    if(!input||!*input||!output||!*output||!mainconf||!*mainconf||rebinFactor<=0)throw std::runtime_error("explicit input, new output, mainconf and rebin factor required");
    if(!gSystem->AccessPathName(output))throw std::runtime_error("refusing overwrite of QA PDF");
    std::string mainPath=mainconf;if(mainPath[0]!='/')mainPath=std::string(gSystem->WorkingDirectory())+"/"+mainPath;
    const size_t configRoot=mainPath.rfind("/config/");if(configRoot==std::string::npos)throw std::runtime_error("mainconf must reside below config/");
    const std::map<std::string,std::string> main=Flat(mainPath);std::string makerPath=Required(main,"maker");
    if(makerPath[0]!='/')makerPath=mainPath.substr(0,configRoot)+"/config/"+makerPath;
    const std::map<std::string,std::string> maker=Flat(makerPath);
    if(Required(maker,"lambdaHe3ControlEnabled")!="true")throw std::runtime_error("QA requires explicit helium control profile");
    const int first=Integer(maker,"cfCent9Min"),last=Integer(maker,"cfCent9Max");
    if(first<0||last>8||last<first)throw std::runtime_error("invalid QA centrality selection");
    TFile* file=TFile::Open(input,"READ");if(!file||file->IsZombie())throw std::runtime_error("cannot read CF-input ROOT");
    const char* species[]={"he4","he3","he3_as_he4"};const char* legacy[]={"4He","3He","3HeAs4He"};
    const char* titles[]={"Lambda-He4","Lambda-He3","He3-as-He4 (same selected pairs)"};
    const char* regions[]={"_signal","_leftSB","_rightSB"};const char* legendNames[]={"#Lambda signal #pm3#sigma","Left #Lambda sideband","Right #Lambda sideband"};
    const int colors[]={kBlack,kBlue+1,kRed+1};
    gStyle->SetOptStat(0);gStyle->SetPaperSize(30,20);
    TCanvas canvas("helium_cf_input_qa","Helium raw CF inputs",1600,1050);canvas.Print(TString(output)+"[");canvas.Divide(3,2);
    for(int mode=0;mode<2;++mode)for(int target=0;target<3;++target){
      canvas.cd(1+target+3*mode);gPad->SetLeftMargin(.15);gPad->SetTopMargin(.23);gPad->SetBottomMargin(.13);
      TH1D* projected[3]={0,0,0};double maximum=1.;
      for(int region=0;region<3;++region){
        TH2* source=Need2(*file,Form("hKstar%sVsCent_lambda_%s%s",mode?"ME":"SE",species[target],regions[region]));
        if(source->GetNbinsX()%rebinFactor)throw std::runtime_error("QA rebin factor does not divide saved kstar bins");
        projected[region]=source->ProjectionX(Form("qa_%d_%d_%d",mode,target,region),source->GetYaxis()->FindBin(first),source->GetYaxis()->FindBin(last),"e");projected[region]->SetDirectory(0);projected[region]->Rebin(rebinFactor);
        maximum=std::max(maximum,projected[region]->GetMaximum());
      }
      TLegend* legend=new TLegend(.5,.775,.965,.91);legend->SetBorderSize(0);legend->SetFillStyle(0);legend->SetTextSize(.023);
      for(int region=0;region<3;++region){
        TH1D* hist=projected[region];hist->SetTitle(Form(";k* (GeV/c);pairs / %.0f MeV/c",1000*hist->GetBinWidth(1)));
        hist->SetLineColor(colors[region]);hist->SetLineWidth(2);hist->SetMinimum(0.);hist->SetMaximum(1.3*maximum);hist->Draw(region?"hist same":"hist");
        legend->AddEntry(hist,Form("%s (%.0f)",legendNames[region],hist->Integral(0,hist->GetNbinsX()+1)),"l");
      }
      const std::string title=std::string(titles[target])+" | "+(mode?"ME":"SE");
      Label(title.c_str(),first,last);legend->Draw();
    }
    canvas.Print(output);canvas.Clear();canvas.Divide(3,2);
    for(int mode=0;mode<2;++mode)for(int target=0;target<3;++target){
      canvas.cd(1+target+3*mode);gPad->SetLeftMargin(.15);gPad->SetRightMargin(.13);gPad->SetTopMargin(.17);gPad->SetBottomMargin(.13);
      TH2* combined=0;
      for(int cent=first;cent<=last;++cent){
        TH2* source=Need2(*file,Form("%s/hKstarMt_%s%s_CentBin%d",mode?"mix":"true",mode?"Mixed_":"",legacy[target],cent));
        if(!combined){combined=static_cast<TH2*>(source->Clone(Form("mt_%d_%d",mode,target)));combined->SetDirectory(0);}else combined->Add(source);
      }
      int lower=combined->GetNbinsY(),upper=1;bool any=false;
      for(int x=1;x<=combined->GetNbinsX();++x)for(int y=1;y<=combined->GetNbinsY();++y)if(combined->GetBinContent(x,y)>0){lower=std::min(lower,y);upper=std::max(upper,y);any=true;}
      if(any)combined->GetYaxis()->SetRange(std::max(1,lower-1),std::min(combined->GetNbinsY(),upper+1));
      combined->SetTitle(";k* (GeV/c);m_{T} (GeV/c^{2});raw pairs / saved bin");combined->Draw("COLZ");
      const std::string title=std::string(titles[target])+" | "+(mode?"ME":"SE")+" signal";
      Label(title.c_str(),first,last);
    }
    canvas.Print(output);canvas.Print(TString(output)+"]");file->Close();delete file;
    std::cout<<"HELIUM_CF_INPUT_QA_COMPLETE pages=2 mainconf="<<mainPath<<" rebin="<<rebinFactor<<" (raw inputs, no CF subtraction)"<<std::endl;gSystem->Exit(0);
  }catch(const std::exception& error){std::cerr<<"HELIUM_CF_INPUT_QA_FAILED: "<<error.what()<<std::endl;gSystem->Exit(1);}
#endif
}
