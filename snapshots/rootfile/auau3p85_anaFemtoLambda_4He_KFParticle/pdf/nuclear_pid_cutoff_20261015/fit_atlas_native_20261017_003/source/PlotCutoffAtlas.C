// Presentation-only atlas: reuse the established spectrum drawing without refitting.
void PlotCutoffAtlas(const char*,const char*,const char*,const char*,const char*);
#include "PlotCutoffSpectra.C"
#ifndef __CINT__
namespace CutoffAtlas {
struct Quality { bool fit; bool response; std::string species; };
typedef std::map<std::string,Quality> Audit;
std::vector<std::string> Fields(const std::string& line){
  std::vector<std::string> out;size_t begin=0;
  while(true){size_t pos=line.find('\t',begin);out.push_back(line.substr(begin,pos==std::string::npos?pos:pos-begin));if(pos==std::string::npos)break;begin=pos+1;}
  return out;
}
bool Flag(const std::string& value){CutoffSpectra::Require(value=="True"||value=="False"||value=="1"||value=="0","Malformed explicit audit flag");return value=="True"||value=="1";}
std::string LabelKey(const std::string& label){
  if(label=="Gaussian + linear")return "gaussian_linear_nominal";
  if(label=="Gaussian + exponential")return "gaussian_exponential_nominal";
  if(label=="Log-normal + linear")return "lognormal_linear_nominal";
  if(label=="Log-normal + exponential")return "lognormal_exponential_nominal";
  throw std::runtime_error("Unknown explicit saved model label: "+label);
}
std::string Key(const std::string& label,int index){std::ostringstream out;out<<label<<'/'<<index;return out.str();}
Audit ReadAudit(const char* path){
  using CutoffSpectra::Require;std::ifstream in(path);Require(bool(in),"Missing explicit atlas quality table");
  std::string line;Require(bool(std::getline(in,line)),"Empty atlas quality table");std::vector<std::string> head=Fields(line);std::map<std::string,unsigned> col;
  for(unsigned i=0;i<head.size();++i)Require(col.insert(std::make_pair(head[i],i)).second,"Duplicate atlas audit column");
  const char* keys[]={"model","index","numericalValid","fitResponseQuality","species"};for(unsigned i=0;i<5;++i)Require(col.count(keys[i])==1,"Missing atlas audit column");
  Audit out;while(std::getline(in,line)){std::vector<std::string> row=Fields(line);Require(row.size()==head.size(),"Malformed atlas audit row");std::istringstream number(row[col["index"]]);int index;Require(bool(number>>index)&&index>=0,"Invalid atlas audit index");Quality q;q.fit=Flag(row[col["numericalValid"]]);q.response=Flag(row[col["fitResponseQuality"]]);q.species=row[col["species"]];Require(q.species=="d"||q.species=="t"||q.species=="he3"||q.species=="he4","Unknown atlas species");Require(out.insert(std::make_pair(Key(row[col["model"]],index),q)).second,"Duplicate model/index audit row (supply one species only)");}
  Require(in.eof()&&!out.empty(),"Unreadable/empty atlas audit");return out;
}
void Annotate(const Quality& quality){
  // Keep the extra status in a dedicated header band above the plot/exponent.
  gPad->SetTopMargin(.23);gPad->Modified();gPad->Update();
  TLatex text;text.SetNDC();text.SetTextFont(42);text.SetTextSize(.026);text.SetTextColor(!quality.fit?kRed+1:(!quality.response?kOrange+7:kGreen+3));
  const char* species=quality.species=="he3"?"^{3}He":(quality.species=="he4"?"^{4}He":quality.species.c_str());
  text.DrawLatex(.13,.835,Form("target: %s; fit: %s; response: %s",species,quality.fit?"valid":"FAILED",quality.response?"pass":"FAIL"));
  gPad->Modified();gPad->Update();
}
}
#endif
void PlotCutoffAtlas(const char* modelsPath,const char* selectionsPath,const char* settingsPath,const char* auditPath,const char* outputPdf){
#ifndef __CINT__
  using namespace CutoffSpectra;
  try{
    Require(modelsPath&&selectionsPath&&settingsPath&&auditPath&&outputPdf,"Explicit atlas plot inputs required");Require(gSystem->AccessPathName(outputPdf),"Existing atlas PDF will not be overwritten");
    Settings cfg=ReadSettings(settingsPath);double low=Number(cfg,"dedxMin"),high=Number(cfg,"dedxMax"),logMin=Number(cfg,"logMinimum"),group=Number(cfg,"displayDedxRebinFactor");
    Require(high>low&&logMin>0&&group>=1&&group==std::floor(group),"Invalid explicit atlas view settings");int rebin=int(group);
    std::vector<Model> models=ReadModels(modelsPath);std::vector<Selection> selections=ReadSelections(selectionsPath);CutoffAtlas::Audit audit=CutoffAtlas::ReadAudit(auditPath);
    gStyle->SetOptStat(0);gStyle->SetPaperSize(26,20);TCanvas canvas("nuclear_cutoff_atlas","Saved nuclear fit atlas",1300,1000);canvas.Print((std::string(outputPdf)+"[").c_str());int serial=0;
    for(unsigned i=0;i<selections.size();++i){const Selection& s=selections[i];double fitLo=high,fitHi=low;
      for(unsigned m=0;m<models.size();++m){VerifyRaw(*Load(*models[0].file,"raw_native",s.index),*Load(*models[m].file,"raw_native",s.index));Require(audit.count(CutoffAtlas::Key(CutoffAtlas::LabelKey(models[m].label),s.index))==1,"Missing per-panel saved quality audit");TH1* mask=Load(*models[m].file,"prediction_fit_mask",s.index);for(int b=1;b<=mask->GetNbinsX();++b)if(mask->GetBinContent(b)>0){fitLo=std::min(fitLo,mask->GetBinLowEdge(b));fitHi=std::max(fitHi,mask->GetXaxis()->GetBinUpEdge(b));}}
      for(int view=0;view<2;++view)for(unsigned pair=0;pair<models.size();pair+=2){double a=view==0?low:std::max(low,fitLo),b=view==0?high:std::min(high,fitHi);if(!(b>a)){a=low;b=high;}canvas.Clear();canvas.Divide(2,2);
        for(unsigned m=0;m<2&&pair+m<models.size();++m)for(int logarithmic=0;logarithmic<2;++logarithmic){const Model& model=models[pair+m];canvas.cd(2*m+logarithmic+1);DrawPanel(*model.file,model.label,s,rebin,a,b,logarithmic!=0,logMin,serial++);CutoffAtlas::Annotate(audit.find(CutoffAtlas::Key(CutoffAtlas::LabelKey(model.label),s.index))->second);}
        canvas.Print(outputPdf);
      }
    }
    canvas.Print((std::string(outputPdf)+"]").c_str());for(unsigned m=0;m<models.size();++m){models[m].file->Close();delete models[m].file;}
    std::cout<<"NUCLEAR_CUTOFF_ATLAS_COMPLETE pages="<<2*((models.size()+1)/2)*selections.size()<<std::endl;
  }catch(const std::exception& e){std::cerr<<"PlotCutoffAtlas: "<<e.what()<<std::endl;gSystem->Exit(2);}
#endif
}
