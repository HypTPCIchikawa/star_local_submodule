// Saved fixed-PID peak-anchored diagnostics only: no fitting or normalization.
void PlotPeakAnchored(const char*,const char*,const char*,const char*,const char*);
#include "PlotCutoffSpectra.C"
#ifndef __CINT__
#include "TLine.h"
namespace PeakAnchoredViews {
using CutoffSpectra::Require;
struct Row {
  int index; bool numerical;
  std::string label,response,jointChi2,jointNdf,jointRatio,heChi2,heBins,hePull,tChi2,tBins,tPull,reason;
  double fitLo,fitHi,heLo,heHi,tLo,tHi,muT,sigmaT,muHe3,sigmaHe3,muHe4,sigmaHe4;
};
typedef std::map<std::string,Row> Rows;
std::vector<std::string> Fields(const std::string& text) {
  std::vector<std::string> out;size_t start=0;
  while(true){const size_t end=text.find('\t',start);out.push_back(text.substr(start,end==std::string::npos?end:end-start));if(end==std::string::npos)break;start=end+1;}return out;
}
double Number(const std::string& text) {
  double value;std::istringstream in(text);Require(bool(in>>value)&&std::isfinite(value),"Invalid peak-view geometry");in>>std::ws;Require(in.eof(),"Malformed peak-view geometry");return value;
}
std::string Key(const std::string& label,int index) {std::ostringstream s;s<<label<<'/'<<index;return s.str();}
Rows Read(const char* path) {
  std::ifstream in(path);Require(bool(in),"Missing explicit peak metrics TSV");std::string line;
  Require(bool(std::getline(in,line)),"Empty peak metrics TSV");const std::vector<std::string> header=Fields(line);
  std::map<std::string,unsigned> c;for(unsigned j=0;j<header.size();++j)Require(c.insert(std::make_pair(header[j],j)).second,"Repeated peak metrics column");
  const char* required[]={"label","index","response","numerical","fitLo","fitHi","heLo","heHi","tLo","tHi","muT","sigmaT","muHe3","sigmaHe3","muHe4","sigmaHe4","jointChi2","jointNdf","jointRatio","heChi2","heBins","hePull","tChi2","tBins","tPull","reason"};
  for(unsigned j=0;j<sizeof(required)/sizeof(required[0]);++j)Require(c.count(required[j])==1,"Missing peak metrics column");Rows out;
  while(std::getline(in,line)){
    const std::vector<std::string> f=Fields(line);Require(f.size()==header.size(),"Malformed peak metrics TSV row");Row r;
    const double index=Number(f[c["index"]]);Require(index>=0&&index==std::floor(index),"Invalid peak view index");r.index=int(index);
    r.label=f[c["label"]];r.response=f[c["response"]];Require(f[c["numerical"]]=="0"||f[c["numerical"]]=="1","Invalid saved numerical flag");r.numerical=f[c["numerical"]]=="1";
    r.fitLo=Number(f[c["fitLo"]]);r.fitHi=Number(f[c["fitHi"]]);r.heLo=Number(f[c["heLo"]]);r.heHi=Number(f[c["heHi"]]);r.tLo=Number(f[c["tLo"]]);r.tHi=Number(f[c["tHi"]]);
    r.muT=Number(f[c["muT"]]);r.sigmaT=Number(f[c["sigmaT"]]);r.muHe3=Number(f[c["muHe3"]]);r.sigmaHe3=Number(f[c["sigmaHe3"]]);r.muHe4=Number(f[c["muHe4"]]);r.sigmaHe4=Number(f[c["sigmaHe4"]]);
    r.jointChi2=f[c["jointChi2"]];r.jointNdf=f[c["jointNdf"]];r.jointRatio=f[c["jointRatio"]];r.heChi2=f[c["heChi2"]];r.heBins=f[c["heBins"]];r.hePull=f[c["hePull"]];r.tChi2=f[c["tChi2"]];r.tBins=f[c["tBins"]];r.tPull=f[c["tPull"]];r.reason=f[c["reason"]];
    Require(r.fitHi>r.fitLo&&r.heHi>r.heLo&&r.tHi>r.tLo&&r.sigmaT>0&&r.sigmaHe3>0&&r.sigmaHe4>0,"Invalid saved peak/view ranges");
    Require(out.insert(std::make_pair(Key(r.label,r.index),r)).second,"Repeated peak model/index metrics");
  }
  Require(in.eof()&&!out.empty(),"Unreadable peak metrics TSV");return out;
}
void Mark(double center,int color,const char* name,double low,double high,double ymin,double ymax) {
  if(center<low||center>high)return;
  TLine* line=new TLine(center,ymin,center,ymax);line->SetLineColor(color);line->SetLineStyle(3);line->SetLineWidth(2);line->Draw("SAME");
  std::cout<<"PEAK_ANCHORED_REFERENCE name="<<name<<" center="<<center<<" low="<<low<<" high="<<high<<std::endl;
}
void Annotate(const Row& r,bool helium,double low,double high) {
  gPad->SetTopMargin(.40);TIter next(gPad->GetListOfPrimitives());TObject* object=0;TH1* axis=0;
  while((object=next())){
    TH1* hist=dynamic_cast<TH1*>(object);if(hist&&!axis)axis=hist;
    TLegend* legend=dynamic_cast<TLegend*>(object);if(legend){legend->SetX1NDC(.13);legend->SetX2NDC(.96);legend->SetY1NDC(.603);legend->SetY2NDC(.645);legend->SetNColumns(4);legend->SetTextSize(.018);}
  }
  Require(axis!=0,"Missing raw axis for original-PID reference markers");
  const double ymin=axis->GetMinimum(),ymax=axis->GetMaximum();
  Mark(r.muT,kOrange+7,"t",low,high,ymin,ymax);Mark(r.muHe3,kGreen+2,"he3",low,high,ymin,ymax);Mark(r.muHe4,kMagenta+1,"he4",low,high,ymin,ymax);
  TLatex text;text.SetNDC();text.SetTextFont(42);text.SetTextSize(.024);
  text.DrawLatex(.13,.845,Form("Joint Pearson #chi^{2}/conditional ndf = %s/%s = %s",r.jointChi2.c_str(),r.jointNdf.c_str(),r.jointRatio.c_str()));
  text.DrawLatex(.13,.807,Form("%s ROI #chi^{2}=%s; bins=%s; max|pull|=%s (no regional ndf)",helium?"^{3}He+^{4}He":"t control",(helium?r.heChi2:r.tChi2).c_str(),(helium?r.heBins:r.tBins).c_str(),(helium?r.hePull:r.tPull).c_str()));
  text.DrawLatex(.13,.769,Form("Original PID: #mu_{t}=%.3g, #mu_{3}=%.3g, #mu_{4}=%.3g keV/cm",r.muT,r.muHe3,r.muHe4));
  text.DrawLatex(.13,.731,Form("Joint fit [%.2f, %.2f]; view [%.2f, %.2f]; fixed shapes",r.fitLo,r.fitHi,low,high));
  text.SetTextColor(r.numerical?kOrange+7:kRed+1);text.DrawLatex(.13,.693,Form("Numerical: %s; isotope support / purity: not certified",r.numerical?"pass":"FAIL"));
  text.SetTextColor(kBlack);text.SetTextSize(.020);
  const std::string reason=r.reason.empty()?"No fit failure; response-model validation remains required":r.reason;
  text.DrawLatex(.13,.660,reason.substr(0,125).c_str());gPad->Modified();gPad->Update();
}
}
#endif
void PlotPeakAnchored(const char* modelsPath,const char* slicesPath,const char* settingsPath,const char* metricsPath,const char* pdf) {
#ifndef __CINT__
  using namespace CutoffSpectra;
  try {
    Require(modelsPath&&slicesPath&&settingsPath&&metricsPath&&pdf,"Explicit peak diagnostic inputs required");Require(gSystem->AccessPathName(pdf),"Existing peak PDF will not be overwritten");
    const Settings cfg=ReadSettings(settingsPath);const double group=Number(cfg,"displayDedxRebinFactor"),floor=Number(cfg,"logMinimum");Require(group>=1&&group==std::floor(group)&&floor>0,"Invalid explicit peak display settings");
    std::vector<Model> models=ReadModels(modelsPath);const std::vector<Selection> selected=ReadSelections(slicesPath);const PeakAnchoredViews::Rows metrics=PeakAnchoredViews::Read(metricsPath);
    gStyle->SetOptStat(0);gStyle->SetPaperSize(26,20);TCanvas canvas("peak_anchored_views","Fixed PID peak-anchored ALL diagnostics",1400,1100);
    canvas.Print((std::string(pdf)+"[").c_str());int serial=0,pages=0;
    for(unsigned i=0;i<selected.size();++i) {
      const Selection& selection=selected[i];
      for(unsigned m=0;m<models.size();++m){VerifyRaw(*Load(*models[0].file,"raw_native",selection.index),*Load(*models[m].file,"raw_native",selection.index));Require(metrics.count(PeakAnchoredViews::Key(models[m].label,selection.index)),"Missing peak model/index metric row");}
      for(int role=0;role<2;++role)for(unsigned pair=0;pair<models.size();pair+=2) {
        const bool helium=role==0;canvas.Clear();canvas.Divide(2,2);
        for(unsigned m=0;m<2&&pair+m<models.size();++m)for(int log=0;log<2;++log) {
          const Model& model=models[pair+m];const PeakAnchoredViews::Row& row=metrics.find(PeakAnchoredViews::Key(model.label,selection.index))->second;
          const double low=helium?row.heLo:row.tLo,high=helium?row.heHi:row.tHi;canvas.cd(2*m+log+1);
          std::cout<<"PEAK_ANCHORED_PANEL serial="<<serial<<" index="<<selection.index<<" label="<<model.label<<" role="<<(helium?"helium":"tControl")<<" log="<<log<<" low="<<low<<" high="<<high<<std::endl;
          DrawPanel(*model.file,model.label,selection,int(group),low,high,log!=0,floor,serial++);PeakAnchoredViews::Annotate(row,helium,low,high);
        }
        canvas.Print(pdf);++pages;
      }
    }
    canvas.Print((std::string(pdf)+"]").c_str());for(unsigned m=0;m<models.size();++m){models[m].file->Close();delete models[m].file;}
    std::cout<<"PEAK_ANCHORED_VIEWS_COMPLETE pages="<<pages<<" panels="<<serial<<std::endl;
  }catch(const std::exception& error){std::cerr<<"PlotPeakAnchored: "<<error.what()<<std::endl;gSystem->Exit(2);}
#endif
}
