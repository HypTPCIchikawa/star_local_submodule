// Presentation only: reuse saved counts and predictions without fitting/scaling.
void PlotQualityDiagnostics(const char*,const char*,const char*,const char*,const char*);
#include "PlotCutoffSpectra.C"
#ifndef __CINT__
namespace QualityDiagnostics {
struct Row {
  int index;
  std::string label,target,chi2,ndf,ratio,reason;
  bool numerical,supported,supportKnown;
  double fitLo,fitHi,focusLo,focusHi;
};
typedef std::map<std::string,Row> Rows;
std::vector<std::string> Fields(const std::string& s){
  std::vector<std::string> v;size_t a=0;
  while(true){size_t b=s.find('\t',a);v.push_back(s.substr(a,b==std::string::npos?b:b-a));if(b==std::string::npos)break;a=b+1;}
  return v;
}
std::string Key(const std::string& label,int i){std::ostringstream s;s<<label<<'/'<<i;return s.str();}
double Num(const std::string& s){
  std::istringstream in(s);double v;
  CutoffSpectra::Require(bool(in>>v)&&std::isfinite(v),"Invalid view geometry");return v;
}
bool Flag(const std::string& s){CutoffSpectra::Require(s=="0"||s=="1","Invalid view flag");return s=="1";}
Rows Read(const char* path){
  using CutoffSpectra::Require;
  std::ifstream in(path);Require(bool(in),"Missing quality TSV");std::string s;
  Require(bool(std::getline(in,s)),"Empty quality TSV");std::vector<std::string> h=Fields(s);
  std::map<std::string,unsigned> c;
  for(unsigned i=0;i<h.size();++i)Require(c.insert(std::make_pair(h[i],i)).second,"Repeated quality column");
  const char* k[]={"label","index","target","chi2","ndf","chi2Ndf","numerical","supported","reason","fitLo","fitHi","focusLo","focusHi"};
  for(unsigned i=0;i<13;++i)Require(c.count(k[i])==1,"Missing quality column");
  Rows out;
  while(std::getline(in,s)){
    std::vector<std::string> f=Fields(s);Require(f.size()==h.size(),"Malformed quality row");Row r;
    r.label=f[c["label"]];double index=Num(f[c["index"]]);
    Require(index>=0&&index==std::floor(index),"Invalid saved index");r.index=int(index);
    r.target=f[c["target"]];r.chi2=f[c["chi2"]];r.ndf=f[c["ndf"]];r.ratio=f[c["chi2Ndf"]];
    r.numerical=Flag(f[c["numerical"]]);r.supported=Flag(f[c["supported"]]);r.reason=f[c["reason"]];
    r.supportKnown=c.count("supportKnown")?Flag(f[c["supportKnown"]]):true;
    r.fitLo=Num(f[c["fitLo"]]);r.fitHi=Num(f[c["fitHi"]]);
    r.focusLo=Num(f[c["focusLo"]]);r.focusHi=Num(f[c["focusHi"]]);
    Require(r.fitHi>r.fitLo&&r.focusHi>r.focusLo,"Invalid focus/fit extent");
    Require(out.insert(std::make_pair(Key(r.label,r.index),r)).second,"Duplicate quality model/index");
  }
  Require(in.eof()&&!out.empty(),"Unreadable quality TSV");return out;
}
void Annotate(const Row& r,bool focus){
  gPad->SetTopMargin(.34);TIter next(gPad->GetListOfPrimitives());TObject* o=0;
  while((o=next())){
    TLegend* l=dynamic_cast<TLegend*>(o);
    if(l){l->SetX1NDC(.69);l->SetX2NDC(.96);l->SetY1NDC(.40);l->SetY2NDC(.65);l->SetTextSize(.025);}
  }
  TLatex t;t.SetNDC();t.SetTextFont(42);t.SetTextSize(.025);
  t.DrawLatex(.13,.845,Form("target %s; Pearson #chi^{2}/conditional ndf = %s/%s = %s",r.target.c_str(),r.chi2.c_str(),r.ndf.c_str(),r.ratio.c_str()));
  t.SetTextColor(!r.numerical?kRed+1:(!r.supportKnown||!r.supported?kOrange+7:kGreen+3));
  t.DrawLatex(.13,.805,Form("numerical: %s; non-#chi^{2} support: %s",r.numerical?"pass":"FAIL",!r.supportKnown?"not evaluated":(r.supported?"pass":"FAIL")));
  t.SetTextColor(kBlack);t.SetTextSize(.022);std::string reason=r.reason.empty()?"none":r.reason;
  for(unsigned a=0,line=0;a<reason.size()&&line<2;a+=104,++line){
    std::string part=reason.substr(a,104);if(line==1&&a+104<reason.size())part+="...";
    t.DrawLatex(.13,.765-.034*line,(std::string(line?"        ":"reason: ")+part).c_str());
  }
  t.DrawLatex(.13,.682,Form("%s; saved fit [%.2f, %.2f]; no refit/renormalization",focus?"Nuclear-neighborhood view":"Wide observed distribution",r.fitLo,r.fitHi));
  gPad->Modified();gPad->Update();
}
}
#endif
void PlotQualityDiagnostics(const char* modelsPath,const char* slicesPath,const char* settingsPath,const char* qualityPath,const char* pdf){
#ifndef __CINT__
  using namespace CutoffSpectra;
  try{
    Require(modelsPath&&slicesPath&&settingsPath&&qualityPath&&pdf,"Explicit view inputs required");
    Require(gSystem->AccessPathName(pdf),"Existing PDF will not be overwritten");
    Settings cfg=ReadSettings(settingsPath);
    double low=Number(cfg,"dedxMin"),high=Number(cfg,"dedxMax"),floor=Number(cfg,"logMinimum"),group=Number(cfg,"displayDedxRebinFactor");
    Require(high>low&&floor>0&&group>=1&&group==std::floor(group),"Invalid display setting");
    std::vector<Model> models=ReadModels(modelsPath);std::vector<Selection> selected=ReadSelections(slicesPath);
    QualityDiagnostics::Rows quality=QualityDiagnostics::Read(qualityPath);
    gStyle->SetOptStat(0);gStyle->SetPaperSize(26,20);
    TCanvas canvas("quality_diagnostics","Saved nuclear local fit diagnostics",1400,1100);
    canvas.Print((std::string(pdf)+"[").c_str());int serial=0;
    for(unsigned i=0;i<selected.size();++i){
      const Selection& s=selected[i];double a=high,b=low;
      for(unsigned m=0;m<models.size();++m){
        VerifyRaw(*Load(*models[0].file,"raw_native",s.index),*Load(*models[m].file,"raw_native",s.index));
        QualityDiagnostics::Rows::const_iterator it=quality.find(QualityDiagnostics::Key(models[m].label,s.index));
        Require(it!=quality.end(),"Missing model/index quality row");
        a=std::min(a,it->second.focusLo);b=std::max(b,it->second.focusHi);
      }
      a=std::max(low,a);b=std::min(high,b);Require(b>a,"Empty nuclear neighborhood");
      for(int focus=0;focus<2;++focus)for(unsigned pair=0;pair<models.size();pair+=2){
        canvas.Clear();canvas.Divide(2,2);
        for(unsigned m=0;m<2&&pair+m<models.size();++m)for(int log=0;log<2;++log){
          const Model& model=models[pair+m];canvas.cd(2*m+log+1);
          DrawPanel(*model.file,model.label,s,int(group),focus?a:low,focus?b:high,log!=0,floor,serial++);
          QualityDiagnostics::Annotate(quality.find(QualityDiagnostics::Key(model.label,s.index))->second,focus!=0);
        }
        canvas.Print(pdf);
      }
    }
    canvas.Print((std::string(pdf)+"]").c_str());
    for(unsigned m=0;m<models.size();++m){models[m].file->Close();delete models[m].file;}
    std::cout<<"QUALITY_DIAGNOSTICS_COMPLETE pages="<<2*((models.size()+1)/2)*selected.size()<<std::endl;
  }catch(const std::exception& e){std::cerr<<"PlotQualityDiagnostics: "<<e.what()<<std::endl;gSystem->Exit(2);}
#endif
}
