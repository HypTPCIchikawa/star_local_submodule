// DISPLAY ONLY. Reuse an immutable renderer snapshot; no fitting or PID writes.
#include "FrozenRenderIterativeHelium.C"
#include "TGaxis.h"

#ifndef __CINT__
namespace WideHeliumViews {
using namespace IterativeHeliumViews;
void Range(Slice& s,double lo,double hi,const Fields& c) {
  Require(lo>=s.native->GetXaxis()->GetXmin()-Number(c,"axisRelativeTolerance")&&hi<=s.native->GetXaxis()->GetXmax()+Number(c,"axisRelativeTolerance")&&hi>lo,"Physical display outside native axis");
  TAxis axis(*s.raw->GetXaxis());axis.SetRangeUser(lo,hi);
  s.firstDisplay=std::max(1,axis.GetFirst());s.lastDisplay=std::min(s.raw->GetNbinsX(),axis.GetLast());
  Require(s.lastDisplay>=s.firstDisplay,"Empty physical display");
}
void NumericAudit(const std::vector<Slice>& slices,const Fields& c,const char* out) {
  std::ofstream audit((std::string(out)+"/wide_numeric_audit.tsv").c_str());
  std::ofstream mask((std::string(out)+"/visible_native_mask.tsv").c_str());
  Require(audit.good()&&mask.good(),"Cannot write wide numerical audit");audit<<std::setprecision(17);
  audit<<"index\trigidity\tnativeNormalCount\trawNormalCount\tnativeNormalVariance\trawNormalVariance\tunderflow\toverflow\tunderflowVariance\toverflowVariance\tfirstVisibleNative\tlastVisibleNative\tfirstVisibleFitBin\tlastVisibleFitBin\tclosureRelativeError\n";
  mask<<"index\tnativeBin\tphysicalLo\tphysicalHi\trebinnedBin\tvisibleFull\tvisibleSupplemental\tusedInFit\n";mask<<std::setprecision(17);
  const int factor=Integer(c,"dedxRebinFactor");
  for(unsigned i=0;i<slices.size();++i) {
    const Slice& s=slices[i];double nv=0,rv=0;
    for(int b=1;b<=s.native->GetNbinsX();++b)nv+=std::pow(s.native->GetBinError(b),2);
    for(int b=1;b<=s.raw->GetNbinsX();++b)rv+=std::pow(s.raw->GetBinError(b),2);
    Require(s.firstDisplay==1&&s.lastDisplay==s.raw->GetNbinsX(),"Full frame does not include all normal native bins");
    audit<<s.index<<'\t'<<s.r<<'\t'<<s.native->Integral(1,s.native->GetNbinsX())<<'\t'<<s.raw->Integral(1,s.raw->GetNbinsX())<<'\t'<<nv<<'\t'<<rv<<'\t'<<s.raw->GetBinContent(0)<<'\t'<<s.raw->GetBinContent(s.raw->GetNbinsX()+1)<<'\t'<<std::pow(s.raw->GetBinError(0),2)<<'\t'<<std::pow(s.raw->GetBinError(s.raw->GetNbinsX()+1),2)<<"\t1\t"<<s.native->GetNbinsX()<<'\t'<<s.firstDisplay<<'\t'<<s.lastDisplay<<'\t'<<s.maximumClosureRelativeError<<'\n';
    for(int b=1;b<=s.native->GetNbinsX();++b) {
      const int fitbin=(b-1)/factor+1;
      const double centre=s.raw->GetBinCenter(fitbin);
      const bool supplemental=s.diagnostic&&centre>=Number(c,"supplementalPhysicalMin")&&centre<Number(c,"supplementalPhysicalMax");
      mask<<s.index<<'\t'<<b<<'\t'<<s.native->GetXaxis()->GetBinLowEdge(b)<<'\t'<<s.native->GetXaxis()->GetBinUpEdge(b)<<'\t'<<fitbin<<"\t1\t"<<int(supplemental)<<'\t'<<int(fitbin>=s.firstFit&&fitbin<=s.lastFit)<<'\n';
    }
  }
  audit.close();mask.close();Require(audit.good()&&mask.good(),"Cannot finish wide numerical audit");
}
void Diagnostics(const std::vector<Slice>& slices,const Fields& c,const TString& book,const char* output,const char* label,bool supplemental) {
  TCanvas page(supplemental?"wide_nuclear_components":"wide_full_components","Saved same-bin views",Integer(c,"diagnosticCanvasWidth"),Integer(c,"diagnosticCanvasHeight"));page.Print(book+"[");
  for(unsigned i=0;i<slices.size();++i)if(slices[i].diagnostic) {
    Slice s=slices[i];if(supplemental)Range(s,Number(c,"supplementalPhysicalMin"),Number(c,"supplementalPhysicalMax"),c);
    page.Clear();page.Divide(2,2);
    for(int p=0;p<4;++p){page.cd(p+1);Draw(s,c,p>=2,p%2==0,label);}
    page.Print(book);page.Print(TString(output)+Form(supplemental?"/wide_nuclear_R_%03d.png":"/wide_full_R_%03d.png",int(s.r*1000+.5)));
  }page.Print(book+"]");
}
}
#endif

void WideIterativeHelium(const char* inputRoot,const char* manifest,const char* settings,const char* output,const char* phaseLabel) {
#ifndef __CINT__
  using namespace IterativeHeliumViews;using namespace WideHeliumViews;
  try {
    Fields c=ReadSettings(settings);ValidateSettings(c);
    Require(Integer(c,"axisMaxDigits")>0,"axisMaxDigits must be explicit positive integer");
    Require(Number(c,"fullPhysicalMin")==0&&Number(c,"fullPhysicalMax")==120,"Authorized full frame is native 0-120 keV/cm");
    Require(Number(c,"supplementalPhysicalMin")==45&&Number(c,"supplementalPhysicalMax")==120,"Supplemental nuclear frame must be 45-120");
    const std::vector<Fields> rows=ReadManifest(manifest);Require(rows.size()==unsigned(Integer(c,"expectedSlices")),"Unexpected slice count");
    TFile input(inputRoot,"READ");Require(!input.IsZombie(),"Cannot read saved phase ROOT");std::vector<Slice> slices;
    for(unsigned i=0;i<rows.size();++i) {
      Slice s=LoadSlice(&input,rows[i],c);Require(s.index==int(i),"Wide slices must start at zero and be contiguous");
      Require(Near(s.native->GetXaxis()->GetXmin(),Number(c,"fullPhysicalMin"),Number(c,"axisRelativeTolerance"))&&Near(s.native->GetXaxis()->GetXmax(),Number(c,"fullPhysicalMax"),Number(c,"axisRelativeTolerance")),"Native axis differs from explicit full frame");
      Range(s,Number(c,"fullPhysicalMin"),Number(c,"fullPhysicalMax"),c);slices.push_back(s);
    }
    gStyle->SetOptStat(0);gStyle->SetOptFit(0);TGaxis::SetMaxDigits(Integer(c,"axisMaxDigits"));gStyle->SetPaperSize(Number(c,"paperWidthCm"),Number(c,"paperHeightCm"));
    TCanvas page("wide_iterative_pages","Native full physical range",Integer(c,"linearCanvasWidth"),Integer(c,"linearCanvasHeight"));
    for(int physical=0;physical<2;++physical)for(int logarithmic=0;logarithmic<2;++logarithmic) {
      const TString book=TString(output)+Form("/wide_%s_%s.pdf",physical?"raw_dedx":"nsigma",logarithmic?"log":"linear");Require(gSystem->AccessPathName(book),"Refusing overwrite of wide PDF");page.Print(book+"[");
      for(unsigned start=0;start<slices.size();start+=4){page.Clear();page.Divide(2,2);for(unsigned j=0;j<4&&start+j<slices.size();++j){page.cd(j+1);Draw(slices[start+j],c,physical,logarithmic,phaseLabel);}page.Print(book);if(start==0)page.Print(TString(output)+Form("/wide_%s_%s_first.png",physical?"raw_dedx":"nsigma",logarithmic?"log":"linear"));}page.Print(book+"]");
    }
    const std::string fullLabel=std::string("FULL 0-120 / ")+phaseLabel,nuclearLabel=std::string("SUPPLEMENT 45-120 / ")+phaseLabel;
    Diagnostics(slices,c,TString(output)+"/wide_low_rigidity_components.pdf",output,fullLabel.c_str(),false);
    Diagnostics(slices,c,TString(output)+"/supplemental_nuclear45_120_components.pdf",output,nuclearLabel.c_str(),true);
    NumericAudit(slices,c,output);WriteAudit(slices,c,output,inputRoot,phaseLabel);
    std::cout<<"WIDE_ITERATIVE_HELIUM_COMPLETE slices="<<slices.size()<<"; no refit; data only outside original fit mask"<<std::endl;
  }catch(const std::exception& e){std::cerr<<"WIDE_ITERATIVE_HELIUM_ERROR "<<e.what()<<std::endl;gSystem->Exit(1);}
#endif
}
