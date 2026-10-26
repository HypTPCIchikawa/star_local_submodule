void PlotNuclearParameters(const char*,const char*,const char*,const char*,double,double);
#ifndef __CINT__
#include "PlotNuclearCutoff.C"
#include "TGraphErrors.h"
#include <set>

void PlotNuclearParameters(const char* cachePath,const char* parameterPath,const char* currentPath,
                          const char* outputPdf,double rMin,double rMax) {
  using namespace NuclearCutoffPlot;
  try {
    if(!(rMax>rMin)||!gSystem->AccessPathName(outputPdf))throw std::runtime_error("Require explicit range and new output PDF");
    const std::vector<Row> rows=Read(parameterPath),current=Read(currentPath);
    const char* species[]={"d","t","he3","he4"};
    const char* labels[]={"d","t","^{3}He","^{4}He"};
    const char* models[]={"gaussian_linear","gaussian_exponential","lognormal_linear","lognormal_exponential"};
    const char* modelLabels[]={"Gaussian + linear","Gaussian + exponential","Log-normal + linear","Log-normal + exponential"};
    const int colors[]={kBlue+1,kAzure+7,kOrange+7,kMagenta+1};
    TFile cache(cachePath,"READ");TH2* all=dynamic_cast<TH2*>(cache.Get("ALL"));
    if(cache.IsZombie()||!all)throw std::runtime_error("Missing ALL cache");
    gStyle->SetOptStat(0);gStyle->SetPaperSize(26,20);
    TCanvas canvas("nuclear_parameters","ALL fitted nuclear responses",1300,1000);
    canvas.Print((std::string(outputPdf)+"[").c_str());
    for(int s=0;s<4;++s) {
      canvas.Clear();canvas.Divide(2,1);
      for(int parameter=0;parameter<2;++parameter) {
        const char* key=parameter?"sigma":"mean";
        const char* errKey=parameter?"sigmaError":"meanError";
        double maximum=0;
        for(unsigned i=0;i<current.size();++i)if(Text(current[i],"species")==species[s])maximum=std::max(maximum,Number(current[i],key));
        for(unsigned i=0;i<rows.size();++i)if(Text(rows[i],"species")==species[s]&&Number(rows[i],"targetFreeSupported")>.5)maximum=std::max(maximum,Number(rows[i],key));
        canvas.cd(parameter+1);gPad->SetTopMargin(.13);gPad->SetBottomMargin(.14);gPad->SetGridy();
        TH1* frame=gPad->DrawFrame(rMin,0,rMax,std::max(1.,maximum*1.12));
        frame->SetTitle(parameter?";R = |p/q| [GeV/c];Response width [keV/cm]":";R = |p/q| [GeV/c];Response center [keV/cm]");
        TGraph* reference=new TGraph;int n=0;
        for(unsigned i=0;i<current.size();++i)if(Text(current[i],"species")==species[s])reference->SetPoint(n++,Number(current[i],"r"),Number(current[i],key));
        reference->SetLineColor(kGray+2);reference->SetLineWidth(2);reference->Draw("L SAME");
        TLegend* legend=new TLegend(.40,.63,.88,.85);legend->SetBorderSize(0);legend->SetFillStyle(0);legend->SetTextSize(.023);
        legend->AddEntry(reference,"Current PID function (unchanged)","l");
        for(int model=0;model<4;++model) {
          TGraphErrors* graph=new TGraphErrors;int count=0;
          for(unsigned i=0;i<rows.size();++i) {
            if(Text(rows[i],"species")!=species[s]||Text(rows[i],"label").find(models[model])==std::string::npos||Number(rows[i],"targetFreeSupported")<.5)continue;
            const double r=Number(rows[i],"r"),value=Number(rows[i],key),error=Number(rows[i],errKey);
            if(r<rMin||r>rMax||!std::isfinite(value))continue;
            graph->SetPoint(count,r,value);graph->SetPointError(count,0,std::isfinite(error)&&error>0?error:0);++count;
          }
          graph->SetMarkerStyle(20+model);graph->SetMarkerSize(.42);graph->SetMarkerColor(colors[model]);graph->SetLineColor(colors[model]);
          if(graph->GetN())graph->Draw("P SAME");legend->AddEntry(graph,modelLabels[model],"p");
        }
        legend->Draw();Label(labels[s],"Own-fit support only; Gaussian: mean/std; log-normal: mode/curvature");
      }
      canvas.Print(outputPdf);
    }
    canvas.Clear();canvas.Divide(2,2);
    for(int s=0;s<4;++s) {
      canvas.cd(s+1);gPad->SetLogz();gPad->SetRightMargin(.16);gPad->SetTopMargin(.13);
      TH2* h=dynamic_cast<TH2*>(all->Clone((std::string("fitted_all_")+species[s]).c_str()));h->SetDirectory(0);
      h->GetXaxis()->SetRangeUser(rMin,rMax);h->SetTitle(";R = |p/q| [GeV/c];dE/dx [keV/cm]");h->Draw("COLZ");
      TGraph* guide=new TGraph;int count=0;
      for(unsigned i=0;i<current.size();++i)if(Text(current[i],"species")==species[s])guide->SetPoint(count++,Number(current[i],"r"),Number(current[i],"mean"));
      guide->SetLineColor(kRed+1);guide->SetLineStyle(2);guide->Draw("L SAME");
      TGraph* peaks=new TGraph;count=0;
      for(unsigned i=0;i<rows.size();++i)if(Text(rows[i],"species")==species[s]&&Text(rows[i],"label").find("gaussian_linear")!=std::string::npos&&Number(rows[i],"targetFreeSupported")>.5)
        peaks->SetPoint(count++,Number(rows[i],"r"),Number(rows[i],"mean"));
      peaks->SetMarkerStyle(20);peaks->SetMarkerSize(.42);peaks->SetMarkerColor(kBlack);if(peaks->GetN())peaks->Draw("P SAME");
      Label(labels[s],"ALL; black: Gaussian+linear fitted center; dashed: current PID");
    }
    canvas.Print(outputPdf);canvas.Print((std::string(outputPdf)+"]").c_str());
  } catch(const std::exception& error) {Error("PlotNuclearParameters","%s",error.what());gSystem->Exit(1);}
}
#endif
