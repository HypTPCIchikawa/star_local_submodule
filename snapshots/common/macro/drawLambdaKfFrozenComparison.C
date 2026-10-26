// ROOT5 ACLiC-only visualization of independently validated frozen KF cuts.
// Reads comparator CSVs; no fits, new selections, rebinning, or mass constraint.
// drawLambdaKfFrozenComparison(comparisonDirectory, NEW_outputPrefix)
#include "drawLambdaPurityStudy.C"

Int_t drawLambdaKfFrozenComparison(const char*, const char*);
#ifndef __CINT__
namespace LambdaFrozenComparison {
using LambdaPurityQA::Row;
using LambdaPurityQA::Require;
using LambdaPurityQA::Get;
using LambdaPurityQA::N;
using LambdaPurityQA::ReadCsv;
using LambdaPurityQA::Edge;

static Row Manifest(const std::string& path) {
  std::ifstream input(path.c_str()); Require(bool(input), "cannot read plot manifest");
  Row values; std::string line;
  while (std::getline(input,line)) {
    if (line.empty() || line[0]=='#') continue;
    const size_t tab=line.find('\t'); Require(tab!=std::string::npos,"manifest requires TAB");
    std::string value=line.substr(tab+1);
    if(!value.empty() && value[value.size()-1]=='\r') value.erase(value.size()-1);
    Require(values.insert(std::make_pair(line.substr(0,tab),value)).second,"duplicate manifest key");
  }
  Require(input.eof(),"manifest read error");
  Require(Get(values,"analysisRole")=="independent_frozen_validation","not frozen independent validation");
  return values;
}

static void Style(TH1D* h, int index) {
  h->SetLineColor(index ? kAzure+2 : kGray+2);
  h->SetMarkerColor(index ? kAzure+2 : kGray+2);
  h->SetMarkerStyle(index ? 20 : 24);
  h->SetMarkerSize(index ? .8 : .8);
  h->SetLineWidth(index ? 2 : 2);
  h->SetFillStyle(0);
}

static void Axis(TH1D* h, const char* ylabel) {
  h->SetTitle(""); h->GetXaxis()->SetTitle("M(p#pi^{-}) [GeV/c^{2}]");
  h->GetYaxis()->SetTitle(ylabel);
  h->GetXaxis()->SetTitleSize(.047); h->GetYaxis()->SetTitleSize(.047);
  h->GetXaxis()->SetLabelSize(.039); h->GetYaxis()->SetLabelSize(.039);
  h->GetXaxis()->SetTitleOffset(1.10); h->GetYaxis()->SetTitleOffset(1.50);
  h->GetXaxis()->SetNdivisions(505); h->GetYaxis()->SetNdivisions(505);
}
} // namespace LambdaFrozenComparison

Int_t drawLambdaKfFrozenComparison(const char* comparisonDirectory, const char* outputPrefix) {
  using namespace LambdaFrozenComparison;
  try {
    Require(comparisonDirectory && *comparisonDirectory && outputPrefix && *outputPrefix,"missing arguments");
    const std::string dir(comparisonDirectory), stem(outputPrefix);
    const char* ext[]={".pdf",".png",".root","_audit.txt"};
    for(int i=0;i<4;++i) Require(gSystem->AccessPathName((stem+ext[i]).c_str()),"refusing overwrite "+stem+ext[i]);
    const Row config=Manifest(dir+"/plot_manifest.tsv");
    const double lo=N(config,"massMin"), hi=N(config,"massMax"), width=N(config,"massBinWidth");
    const double displayLo=N(config,"displayMin"), displayHi=N(config,"displayMax");
    const double windowLo=N(config,"windowLow"), windowHi=N(config,"windowHigh"), events=N(config,"inputEvents");
    Require(lo<hi && width>0 && displayLo>=lo && displayHi<=hi && displayLo<displayHi,"invalid histogram/display ranges");
    Require(events>0 && events==std::floor(events),"invalid input-event count");
    const int nbins=int(std::floor((hi-lo)/width+.5));
    Require(nbins>0 && std::fabs(lo+nbins*width-hi)<1.e-9,"invalid uniform binning");
    const char* ids[]={"baseline","best50"};
    TH1D* hist[2];
    std::vector<bool> seen[2];
    double norms[2]={-1.,-1.};
    for(int j=0;j<2;++j) {
      hist[j]=new TH1D(Form("FrozenNativeCounts_%s",ids[j]),"",nbins,lo,hi);
      hist[j]->SetDirectory(0); hist[j]->Sumw2(); seen[j].assign(nbins,false);
    }
    std::vector<Row> rows=ReadCsv(dir+"/histograms.csv");
    for(size_t i=0;i<rows.size();++i) {
      const Row& row=rows[i]; Require(Get(row,"split")=="full","unexpected non-full histogram");
      const std::string id=Get(row,"choice_id");
      const int j=id==ids[0] ? 0 : (id==ids[1] ? 1 : -1);
      Require(j>=0,"unexpected histogram selection");
      const int bin=int(N(row,"bin"));
      Require(bin==N(row,"bin") && bin>=1 && bin<=nbins && !seen[j][bin-1],"invalid/duplicate histogram bin");
      Require(std::fabs(N(row,"low")-hist[j]->GetBinLowEdge(bin))<1.e-9 &&
              std::fabs(N(row,"high")-hist[j]->GetBinLowEdge(bin+1))<1.e-9,"histogram edges differ");
      const double count=N(row,"count"), error=N(row,"error"), norm=N(row,"normalization_count");
      Require(count>=0 && count==std::floor(count) && error>=0 &&
              std::fabs(error*error-count)<1.e-8*std::max(1.,count),"non-Poisson or weighted input histogram");
      if(norms[j]<0) norms[j]=norm; else Require(norms[j]==norm,"inconsistent normalization metadata");
      hist[j]->SetBinContent(bin,count); hist[j]->SetBinError(bin,error); seen[j][bin-1]=true;
    }
    for(int j=0;j<2;++j) {
      for(int b=0;b<nbins;++b) Require(seen[j][b],"missing native histogram bin");
      Require(hist[j]->Integral(1,nbins)==norms[j] && norms[j]>0,"invalid native normalization count");
      Edge(hist[j],displayLo); Edge(hist[j],displayHi); Edge(hist[j],windowLo); Edge(hist[j],windowHi);
    }
    Row metrics[2]; bool got[2]={false,false};
    rows=ReadCsv(dir+"/chosen_metrics.csv");
    for(size_t i=0;i<rows.size();++i) {
      const Row& row=rows[i]; Require(Get(row,"split")=="full","unexpected non-full metric");
      if(Get(row,"sideband")!="default") continue;
      const std::string id=Get(row,"choice_id");
      const int j=id==ids[0] ? 0 : (id==ids[1] ? 1 : -1);
      Require(j>=0 && !got[j],"unknown/duplicate default-sideband metric");
      metrics[j]=row; got[j]=true;
      Require(hist[j]->Integral(Edge(hist[j],windowLo),Edge(hist[j],windowHi)-1)==N(row,"Nwin"),"fixed-window metric disagrees with histogram");
      Require(std::fabs(N(row,"S")+N(row,"B")-N(row,"Nwin"))<1.e-8,"inconsistent signal/background counts");
      Require(N(row,"Nwin")>0 && std::fabs(N(row,"purity")-N(row,"S")/N(row,"Nwin"))<1.e-10,"inconsistent purity");
    }
    Require(got[0] && got[1] && N(metrics[0],"S")>0,"missing/nonpositive baseline metrics");
    Require(std::fabs(N(metrics[1],"signalRetention")-N(metrics[1],"S")/N(metrics[0],"S"))<1.e-10,"incorrect common-baseline retention");
    gSystem->mkdir(gSystem->DirName(stem.c_str()),kTRUE);
    gROOT->SetBatch(kTRUE); gStyle->SetOptStat(0); gStyle->SetTextFont(42); gStyle->SetEndErrorSize(2);
    TFile output((stem+".root").c_str(),"CREATE"); Require(!output.IsZombie(),"cannot create output ROOT file");
    for(int j=0;j<2;++j) hist[j]->Write();
    TCanvas* canvas=new TCanvas("LambdaKfFrozenComparison","Frozen KF Lambda comparison",1800,1040);
    TH1D* legendHist[2]={0,0};
    for(int panel=0;panel<2;++panel) {
      canvas->cd();
      TPad* pad=new TPad(Form("FrozenComparisonPad%d",panel),"",panel*.5,.055,(panel+1)*.5,.785);
      pad->SetLeftMargin(.155); pad->SetRightMargin(.045); pad->SetBottomMargin(.135);
      pad->SetTopMargin(.075); pad->SetTicks(1,1); pad->SetGridy(); pad->Draw(); pad->cd();
      TH1D* drawn[2]; double maximum=0;
      for(int j=0;j<2;++j) {
        drawn[j]=(TH1D*)hist[j]->Clone(Form("FrozenDisplay_P%d_%s",panel,ids[j]));
        drawn[j]->SetDirectory(0); if(panel) drawn[j]->Scale(1./norms[j]);
        Style(drawn[j],j);
        if(panel) Require(std::fabs(drawn[j]->Integral(1,nbins)-1.)<1.e-12,"full-range normalized integral differs from one");
        for(int b=1;b<=nbins;++b) {
          const double scale=panel ? norms[j] : 1.;
          Require(std::fabs(drawn[j]->GetBinContent(b)-hist[j]->GetBinContent(b)/scale)<1.e-10,"drawn content changed");
          Require(std::fabs(drawn[j]->GetBinError(b)-hist[j]->GetBinError(b)/scale)<1.e-10,"drawn error changed");
        }
        drawn[j]->GetXaxis()->SetRange(Edge(drawn[j],displayLo),Edge(drawn[j],displayHi)-1);
        for(int b=Edge(drawn[j],displayLo);b<Edge(drawn[j],displayHi);++b)
          maximum=std::max(maximum,drawn[j]->GetBinContent(b)+drawn[j]->GetBinError(b));
        if(!panel) legendHist[j]=drawn[j];
        output.cd(); drawn[j]->Write(); pad->cd();
      }
      TH1D* frame=(TH1D*)drawn[0]->Clone(Form("FrozenFrame%d",panel));
      frame->SetDirectory(0); frame->Reset("ICES"); frame->SetMinimum(0); frame->SetMaximum(std::max(1.e-8,maximum*1.18));
      Axis(frame,Form(panel ? "Fraction per %.3g MeV/c^{2} bin" : "Candidates per %.3g MeV/c^{2} bin",width*1000));
      frame->Draw("AXIS");
      for(int j=0;j<2;++j) { drawn[j]->Draw("HIST SAME"); drawn[j]->Draw("E1 SAME"); }
      TLatex label; label.SetNDC(); label.SetTextFont(42); label.SetTextSize(.042);
      label.DrawLatex(.155,.955,panel ? "Unit-area normalized shape" : "Raw candidate counts");
      pad->RedrawAxis();
    }
    canvas->cd();
    TLatex text; text.SetNDC(); text.SetTextFont(42); text.SetTextSize(.029);
    text.DrawLatex(.045,.963,Get(config,"title").c_str());
    TLegend* legend=new TLegend(.045,.825,.96,.930);
    legend->SetBorderSize(0); legend->SetFillStyle(0); legend->SetTextFont(42); legend->SetTextSize(.025); legend->SetMargin(.06);
    legend->AddEntry(legendHist[0],Form("Default: N_{win} = %.0f, purity = %.2f%%",N(metrics[0],"Nwin"),100*N(metrics[0],"purity")),"lep");
    legend->AddEntry(legendHist[1],Form("Optimized (frozen 10k): N_{win} = %.0f, purity = %.2f%%, S/S_{default} = %.1f%%",N(metrics[1],"Nwin"),100*N(metrics[1],"purity"),100*N(metrics[1],"signalRetention")),"lep");
    legend->Draw(); text.SetTextSize(.021);
    text.DrawLatex(.045,.795,"Same input events for both selections; independent of the cut-discovery sample. No retuning.");
    text.SetTextSize(.019);
    text.DrawLatex(.045,.022,Form("Fixed mass window [%.3f, %.3f); unit-area normalization [%.2f, %.2f). Sideband purity estimates; no mass constraint, rebinning or smoothing.",windowLo,windowHi,lo,hi));
    LambdaPurityQA::Save(canvas,stem,output);
    std::ostringstream audit;
    audit<<"Independent frozen-cut validation visualization; actual reconstruction selections.\nsource="<<dir
         <<"\ninputEvents="<<events<<"\ncomplete-range-normalization="<<lo<<":"<<hi
         <<"\ndisplay="<<displayLo<<":"<<displayHi
         <<"\nAll native bins, count errors, normalization sums, and fixed-window metrics checked.\n"
         <<"Raw counts and Poisson error bars preserved; normalized errors divided by original total.\n"
         <<"Normalization-total uncertainty/covariance is not added to the shape error bars.\n"
         <<"Purity is a fixed-sideband estimate, not a truth label or proof of zero background.\n";
    output.cd(); TObjString(audit.str().c_str()).Write("FrozenComparisonAudit"); output.Flush();
    Require(!output.TestBit(TFile::kWriteError),"ROOT write error"); output.Close();
    std::ofstream log((stem+"_audit.txt").c_str()); log<<audit.str(); log.close(); Require(bool(log),"audit write error");
    std::cout<<"PASS frozen comparison plot: "<<stem<<"; no cuts changed."<<std::endl;
    return 0;
  } catch(const std::exception& error) {
    std::cerr<<"drawLambdaKfFrozenComparison ERROR: "<<error.what()<<std::endl; return 1;
  }
}
#endif // __CINT__
