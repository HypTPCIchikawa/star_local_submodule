// Sideband-model sensitivity only. No selection or misleading holdout overlay.
#include "checkLambdaPurityRefinement.C"
Int_t checkLambdaKfPurityBackground(const char* directory, const char* prefix);
#ifndef __CINT__
Int_t checkLambdaKfPurityBackground(const char* directory, const char* prefix) {
  using namespace LambdaPurityRefinement;
  try {
    Require(directory && *directory && prefix && *prefix,"missing arguments");
    std::string dir(directory),stem(prefix);
    Config c=Manifest(dir+"/plot_manifest.tsv");
    std::vector<Row> rows=ReadCsv(dir+"/chosen_selections.csv");
    std::vector<Choice> choices;
    std::map<std::string,int> ids;
    const int bins=int(std::floor((c.hi-c.lo)/c.bw+.5));
    Require(bins>0 && std::fabs(c.lo+bins*c.bw-c.hi)<1.e-9,"invalid bins");
    std::vector<TH1D*> hist;
    for(size_t i=0;i<rows.size();++i) {
      Choice choice;choice.id=Get(rows[i],"choice_id");choice.label=Get(rows[i],"label");choice.role=Get(rows[i],"role");
      Require(ids.insert(std::make_pair(choice.id,int(i))).second,"duplicate choice");
      choices.push_back(choice);
      TH1D* h=new TH1D(Form("KfBackgroundInput_%u",(unsigned)i),"",bins,c.lo,c.hi);
      h->SetDirectory(0);h->Sumw2();hist.push_back(h);
    }
    Require(!hist.empty(),"no histograms");
    std::vector<std::vector<bool> > seen(hist.size(),std::vector<bool>(bins,false));
    std::vector<double> norms(hist.size(),-1.);
    rows=ReadCsv(dir+"/histograms.csv");
    for(size_t i=0;i<rows.size();++i) {
      const Row&r=rows[i];if(Get(r,"split")!="full")continue;
      Require(ids.count(Get(r,"choice_id")),"unknown histogram");
      int j=ids[Get(r,"choice_id")],b=int(N(r,"bin"));
      Require(b==N(r,"bin")&&b>=1&&b<=bins&&!seen[j][b-1],"invalid/duplicate bin");seen[j][b-1]=true;
      Require(std::fabs(N(r,"low")-hist[j]->GetBinLowEdge(b))<1.e-9&&std::fabs(N(r,"high")-hist[j]->GetBinLowEdge(b+1))<1.e-9,"different edges");
      double n=N(r,"count"),e=N(r,"error");Require(n>=0&&std::fabs(n-std::floor(n+.5))<1.e-8&&e>=0&&std::fabs(e*e-n)<1.e-7*std::max(1.,n),"invalid count/error");
      hist[j]->SetBinContent(b,n);hist[j]->SetBinError(b,e);
      double norm=N(r,"normalization_count");if(norms[j]<0)norms[j]=norm;else Require(norms[j]==norm,"inconsistent norm");
    }
    for(size_t j=0;j<hist.size();++j) {
      for(int b=0;b<bins;++b)Require(seen[j][b],"missing bin");
      Require(hist[j]->Integral(1,bins)==norms[j],"normalization mismatch");
      const char*ext[]={".pdf",".png"};
      for(int k=0;k<2;++k)Require(gSystem->AccessPathName(Form("%s_background_%02u%s",prefix,(unsigned)j,ext[k])),"refusing diagnostic overwrite");
    }
    Require(gSystem->AccessPathName((stem+".root").c_str())&&gSystem->AccessPathName((stem+"_background.csv").c_str()),"refusing output overwrite");
    gROOT->SetBatch(kTRUE);gStyle->SetOptStat(0);
    TFile out((stem+".root").c_str(),"CREATE");Require(!out.IsZombie(),"cannot create output");
    Checks(choices,hist,c,stem,out);
    TObjString("Sideband-only nonnegative Bernstein degree0/1/2 sensitivity check. No cut ranking/reselection. Purity and asymptotic profile limits depend on background family; no peaking-background or feeddown exclusion.").Write("BackgroundScope");
    out.Flush();Require(!out.TestBit(TFile::kWriteError),"ROOT write failure");out.Close();
    std::cout<<"PASS background model sensitivity checks: "<<choices.size()<<" selections"<<std::endl;
    return 0;
  } catch(const std::exception&e) {std::cerr<<"FAIL background checks: "<<e.what()<<std::endl;return 1;}
}
#endif
