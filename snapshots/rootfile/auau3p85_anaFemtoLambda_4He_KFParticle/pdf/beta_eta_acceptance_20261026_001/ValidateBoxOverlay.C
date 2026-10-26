#include <TFile.h>
#include <TCanvas.h>
#include <TH2D.h>
#include <TEnv.h>
#include <TROOT.h>
#include <TColor.h>
#include <TSystem.h>
#include <TString.h>
#include <fstream>
#include <iostream>
#include <cmath>
void ValidateBoxOverlay(const char* configPath){
 TEnv e;if(e.ReadFile(configPath,kEnvLocal)!=0){gSystem->Exit(1);return;}
 TString out=e.GetValue("Output","");TH1::AddDirectory(kFALSE);TFile source(e.GetValue("Cache",""),"READ"),overlay(out+"/beta_eta_box_overlay.root","READ");if(source.IsZombie()||overlay.IsZombie()){gSystem->Exit(2);return;}
 TCanvas* canvas=dynamic_cast<TCanvas*>(overlay.Get("beta_eta_box_overlay_canvas"));if(!canvas){gSystem->Exit(3);return;}
 int mismatches=0;std::ofstream report((out+"/overlay_validation.tsv").Data());
 for(int i=0;i<3;++i){TString key=e.GetValue(TString::Format("Sample.%d.Key",i),"");TH2D* a=dynamic_cast<TH2D*>(source.Get(key)),*b=dynamic_cast<TH2D*>(overlay.Get(key));if(!a||!b||a->GetEntries()!=b->GetEntries()||a->GetNbinsX()!=b->GetNbinsX()||a->GetNbinsY()!=b->GetNbinsY()){gSystem->Exit(4);return;}
  for(int x=0;x<=a->GetNbinsX()+1;++x)for(int y=0;y<=a->GetNbinsY()+1;++y){int bin=a->GetBin(x,y);if(a->GetBinContent(x,y)!=b->GetBinContent(x,y)||static_cast<TH1*>(a)->GetBinError(bin)!=static_cast<TH1*>(b)->GetBinError(bin))++mismatches;}
  double xmin=b->GetXaxis()->GetBinLowEdge(b->GetXaxis()->GetFirst()),xmax=b->GetXaxis()->GetBinUpEdge(b->GetXaxis()->GetLast()),ymin=b->GetYaxis()->GetBinLowEdge(b->GetYaxis()->GetFirst()),ymax=b->GetYaxis()->GetBinUpEdge(b->GetYaxis()->GetLast());
  if(std::fabs(xmin-e.GetValue("EtaMin",-999.))>1e-12||std::fabs(xmax-e.GetValue("EtaMax",-999.))>1e-12||std::fabs(ymin-e.GetValue("BetaMin",-999.))>1e-12||std::fabs(ymax-e.GetValue("BetaMax",-999.))>1e-12||b->GetMaximum()!=e.GetValue("CommonMaximum",-999.)||b->GetFillStyle()!=e.GetValue("OverlayFillStyle",-999))++mismatches;
  TColor* color=gROOT->GetColor(b->GetLineColor());report<<"HIST\t"<<key<<"\t"<<a->GetEntries()<<"\t"<<xmin<<"\t"<<xmax<<"\t"<<ymin<<"\t"<<ymax<<"\t"<<(color?"color-loaded":"color-not-loaded")<<"\n";
 }
 canvas->Draw();canvas->Update();canvas->Print(out+"/box_overlay_canvas_readback.png");
 report<<"MISMATCHES\t"<<mismatches<<"\nDONE\t1\n";report.close();gSystem->Exit(mismatches?5:0);
}
