#include <TFile.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TStyle.h>
#include <TColor.h>
#include <TEnv.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TNamed.h>
#include <TSystem.h>
#include <TString.h>
#include <fstream>
#include <sstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
std::string ReRequired(TEnv& e,const char* key){std::string s=e.GetValue(key,"");if(s.empty())throw std::runtime_error(std::string("Missing render setting ")+key);return s;}
void ReText(double x,double y,const char* text,double pixels){TLatex latex;latex.SetNDC();latex.SetTextFont(43);latex.SetTextSize(pixels);latex.DrawLatex(x,y,text);}
void ReAxes(TH2D* h,double xmin,double xmax,double ymin,double ymax,double titlepx,double labelpx){
 h->SetTitle("");h->SetStats(kFALSE);h->GetXaxis()->SetRangeUser(xmin,xmax);h->GetYaxis()->SetRangeUser(ymin,ymax);
 h->GetXaxis()->SetTitle("#eta_{lab}");h->GetYaxis()->SetTitle("#beta = |p|/E");h->GetZaxis()->SetTitle("Candidates / bin");
 for(int a=0;a<3;++a){TAxis* axis=a==0?h->GetXaxis():(a==1?h->GetYaxis():h->GetZaxis());axis->SetTitleFont(43);axis->SetLabelFont(43);axis->SetTitleSize(a==2?18:titlepx);axis->SetLabelSize(a==2?15:labelpx);axis->SetNdivisions(505);axis->CenterTitle(kTRUE);}
 h->GetXaxis()->SetNdivisions(510);h->GetXaxis()->SetTitleOffset(1.2);h->GetYaxis()->SetTitleOffset(1.35);h->GetZaxis()->SetTitleOffset(1.5);
}
void RenderBetaEta(const char* configPath){
 try{
  TEnv e;if(e.ReadFile(configPath,kEnvLocal)!=0)throw std::runtime_error("Cannot read explicit render config");
  const std::string cache=ReRequired(e,"Cache"),out=ReRequired(e,"Output"),source=ReRequired(e,"SourceCaption");
  const double xmin=e.GetValue("EtaMin",-999.),xmax=e.GetValue("EtaMax",-999.),ymin=e.GetValue("BetaMin",-999.),ymax=e.GetValue("BetaMax",-999.),zmax=e.GetValue("CommonMaximum",-999.);
  const int titlepx=e.GetValue("AxisTitlePixels",-1),labelpx=e.GetValue("AxisLabelPixels",-1),panelpx=e.GetValue("PanelTitlePixels",-1),fill=e.GetValue("OverlayFillStyle",-1);
  if(xmin>=xmax||ymin>=ymax||zmax<=0||titlepx<=0||labelpx<=0||panelpx<=0||fill<0)throw std::runtime_error("Incomplete plot geometry settings");
  TH1::AddDirectory(kFALSE);TFile f(cache.c_str(),"READ");if(f.IsZombie()||f.TestBit(TFile::kRecovered))throw std::runtime_error("Bad reviewed histogram cache");
  std::vector<TH2D*> h;std::vector<std::string> keys,titles,selection;std::vector<int> colors;
  for(int i=0;i<3;++i){keys.push_back(ReRequired(e,TString::Format("Sample.%d.Key",i)));titles.push_back(ReRequired(e,TString::Format("Sample.%d.Title",i)));selection.push_back(ReRequired(e,TString::Format("Sample.%d.Selection",i)));colors.push_back(TColor::GetColor(ReRequired(e,TString::Format("Sample.%d.Color",i)).c_str()));TH2D* original=dynamic_cast<TH2D*>(f.Get(keys.back().c_str()));if(!original)throw std::runtime_error("Missing reviewed TH2D");h.push_back((TH2D*)original->Clone(TString::Format("overlay_%d",i)));h.back()->SetDirectory(0);delete original;}
  f.Close();gStyle->SetOptStat(0);gStyle->SetCanvasColor(0);gStyle->SetPadColor(0);gStyle->SetFrameFillColor(0);
  Double_t stops[5]={0,.25,.5,.75,1},red[5]={.267,.230,.128,.369,.993},green[5]={.005,.322,.567,.789,.906},blue[5]={.329,.546,.551,.383,.144};
  int first=TColor::CreateGradientColorTable(5,stops,red,green,blue,255),palette[255];for(int i=0;i<255;++i)palette[i]=first+i;gStyle->SetPalette(255,palette);gStyle->SetNumberContours(255);
  {TCanvas c("beta_eta_overview_revised","beta eta overview",1800,700);c.Divide(3,1,0,0);
   for(int i=0;i<3;++i){c.cd(i+1);gPad->SetLeftMargin(.13);gPad->SetRightMargin(.17);gPad->SetBottomMargin(.18);gPad->SetTopMargin(.21);gPad->SetTicks(1,1);gPad->SetLogz();
    ReAxes(h[i],xmin,xmax,ymin,ymax,titlepx,labelpx);h[i]->SetMinimum(.5);h[i]->SetMaximum(zmax);h[i]->Draw("COLZ");
    ReText(.13,.953,titles[i].c_str(),panelpx);ReText(.13,.905,selection[i].c_str(),15);
    const double flow=h[i]->GetEntries()-h[i]->Integral(1,h[i]->GetNbinsX(),1,h[i]->GetNbinsY());ReText(.13,.850,(flow>0?TString::Format("N = %.0f (flow: %.0f)",h[i]->GetEntries(),flow):TString::Format("N = %.0f",h[i]->GetEntries())).Data(),16);
    ReText(.13,.044,source.c_str(),13);ReText(.13,.018,"Kinematic #beta = |p|/E; lab frame",13);gPad->Modified();gPad->Update();
   }
   c.Print((out+"/page1_revised.pdf").c_str());c.Print((out+"/overview.png").c_str());
  }
  {TCanvas c("beta_eta_box_overlay","beta eta BOX overlay",1100,850);c.cd();gPad->SetLeftMargin(.13);gPad->SetRightMargin(.055);gPad->SetBottomMargin(.16);gPad->SetTopMargin(.20);gPad->SetTicks(1,1);gPad->SetLogz(0);
   TH2D* frame=(TH2D*)h[0]->Clone("box_overlay_frame");frame->Reset("ICES");frame->SetDirectory(0);ReAxes(frame,xmin,xmax,ymin,ymax,24,20);frame->Draw("AXIS");
   for(int i=0;i<3;++i){ReAxes(h[i],xmin,xmax,ymin,ymax,24,20);h[i]->SetMinimum(0);h[i]->SetMaximum(zmax);h[i]->SetLineColor(colors[i]);h[i]->SetFillColor(colors[i]);h[i]->SetLineWidth(1);h[i]->SetFillStyle(fill);h[i]->Draw("BOX SAME");}
   frame->Draw("AXIS SAME");TLegend legend(.15,.615,.66,.79);legend.SetTextFont(43);legend.SetTextSize(18);legend.SetFillStyle(0);legend.SetBorderSize(0);
   legend.AddEntry(h[0],"#Lambda signal","l");legend.AddEntry(h[1],"^{4}He: 0.61 < R < 2, TOF not required","l");legend.AddEntry(h[2],"^{4}He: 2 #leq R < 5, TOF + m^{2}","l");legend.Draw();
   ReText(.13,.955,"#beta vs #eta: three samples, BOX overlay",26);ReText(.13,.909,"Same #eta and #beta ranges as page 1",17);ReText(.13,.869,"Common linear raw-count scale; no per-sample normalization",16);
   ReText(.13,.046,source.c_str(),14);ReText(.13,.021,"Kinematic #beta = |p|/E; colored box outlines",14);c.Modified();c.Update();
   c.Print((out+"/beta_eta_box_overlay.pdf").c_str());c.Print((out+"/beta_eta_box_overlay.png").c_str());
   {TFile output((out+"/beta_eta_box_overlay.root").c_str(),"RECREATE");if(output.IsZombie())throw std::runtime_error("Cannot create overlay ROOT");for(int i=0;i<3;++i)if(h[i]->Write(keys[i].c_str())<=0)throw std::runtime_error("Cannot save overlay TH2D");if(c.Write("beta_eta_box_overlay_canvas")<=0)throw std::runtime_error("Cannot save overlay canvas");std::ifstream config(configPath);std::ostringstream contents;contents<<config.rdbuf();TNamed metadata("BetaEtaOverlayConfiguration",contents.str().c_str());metadata.Write();TNamed mode("BetaEtaOverlayDrawOption","BOX SAME; raw counts; shared linear maximum; fill style 0");mode.Write();output.Close();}
   delete frame;
  }
  for(unsigned i=0;i<h.size();++i)delete h[i];
  std::cout<<"BETA_ETA_LAYOUT_OVERLAY_DONE"<<std::endl;gSystem->Exit(0);
 }catch(const std::exception& error){std::cerr<<"ERROR: "<<error.what()<<std::endl;gSystem->Exit(1);}
}
