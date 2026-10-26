void DrawIterativeSummary(const char* table,const char* settings,const char* output);
#ifndef __CINT__
#include "RenderIterativeHelium.C"
#include "TGraph.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TMath.h"
#include <limits>
namespace IterativeSummary {
struct Row {int phase,index,ndf,valid;std::string label;double r,mean3,sigma3,mean4,sigma4,deviance,chi2;};
void Panel(const std::vector<Row>& rows,const IterativeHeliumViews::Fields& cfg,int variable,double lo,double hi,const char* title){
  using namespace IterativeHeliumViews;
  gPad->SetLeftMargin(Number(cfg,"summaryMarginLeft"));gPad->SetRightMargin(Number(cfg,"summaryMarginRight"));
  gPad->SetTopMargin(Number(cfg,"summaryMarginTop"));gPad->SetBottomMargin(Number(cfg,"summaryMarginBottom"));
  double low=std::numeric_limits<double>::max(),high=-low;
  const int colors[]={kBlack,kBlue+1,kOrange+7,kRed+1};
  std::map<int,TGraph*> graphs;std::map<int,std::string> labels;
  for(unsigned i=0;i<rows.size();++i){const Row& r=rows[i];if(r.r<lo||r.r>hi)continue;
    double y=variable==0?r.mean3:variable==1?r.sigma3:variable==2?r.mean4:variable==3?r.sigma4:
      variable==4?r.deviance:(r.ndf>0?r.chi2/r.ndf:0.);
    low=std::min(low,y);high=std::max(high,y);
    if(!graphs.count(r.phase)){graphs[r.phase]=new TGraph();labels[r.phase]=r.label;}
    graphs[r.phase]->SetPoint(graphs[r.phase]->GetN(),r.r,y);
  }
  if(graphs.empty())throw std::runtime_error("Empty summary panel");
  double pad=Number(cfg,"summaryRangePadding")*std::max(Number(cfg,"summaryMinimumRange"),high-low);
  TH1D* frame=new TH1D(Form("summary_frame_%d_%g",variable,lo),title,Integer(cfg,"summaryFrameBins"),lo,hi);
  frame->SetDirectory(0);frame->SetMinimum(std::max(0.,low-pad));frame->SetMaximum(high+pad);
  frame->GetXaxis()->SetTitle("R = |p/q| (GeV/c)");frame->GetYaxis()->SetTitle(variable==0||variable==2?"Mean dE/dx (keV/cm)":variable==1||variable==3?"Gaussian width (keV/cm)":variable==4?"Raw data deviance":"Conditional Pearson #chi^{2}/ndf");
  frame->GetXaxis()->SetTitleSize(Number(cfg,"summaryAxisTitleSize"));frame->GetYaxis()->SetTitleSize(Number(cfg,"summaryAxisTitleSize"));
  frame->GetXaxis()->SetLabelSize(Number(cfg,"summaryAxisLabelSize"));frame->GetYaxis()->SetLabelSize(Number(cfg,"summaryAxisLabelSize"));frame->GetYaxis()->SetTitleOffset(Number(cfg,"summaryYTitleOffset"));frame->Draw();
  TLegend* legend=new TLegend(Number(cfg,"summaryLegendXMin"),Number(cfg,"summaryLegendYMin"),Number(cfg,"summaryLegendXMax"),Number(cfg,"summaryLegendYMax"));legend->SetBorderSize(0);legend->SetFillStyle(0);legend->SetTextSize(Number(cfg,"summaryLegendTextSize"));
  for(std::map<int,TGraph*>::iterator i=graphs.begin();i!=graphs.end();++i){TGraph* graph=i->second;graph->SetLineColor(colors[i->first%4]);graph->SetMarkerColor(colors[i->first%4]);graph->SetLineWidth(2);graph->SetMarkerStyle(20+i->first);graph->SetMarkerSize(Number(cfg,"summaryMarkerSize"));graph->Draw(i->first==0?"P SAME":"LP SAME");legend->AddEntry(graph,labels[i->first].c_str(),i->first==0?"p":"lp");}legend->Draw();
  TLatex note;note.SetNDC();note.SetTextSize(Number(cfg,"summaryNoteSize"));note.SetTextColor(kRed+1);note.DrawLatex(.14,.94,"Study only: curve smoothness is NOT isotope validation");
  gPad->Modified();gPad->Update();
  // ROOT pads hold pointers; delete only after their canvas page has printed.
}
}
#endif
void DrawIterativeSummary(const char* table,const char* settings,const char* output){
#ifndef __CINT__
  using namespace IterativeSummary;using namespace IterativeHeliumViews;
  try{
    if(!gSystem->AccessPathName(output))throw std::runtime_error("Existing summary PDF must be preserved");
    const Fields cfg=ReadSettings(settings);std::ifstream input(table);std::vector<Row> rows;std::string line;
    while(std::getline(input,line)){Row r;std::istringstream in(line);if(!(in>>r.phase>>r.label>>r.index>>r.r>>r.mean3>>r.sigma3>>r.mean4>>r.sigma4>>r.deviance>>r.chi2>>r.ndf>>r.valid))throw std::runtime_error("Bad summary row");rows.push_back(r);}
    if(!input.eof()||rows.empty())throw std::runtime_error("Cannot read summary table");
    gStyle->SetOptStat(0);TCanvas canvas("iterative_summary","Iterative helium calibration",Integer(cfg,"summaryCanvasWidth"),Integer(cfg,"summaryCanvasHeight"));
    canvas.Print(TString(output)+"[");const char* titles[]={"^{3}He component mean","^{3}He component width","^{4}He component mean","^{4}He component width"};
    for(int page=0;page<2;++page){canvas.Clear();canvas.Divide(2,2);double lo=Number(cfg,page==0?"summaryRigidityMin":"summaryLowRigidityMin"),hi=Number(cfg,page==0?"summaryRigidityMax":"summaryLowRigidityMax");for(int panel=0;panel<4;++panel){canvas.cd(panel+1);Panel(rows,cfg,panel,lo,hi,titles[panel]);}canvas.Print(output);}
    canvas.Clear();canvas.Divide(2,1);for(int panel=0;panel<2;++panel){canvas.cd(panel+1);Panel(rows,cfg,4+panel,Number(cfg,"summaryRigidityMin"),Number(cfg,"summaryRigidityMax"),panel==0?"Unpenalized data fit comparison":"Conditional residual diagnostic (not global fit NDF)");}canvas.Print(output);canvas.Print(TString(output)+"]");
    std::cout<<"ITERATIVE_SUMMARY_COMPLETE"<<std::endl;
  }catch(const std::exception& e){std::cerr<<e.what()<<std::endl;gSystem->Exit(2);}
#endif
}
