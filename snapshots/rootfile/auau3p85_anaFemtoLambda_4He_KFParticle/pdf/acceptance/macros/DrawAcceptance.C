// Read-only postprocessing of saved candidate-acceptance histograms (ROOT 5).
#include "TCanvas.h"
#include "TFile.h"
#include "TH2D.h"
#include "TLatex.h"
#include "TMath.h"
#include "TPad.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace KFAcceptance {
struct Item {
  std::string id, histogram, title, selection;
  TH2D *source;
  int firstX, lastX, lastY;
  double entries, regular, flow, minimumPositive, maximum, displayed;
  Item():source(0),firstX(0),lastX(0),lastY(0),entries(0),regular(0),flow(0),
         minimumPositive(0),maximum(0),displayed(0) {}
};

std::string JsonString(const std::string &s) {
  std::ostringstream out; out << '"';
  for(std::string::const_iterator p=s.begin();p!=s.end();++p) {
    const unsigned char c=static_cast<unsigned char>(*p);
    if(c=='"' || c=='\\') out << '\\' << *p;
    else if(c=='\n') out << "\\n";
    else if(c=='\r') out << "\\r";
    else if(c=='\t') out << "\\t";
    else if(c<32) { out << "\\u" << std::hex << std::setw(4) << std::setfill('0') << int(c) << std::dec; }
    else out << *p;
  }
  out << '"'; return out.str();
}
void Require(bool okay,const std::string &message) {
  if(!okay) throw std::runtime_error(message);
}
bool SafeId(const std::string &s) {
  if(s.empty()) return false;
  for(unsigned i=0;i<s.size();++i) {
    const char c=s[i];
    if(!((c>='a' && c<='z') || (c>='A' && c<='Z') ||
         (i>0 && ((c>='0' && c<='9') || c=='_')))) return false;
  }
  return true;
}
std::vector<std::string> Fields(const std::string &line) {
  std::vector<std::string> result;
  std::string::size_type first=0;
  while(true) {
    const std::string::size_type end=line.find('\t',first);
    result.push_back(line.substr(first,end==std::string::npos?end:end-first));
    if(end==std::string::npos) break;
    first=end+1;
  }
  return result;
}
std::vector<KFAcceptance::Item> ReadManifest(const char *path) {
  std::ifstream in(path); Require(bool(in),"Cannot open TSV manifest");
  std::string line; Require(bool(std::getline(in,line)),"Empty TSV manifest");
  if(!line.empty() && line[line.size()-1]=='\r') line.erase(line.size()-1);
  Require(line=="id\thistogram\ttitle\tselection","Unexpected TSV manifest header");
  std::vector<KFAcceptance::Item> items; std::set<std::string> ids, histograms;
  while(std::getline(in,line)) {
    if(!line.empty() && line[line.size()-1]=='\r') line.erase(line.size()-1);
    if(line.empty()) continue;
    const std::vector<std::string> fields=Fields(line);
    Require(fields.size()==4,"Expected four TSV fields per row");
    Require(SafeId(fields[0]),"Invalid output id: "+fields[0]);
    Require(fields[0]!="overview" && fields[0]!="acceptance_pt_ylab","Reserved output id: "+fields[0]);
    Require(!fields[1].empty() && !fields[2].empty() && !fields[3].empty(),"Empty required TSV field");
    Require(ids.insert(fields[0]).second && histograms.insert(fields[1]).second,"Duplicate TSV id or histogram");
    KFAcceptance::Item item;
    item.id=fields[0]; item.histogram=fields[1]; item.title=fields[2]; item.selection=fields[3];
    items.push_back(item);
  }
  Require(items.size()==4,"Manifest must list exactly four species histograms");
  return items;
}
bool SameX(const TH2D *left,const TH2D *right) {
  if(left->GetNbinsX()!=right->GetNbinsX()) return false;
  for(int b=1;b<=left->GetNbinsX()+1;++b) {
    const double a=left->GetXaxis()->GetBinLowEdge(b), c=right->GetXaxis()->GetBinLowEdge(b);
    if(!TMath::Finite(a) || !TMath::Finite(c) || std::fabs(a-c)>1.e-12*std::max(1.,std::max(std::fabs(a),std::fabs(c)))) return false;
  }
  return true;
}
void Inspect(KFAcceptance::Item &item,TFile *input) {
  item.source=dynamic_cast<TH2D*>(input->Get(item.histogram.c_str()));
  Require(item.source!=0,"Missing TH2D: "+item.histogram);
  TH2D *h=item.source;
  Require(h->GetNbinsX()>0 && h->GetNbinsY()>0,"Empty histogram axis: "+item.histogram);
  for(int axis=0;axis<2;++axis) {
    const TAxis *a=axis==0?h->GetXaxis():h->GetYaxis();
    for(int b=1;b<=a->GetNbins()+1;++b) {
      Require(TMath::Finite(a->GetBinLowEdge(b)),"Nonfinite axis boundary: "+item.histogram);
      if(b>1) Require(a->GetBinLowEdge(b)>a->GetBinLowEdge(b-1),"Non-increasing axis: "+item.histogram);
    }
  }
  Require(std::fabs(h->GetYaxis()->GetXmin())<1.e-12,"Expected pT axis starting at zero: "+item.histogram);
  item.entries=h->GetEntries();
  Require(TMath::Finite(item.entries) && item.entries>0,"Nonpositive/invalid Entries: "+item.histogram);
  item.firstX=h->GetNbinsX()+1; item.lastX=0; item.lastY=0;
  item.minimumPositive=std::numeric_limits<double>::max();
  for(int x=0;x<=h->GetNbinsX()+1;++x) for(int y=0;y<=h->GetNbinsY()+1;++y) {
    const double n=h->GetBinContent(x,y), e=h->GetBinError(x,y);
    Require(TMath::Finite(n) && TMath::Finite(e) && n>=0 && e>=0,"Invalid or negative count/error: "+item.histogram);
    const bool regular=x>=1 && x<=h->GetNbinsX() && y>=1 && y<=h->GetNbinsY();
    if(regular) {
      item.regular+=n;
      if(n>0) {
        item.firstX=std::min(item.firstX,x); item.lastX=std::max(item.lastX,x);
        item.lastY=std::max(item.lastY,y);
        item.minimumPositive=std::min(item.minimumPositive,n); item.maximum=std::max(item.maximum,n);
      }
    } else item.flow+=n;
  }
  Require(TMath::Finite(item.regular) && TMath::Finite(item.flow) && item.regular>0,
          "No positive regular-bin counts or nonfinite sum: "+item.histogram);
}
void Header(TCanvas *canvas,const char *energy,const char *sample) {
  canvas->cd(); TLatex text; text.SetNDC(); text.SetTextFont(42);
  text.SetTextSize(.026); text.DrawLatex(.04,.966,energy);
  text.SetTextSize(.019); text.DrawLatex(.04,.934,sample);
  text.SetTextSize(.017); text.DrawLatex(.04,.044,"Pre-pair candidate counts, not efficiency; laboratory rapidity, no CM shift");
  text.DrawLatex(.04,.020,"Inclusive accepted events (no centrality slicing); flow bins recorded in audit, not reassigned");
}
void Panel(const KFAcceptance::Item &item,int firstX,int lastX,int lastY,
           const std::string &unique,std::vector<TH2D*> &drawn,bool overview) {
  TPad *pad=dynamic_cast<TPad*>(gPad); Require(pad!=0,"No current drawing pad");
  pad->SetLogz(1); pad->SetLeftMargin(.14); pad->SetRightMargin(.18);
  pad->SetBottomMargin(.15); pad->SetTopMargin(.23);
  TH2D *h=dynamic_cast<TH2D*>(item.source->Clone(unique.c_str()));
  Require(h!=0,"Cannot clone source histogram"); h->SetDirectory(0); drawn.push_back(h);
  h->SetTitle(""); h->SetStats(false);
  h->GetXaxis()->SetRange(firstX,lastX); h->GetYaxis()->SetRange(1,lastY);
  h->GetXaxis()->SetTitle("Y_{lab}"); h->GetYaxis()->SetTitle("p_{T} [GeV/c]");
  h->GetZaxis()->SetTitle("Candidates / bin");
  const double titleSize=overview?.053:.044, labelSize=overview?.045:.037;
  h->GetXaxis()->SetTitleSize(titleSize); h->GetYaxis()->SetTitleSize(titleSize); h->GetZaxis()->SetTitleSize(titleSize*.87);
  h->GetXaxis()->SetLabelSize(labelSize); h->GetYaxis()->SetLabelSize(labelSize); h->GetZaxis()->SetLabelSize(labelSize*.80);
  h->GetYaxis()->SetTitleOffset(1.05); h->GetZaxis()->SetTitleOffset(1.13);
  h->SetMinimum(item.minimumPositive); h->SetMaximum(item.maximum>item.minimumPositive?item.maximum:item.maximum*2.);
  h->Draw("COLZ");
  TLatex label; label.SetNDC(); label.SetTextFont(42);
  label.SetTextSize(overview?.044:.039); label.DrawLatex(.14,.948,item.title.c_str());
  label.SetTextSize(overview?.029:.026); label.DrawLatex(.14,.886,item.histogram.c_str());
  label.SetTextSize(overview?.027:.024); label.DrawLatex(.14,.828,item.selection.c_str());
}

void Run(const char *inputPath,const char *outPath,const char *manifestPath,
         const char *energy,const char *sample,int padding) {
  Require(inputPath && *inputPath && outPath && *outPath && manifestPath && *manifestPath &&
          energy && *energy && sample && *sample,"All path and label arguments are required");
  Require(padding>=0,"paddingBins must be nonnegative");
  std::vector<KFAcceptance::Item> items=ReadManifest(manifestPath);
  const std::string out(outPath);
  std::vector<std::string> outputs;
  outputs.push_back(out+"/acceptance_pt_ylab.pdf"); outputs.push_back(out+"/overview.png");
  outputs.push_back(out+"/acceptance_canvases.root"); outputs.push_back(out+"/acceptance_audit.json");
  for(unsigned i=0;i<items.size();++i) { outputs.push_back(out+"/"+items[i].id+".pdf"); outputs.push_back(out+"/"+items[i].id+".png"); }
  for(unsigned i=0;i<outputs.size();++i)
    Require(gSystem->AccessPathName(outputs[i].c_str()),"Refusing to overwrite: "+outputs[i]);
  TFile *input=TFile::Open(inputPath,"READ");
  Require(input && !input->IsZombie(),"Cannot open input ROOT read-only");
  int firstX=0,lastX=0;
  for(unsigned i=0;i<items.size();++i) {
    Inspect(items[i],input);
    if(i==0) {firstX=items[i].firstX;lastX=items[i].lastX;}
    else { Require(SameX(items[0].source,items[i].source),"YLab axes do not match");
      firstX=std::min(firstX,items[i].firstX); lastX=std::max(lastX,items[i].lastX); }
  }
  const int nx=items[0].source->GetNbinsX();
  firstX=std::max(1,firstX-std::min(padding,nx)); lastX=std::min(nx,lastX+std::min(padding,nx));
  for(unsigned i=0;i<items.size();++i) {
    const int ny=items[i].source->GetNbinsY();
    items[i].lastY=std::min(ny,items[i].lastY+std::min(padding,ny));
    items[i].displayed=items[i].source->Integral(firstX,lastX,1,items[i].lastY);
    Require(std::fabs(items[i].displayed-items[i].regular)<=1.e-10*std::max(1.,items[i].regular),
            "Display would lose regular-bin counts: "+items[i].histogram);
  }
  if(gSystem->AccessPathName(out.c_str())) Require(gSystem->mkdir(out.c_str(),kTRUE)==0,"Cannot create output directory");
  TFile *canvasFile=TFile::Open((out+"/acceptance_canvases.root").c_str(),"NEW");
  Require(canvasFile && !canvasFile->IsZombie(),"Cannot create fresh canvas ROOT");
  gStyle->SetOptStat(0); gStyle->SetPalette(1); gStyle->SetNumberContours(60);
  std::vector<TCanvas*> canvases; std::vector<TH2D*> drawn;
  TCanvas *overview=new TCanvas("acceptance_overview","Candidate acceptance overview",1500,1200);
  canvases.push_back(overview); Header(overview,energy,sample);
  TPad *grid=new TPad("acceptance_grid","",0,.078,1,.913); grid->Draw(); grid->Divide(2,2,0.004,0.004);
  for(unsigned i=0;i<items.size();++i) {
    grid->cd(i+1); Panel(items[i],firstX,lastX,items[i].lastY,"acceptance_overview_"+items[i].id,drawn,true);
  }
  const std::string book=out+"/acceptance_pt_ylab.pdf";
  overview->Print((book+"[").c_str()); overview->Print(book.c_str());
  overview->Print((out+"/overview.png").c_str()); canvasFile->cd();
  Require(overview->Write()>0,"Failed to write overview canvas");
  for(unsigned i=0;i<items.size();++i) {
    TCanvas *canvas=new TCanvas(("acceptance_"+items[i].id).c_str(),items[i].title.c_str(),1200,1000);
    canvases.push_back(canvas); Header(canvas,energy,sample);
    TPad *panel=new TPad(("acceptance_panel_"+items[i].id).c_str(),"",0,.078,1,.913); panel->Draw(); panel->cd();
    Panel(items[i],firstX,lastX,items[i].lastY,"acceptance_single_"+items[i].id,drawn,false);
    canvas->Print(book.c_str()); canvas->Print((out+"/"+items[i].id+".pdf").c_str());
    canvas->Print((out+"/"+items[i].id+".png").c_str()); canvasFile->cd();
    Require(canvas->Write()>0,"Failed to write individual canvas");
  }
  overview->Print((book+"]").c_str()); canvasFile->Close(); delete canvasFile;

  std::ofstream audit((out+"/acceptance_audit.json").c_str()); Require(bool(audit),"Cannot create audit JSON");
  audit<<std::setprecision(17)<<"{\n\"schemaVersion\":1,\"inputRoot\":"<<JsonString(inputPath)
    <<",\"manifest\":"<<JsonString(manifestPath)<<",\"energy_label\":"<<JsonString(energy)
    <<",\"sample_label\":"<<JsonString(sample)<<",\"padding_bins\":"<<padding
    <<",\"centrality_selection\":\"inclusive stored accepted events; no centrality slicing\""
    <<",\"caveat\":\"Pre-pair candidate counts, not efficiency. Laboratory rapidity, no CM shift. Daughter entries may repeat across Lambda candidates. Flow bins are recorded separately and not reassigned.\""
    <<",\"normalization\":\"none\",\"rebin\":false,\"book_pages\":5,\n\"histograms\":[\n";
  for(unsigned i=0;i<items.size();++i) {
    if(i) audit<<",\n"; const KFAcceptance::Item &v=items[i]; const TH2D *h=v.source;
    audit<<"{\"id\":"<<JsonString(v.id)<<",\"histogram\":"<<JsonString(v.histogram)
      <<",\"title\":"<<JsonString(v.title)<<",\"selection\":"<<JsonString(v.selection)
      <<",\"entries\":"<<v.entries<<",\"regular_count\":"<<v.regular<<",\"flow_count\":"<<v.flow
      <<",\"x_bins\":"<<h->GetNbinsX()<<",\"y_bins\":"<<h->GetNbinsY()
      <<",\"x_min\":"<<h->GetXaxis()->GetXmin()<<",\"x_max\":"<<h->GetXaxis()->GetXmax()
      <<",\"y_min\":"<<h->GetYaxis()->GetXmin()<<",\"y_max\":"<<h->GetYaxis()->GetXmax()
      <<",\"display_x\":["<<h->GetXaxis()->GetBinLowEdge(firstX)<<','<<h->GetXaxis()->GetBinUpEdge(lastX)
      <<"],\"display_y\":["<<h->GetYaxis()->GetBinLowEdge(1)<<','<<h->GetYaxis()->GetBinUpEdge(v.lastY)
      <<"],\"display_count\":"<<v.displayed<<",\"minimum_positive\":"<<v.minimumPositive
      <<",\"maximum\":"<<v.maximum<<'}';
  }
  audit<<"\n]}\n"; audit.close(); Require(bool(audit),"Failed writing audit JSON");
  for(unsigned i=0;i<outputs.size();++i) Require(!gSystem->AccessPathName(outputs[i].c_str()),"Expected output missing: "+outputs[i]);
  for(unsigned i=0;i<canvases.size();++i) delete canvases[i];
  for(unsigned i=0;i<drawn.size();++i) delete drawn[i];
  input->Close(); delete input;
  std::cout<<"DRAW_ACCEPTANCE_OK "<<book<<" pages=5 histograms=4"<<std::endl;
}
}

void DrawAcceptance(const char *input,const char *out,const char *manifest,
                    const char *energyLabel,const char *sampleLabel,int paddingBins) {
  try { KFAcceptance::Run(input,out,manifest,energyLabel,sampleLabel,paddingBins); }
  catch(const std::exception &error) {
    std::cerr<<"DrawAcceptance ERROR: "<<error.what()<<std::endl; gSystem->Exit(2);
  }
}
