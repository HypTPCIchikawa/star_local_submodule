void DrawSharedPidFits(const char*,const char*,const char*,const char*);
#ifndef __CINT__
#include "DrawPidRecalibration.C"
#include <TH1.h>
#include <TPad.h>
#include <cstdlib>
namespace SharedPidFitView {
struct Row {unsigned index;std::string id,target,mode,numerical,reason;double r,lo,hi,mean,sigma,chi;};
std::vector<std::string> Split(const std::string& s){std::vector<std::string> out;size_t a=0,b;while((b=s.find('\t',a))!=std::string::npos){out.push_back(s.substr(a,b-a));a=b+1;}out.push_back(s.substr(a));return out;}
double Num(const std::string& s){char* end=0;double n=std::strtod(s.c_str(),&end);if(end==s.c_str()||!end||*end||!TMath::Finite(n))throw std::runtime_error("Invalid metrics number "+s);return n;}
std::vector<Row> Read(const char* path){std::ifstream in(path);if(!in)throw std::runtime_error("Cannot read metrics");std::string line;std::vector<Row> rows;std::map<unsigned,bool> seen;
  while(std::getline(in,line)){if(line.empty()||line[0]=='#')continue;std::vector<std::string> f=Split(line);if(f[0]=="index")continue;if(f.size()!=12)throw std::runtime_error("Expected 12 metrics columns");Row q;double ix=Num(f[0]);if(ix<0||ix!=TMath::Floor(ix))throw std::runtime_error("Invalid index");q.index=unsigned(ix);if(seen[q.index])throw std::runtime_error("Duplicate metrics index");seen[q.index]=true;q.id=f[1];q.target=f[2];q.mode=f[3];q.r=Num(f[4]);q.lo=Num(f[5]);q.hi=Num(f[6]);q.mean=Num(f[7]);q.sigma=Num(f[8]);q.chi=Num(f[9]);q.numerical=f[10];q.reason=f[11];if(!(q.lo<q.hi&&q.sigma>0))throw std::runtime_error("Invalid ROI");rows.push_back(q);}if(rows.empty())throw std::runtime_error("No metrics");return rows;}
TH1* Hist(TFile& file,const char* prefix,unsigned index){TH1* h=dynamic_cast<TH1*>(file.Get(Form("%s_%04u",prefix,index)));if(!h||h->GetDimension()!=1)throw std::runtime_error("Missing saved histogram");return h;}
void SameAxis(const TH1* a,const TH1* b){if(a->GetNbinsX()!=b->GetNbinsX())throw std::runtime_error("Histogram axes differ");for(int k=1;k<=a->GetNbinsX()+1;++k)if(a->GetXaxis()->GetBinLowEdge(k)!=b->GetXaxis()->GetBinLowEdge(k))throw std::runtime_error("Histogram edges differ");}
TGraph* Masked(const TH1* h,const TH1* mask,int color,int style,int width,std::vector<TGraph*>& keep){TGraph* first=0;int start=0;
  for(int b=1;b<=h->GetNbinsX()+1;++b){bool on=b<=h->GetNbinsX()&&mask->GetBinContent(b)==1.;if(on&&!start)start=b;if(!on&&start){TGraph* g=new TGraph(b-start);for(int k=start;k<b;++k)g->SetPoint(k-start,h->GetXaxis()->GetBinCenter(k),h->GetBinContent(k));g->SetLineColor(color);g->SetLineStyle(style);g->SetLineWidth(width);g->Draw("L SAME");keep.push_back(g);if(!first)first=g;start=0;}}
  if(!first){first=new TGraph();first->SetLineColor(color);first->SetLineStyle(style);keep.push_back(first);}return first;}
void Panel(TPad* pad,TFile& file,const Row& q,bool log,double padding,double headroom,double logMin,std::vector<TH1*>& hs,std::vector<TGraph*>& gs){
  pad->cd();pad->SetLeftMargin(.12);pad->SetRightMargin(.035);pad->SetTopMargin(.29);pad->SetBottomMargin(.18);pad->SetLogy(log);
  TH1* raw=Hist(file,"raw_native",q.index);TH1* mask=Hist(file,"fit_mask",q.index);SameAxis(raw,mask);const char* prefixes[]={"prediction_total","prediction_smooth","prediction_d","prediction_t","prediction_he3","prediction_he4"};const char* labels[]={"Total fit","Smooth BG","d","t","^{3}He","^{4}He"};int colors[]={kRed+1,kGray+2,kBlue+1,kOrange+7,kGreen+2,kMagenta+1};TH1* pred[6];for(unsigned j=0;j<6;++j){pred[j]=Hist(file,prefixes[j],q.index);SameAxis(raw,pred[j]);}
  double lo=std::max(raw->GetXaxis()->GetXmin(),q.lo-padding*q.sigma),hi=std::min(raw->GetXaxis()->GetXmax(),q.hi+padding*q.sigma),largest=0.;
  for(int b=1;b<=raw->GetNbinsX();++b)if(raw->GetXaxis()->GetBinCenter(b)>=lo&&raw->GetXaxis()->GetBinCenter(b)<=hi){largest=std::max(largest,raw->GetBinContent(b)+raw->GetBinError(b));for(unsigned j=0;j<6;++j)if(mask->GetBinContent(b)==1.)largest=std::max(largest,pred[j]->GetBinContent(b));}
  if(!(largest>0))throw std::runtime_error("Empty local view");TH1* h=dynamic_cast<TH1*>(raw->Clone(Form("sharedPidView_%u_%d",q.index,int(log))));h->SetDirectory(0);hs.push_back(h);h->SetTitle("");h->GetXaxis()->SetTitle("dE/dx [keV/cm]");h->GetYaxis()->SetTitle(Form("Counts / %.4g keV/cm",h->GetXaxis()->GetBinWidth(1)));h->GetXaxis()->SetRangeUser(lo,hi);h->SetMinimum(log?logMin:0.);h->SetMaximum(largest*(1.+headroom));h->SetMarkerStyle(20);h->SetMarkerSize(.38);h->SetLineColor(kBlack);h->Draw("E1");
  TLegend legend(.52,.735,.97,.975);legend.SetBorderSize(0);legend.SetFillStyle(0);legend.SetTextSize(.029);legend.SetNColumns(2);legend.AddEntry(h,"Raw counts","lep");
  for(unsigned j=0;j<6;++j){TGraph* g=Masked(pred[j],mask,colors[j],j==1?2:1,j==0?2:1,gs);legend.AddEntry(g,labels[j],"l");}legend.DrawClone();
  using PidRecalibrationView::Text;Text(.12,.967,Form("%s; R=%.4g GeV/c",q.target.c_str(),q.r),.036);Text(.12,.915,Form("%s; %s",q.mode.c_str(),log?"log":"linear"),.031);Text(.12,.865,Form("#chi^{2}/ndf=%.4g; quality=%s",q.chi,q.numerical.c_str()),.028);Text(.12,.815,Form("Fit [%.4g,%.4g] keV/cm",q.lo,q.hi),.027);
  std::string note=q.reason;std::replace(note.begin(),note.end(),'\n',' ');if(note.size()>156)note=note.substr(0,153)+"...";Text(.12,.100,note.substr(0,78).c_str(),.024);if(note.size()>78)Text(.12,.065,note.substr(78).c_str(),.024);Text(.12,.027,"Fixed common curves; native counts; no refit / rebin / normalization.",.023);pad->RedrawAxis();
}
void Run(const char* root,const char* metrics,const char* settings,const char* pdf){using PidRecalibrationView::Number;std::map<std::string,std::string> c=PidRecalibrationView::Settings(settings);if(Number(c,"displayRebin")!=1.)throw std::runtime_error("Native binning required");double padding=Number(c,"viewPaddingSigma"),head=Number(c,"linearPadFraction"),logMin=Number(c,"logMin"),cp=Number(c,"maxCasesPerPage");if(padding<0||head<=0||logMin<=0||(cp!=1&&cp!=2))throw std::runtime_error("Invalid view settings");unsigned perPage=unsigned(cp);std::vector<Row> rows=Read(metrics);TFile file(root,"READ");if(file.IsZombie())throw std::runtime_error("Cannot open fit ROOT");gStyle->SetOptStat(0);TCanvas canvas("sharedPidFitCanvas","Common-curve fit diagnostics",1500,perPage==1?850:1250);canvas.Print((std::string(pdf)+"[").c_str());unsigned pages=0;
  for(unsigned i=0;i<rows.size();i+=perPage){canvas.Clear();canvas.Divide(2,perPage,.002,.003);std::vector<TH1*> hs;std::vector<TGraph*> gs;for(unsigned j=0;j<perPage&&i+j<rows.size();++j){Panel(dynamic_cast<TPad*>(canvas.cd(2*j+1)),file,rows[i+j],false,padding,head,logMin,hs,gs);Panel(dynamic_cast<TPad*>(canvas.cd(2*j+2)),file,rows[i+j],true,padding,head,logMin,hs,gs);}canvas.Print(pdf);++pages;canvas.Clear();for(unsigned j=0;j<gs.size();++j)delete gs[j];for(unsigned j=0;j<hs.size();++j)delete hs[j];}
  canvas.Print((std::string(pdf)+"]").c_str());std::cout<<"SHARED_PID_FIT_DRAW_COMPLETE rows="<<rows.size()<<" pages="<<pages<<" rawCountsUnchanged=1 maskApplied=1"<<std::endl;
}
}
#endif
void DrawSharedPidFits(const char* root,const char* metrics,const char* settings,const char* pdf){
#ifndef __CINT__
  try{SharedPidFitView::Run(root,metrics,settings,pdf);}catch(const std::exception& e){std::cerr<<"SHARED_PID_FIT_DRAW_FATAL: "<<e.what()<<std::endl;gSystem->Exit(1);}
#endif
}
