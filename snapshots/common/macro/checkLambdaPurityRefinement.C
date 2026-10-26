// ROOT5 saved-candidate QA. New full-data refinement; no fresh holdout claim.
// Reuses read-only CSV/plot helpers without modifying the original study.
// checkLambdaPurityRefinement(scanDirectory, outputPrefix)
// Manifest background keys: backgroundDegrees (0;1;2), backgroundRegions
// (name:leftLow:leftHigh:rightLow:rightHigh;...), backgroundProfileDeltaNll,
// backgroundProfileMaxFactor (numerical bracket cap), backgroundMinCounts
// (warning only). All mass bounds and background models are manifest inputs.
// Positive Bernstein coefficients define a restricted smooth-background family,
// not a guarantee against peaking backgrounds or unmodeled signal tails.
#include "drawLambdaPurityStudy.C"
#include <TBox.h>

Int_t checkLambdaPurityRefinement(const char*, const char*);
#ifndef __CINT__
namespace LambdaPurityRefinement {
using LambdaPurityQA::Row;
using LambdaPurityQA::Choice;
using LambdaPurityQA::Require;
using LambdaPurityQA::Get;
using LambdaPurityQA::N;
using LambdaPurityQA::Number;
using LambdaPurityQA::ReadCsv;
using LambdaPurityQA::CsvQuote;
using LambdaPurityQA::Edge;
using LambdaPurityQA::Color;
using LambdaPurityQA::Marker;
using LambdaPurityQA::HistStyle;
using LambdaPurityQA::Axis;
using LambdaPurityQA::Save;
typedef std::vector<double> Vector;
typedef std::vector<Vector> Matrix;
struct Region { std::string name; double ll,lh,rl,rh; };
struct Config {
  double lo,hi,bw,displayLo,displayHi,wl,wh,delta,maxFactor,minCounts;
  std::string title,raw;
  std::vector<int> degrees;
  std::vector<Region> regions;
};
struct Opt {
  Vector weights; double nll,kkt; int iterations; bool converged;
  Opt():nll(0),kkt(0),iterations(0),converged(false){}
};
struct Result {
  int degree,iterations; bool valid,boundary,profileValid;
  double side,nwin,b,purity,nll,kkt,upper,zeroUpper,deviance;
  Vector coefficients;
  std::string flags;
  Result():degree(0),iterations(0),valid(false),boundary(false),profileValid(false),
    side(0),nwin(0),b(0),purity(0),nll(0),kkt(0),upper(0),zeroUpper(0),deviance(0){}
};
static std::vector<std::string> Split(const std::string&s,char sep){std::vector<std::string>v;std::stringstream in(s);std::string x;while(std::getline(in,x,sep))v.push_back(x);return v;}
static void Flag(Result&r,const std::string&s){if(!r.flags.empty())r.flags+=';';r.flags+=s;}
static Config Manifest(const std::string&path){
  std::ifstream in(path.c_str());Require(bool(in),"cannot read manifest");Row m;Config c;std::string line;
  while(std::getline(in,line)){c.raw+=line+"\n";if(line.empty()||line[0]=='#')continue;size_t t=line.find('\t');Require(t!=std::string::npos,"manifest requires TAB separator");std::string v=line.substr(t+1);if(!v.empty()&&v[v.size()-1]=='\r')v.erase(v.size()-1);Require(m.insert(std::make_pair(line.substr(0,t),v)).second,"duplicate manifest key");}
  c.lo=N(m,"massMin");c.hi=N(m,"massMax");c.bw=N(m,"massBinWidth");c.displayLo=N(m,"displayMin");c.displayHi=N(m,"displayMax");c.wl=N(m,"windowLow");c.wh=N(m,"windowHigh");c.title=Get(m,"title");
  c.delta=N(m,"backgroundProfileDeltaNll");c.maxFactor=N(m,"backgroundProfileMaxFactor");c.minCounts=N(m,"backgroundMinCounts");
  Require(c.lo<c.wl&&c.wl<c.wh&&c.wh<c.hi&&c.bw>0&&c.displayLo>=c.lo&&c.displayHi<=c.hi&&c.displayLo<c.displayHi,"invalid mass ranges");Require(c.delta>0&&c.maxFactor>1&&c.minCounts>=0,"invalid background numerical settings");
  std::vector<std::string>d=Split(Get(m,"backgroundDegrees"),';');for(size_t k=0;k<d.size();++k){double x=Number(d[k]);Require(x==int(x)&&x>=0&&x<=2,"background degree must be 0,1,2");Require(std::find(c.degrees.begin(),c.degrees.end(),int(x))==c.degrees.end(),"duplicate degree");c.degrees.push_back(int(x));}
  std::vector<std::string>groups=Split(Get(m,"backgroundRegions"),';');
  for(size_t k=0;k<groups.size();++k){std::vector<std::string>f=Split(groups[k],':');Require(f.size()==5,"backgroundRegions schema mismatch");Region r;r.name=f[0];r.ll=Number(f[1]);r.lh=Number(f[2]);r.rl=Number(f[3]);r.rh=Number(f[4]);Require(!r.name.empty()&&r.ll>=c.lo&&r.ll<r.lh&&r.lh<=c.wl&&c.wh<=r.rl&&r.rl<r.rh&&r.rh<=c.hi,"background regions must be outside fixed signal window");c.regions.push_back(r);}
  Require(!c.degrees.empty()&&!c.regions.empty(),"empty background model list");return c;
}
static double Primitive(int degree,int k,double t){
  if(degree==0)return t;
  if(degree==1)return k==0?t-.5*t*t:.5*t*t;
  if(k==0)return t-t*t+t*t*t/3.;
  if(k==1)return t*t-2.*t*t*t/3.;
  return t*t*t/3.;
}
static double Integral(int degree,int k,double low,double high,const Region&r){
  const double span=r.rh-r.ll;return span*(Primitive(degree,k,(high-r.ll)/span)-Primitive(degree,k,(low-r.ll)/span));
}
static double Basis(int degree,int k,double t){if(degree==0)return 1.;if(degree==1)return k==0?1-t:t;if(k==0)return (1-t)*(1-t);if(k==1)return 2*t*(1-t);return t*t;}
static double Nll(const Matrix&a,const Vector&n,const Vector&w,double total){
  double value=0;for(size_t i=0;i<n.size();++i){double mu=0;for(size_t k=0;k<w.size();++k)mu+=total*a[i][k]*w[k];if(mu<0||!TMath::Finite(mu))return 1.e300;if(n[i]>0){if(mu<=0)return 1.e300;value+=mu-n[i]*std::log(mu);}else value+=mu;}return value;
}
static double Derivative(const Matrix&a,const Vector&n,const Vector&w,double total,int p,int q,double shift){
  double g=0;for(size_t i=0;i<n.size();++i){double mu=0;for(size_t k=0;k<w.size();++k)mu+=total*a[i][k]*w[k];const double d=total*(a[i][p]-a[i][q]);mu+=shift*d;
    if(n[i]>0&&mu<=0){if(d>0)return -1.e300;if(d<0)return 1.e300;continue;}g+=d*(n[i]>0?1-n[i]/mu:1);}
  return g;
}
// Convex nonnegative mixture: cyclic exact pair transfers on the simplex.
// Each one-dimensional derivative is monotone. KKT residual is checked after
// convergence; boundaries are physical allowed solutions, not Hessian errors.
static Opt Optimize(const Matrix&a,const Vector&n,double total){
  Opt out;const int dim=int(a[0].size());out.weights.assign(dim,1./dim);
  if(dim==1){out.converged=true;out.nll=Nll(a,n,out.weights,total);return out;}
  if(total==0){double sum=0;for(size_t i=0;i<n.size();++i)sum+=n[i];out.converged=(sum==0);out.nll=sum?1.e300:0;return out;}
  for(int iteration=0;iteration<3000;++iteration){double change=0;
    for(int p=0;p<dim;++p)for(int q=p+1;q<dim;++q){double low=-out.weights[p],high=out.weights[q],shift;
      if(Derivative(a,n,out.weights,total,p,q,low)>=0)shift=low;
      else if(Derivative(a,n,out.weights,total,p,q,high)<=0)shift=high;
      else {for(int b=0;b<55;++b){double mid=(low+high)/2;if(Derivative(a,n,out.weights,total,p,q,mid)<0)low=mid;else high=mid;}shift=(low+high)/2;}
      const double pair=out.weights[p]+out.weights[q];out.weights[p]=std::max(0.,std::min(pair,out.weights[p]+shift));out.weights[q]=pair-out.weights[p];change=std::max(change,std::fabs(shift));
    }
    out.iterations=iteration+1;
    if(change<1.e-11){out.converged=true;break;}
  }
  out.nll=Nll(a,n,out.weights,total);Vector gradients(dim,0.);double common=0,weight=0;
  for(size_t i=0;i<n.size();++i){double mu=0;for(int k=0;k<dim;++k)mu+=total*a[i][k]*out.weights[k];if(n[i]>0&&mu<=0){out.converged=false;return out;}for(int k=0;k<dim;++k)gradients[k]+=total*a[i][k]*(n[i]>0?1-n[i]/mu:1);}
  for(int k=0;k<dim;++k)if(out.weights[k]>1.e-8){common+=out.weights[k]*gradients[k];weight+=out.weights[k];}common/=weight;
  for(int k=0;k<dim;++k)out.kkt=std::max(out.kkt,out.weights[k]>1.e-8?std::fabs(gradients[k]-common):std::max(0.,common-gradients[k]));
  if(out.kkt>1.e-6*(1.+total)||!TMath::Finite(out.nll)||out.nll>=1.e299)out.converged=false;return out;
}
static Result Background(TH1D*h,const Region&r,int degree,const Config&c){
  Result result;result.degree=degree;const int dim=degree+1;
  result.nwin=h->Integral(Edge(h,c.wl),Edge(h,c.wh)-1);
  Matrix a;Vector n;Vector ratio(dim,0),window(dim,0);
  for(int k=0;k<dim;++k){window[k]=Integral(degree,k,c.wl,c.wh,r);Require(window[k]>0,"nonpositive basis window integral");}
  for(int b=Edge(h,r.ll);b<Edge(h,r.rh);++b){const double low=h->GetBinLowEdge(b),high=h->GetBinLowEdge(b+1);if(!(high<=r.lh+1.e-10||low>=r.rl-1.e-10))continue;
    Vector row(dim,0);for(int k=0;k<dim;++k){row[k]=Integral(degree,k,low,high,r)/window[k];Require(row[k]>=0,"negative basis integral");ratio[k]+=row[k];}a.push_back(row);n.push_back(h->GetBinContent(b));result.side+=n.back();}
  Require(!a.empty()&&result.nwin>0,"empty sideband bins or signal window");
  Matrix normalized=a;for(size_t i=0;i<a.size();++i)for(int k=0;k<dim;++k){Require(ratio[k]>0,"basis not constrained by sidebands");normalized[i][k]/=ratio[k];}
  Opt main=Optimize(normalized,n,result.side);result.iterations=main.iterations;result.kkt=main.kkt;result.nll=main.nll;result.valid=main.converged;
  if(!result.valid){Flag(result,"simplex_nonconvergence_or_kkt_failure");return result;}
  result.coefficients.resize(dim);for(int k=0;k<dim;++k){result.coefficients[k]=result.side*main.weights[k]/ratio[k];result.b+=result.coefficients[k];if(main.weights[k]<1.e-7)result.boundary=true;}
  result.purity=1-result.b/result.nwin;
  if(result.side<c.minCounts)Flag(result,"low_sideband_statistics");if(result.boundary)Flag(result,"shape_parameter_on_boundary");if(result.b>result.nwin)Flag(result,"negative_signal_estimate");
  if(result.side==0){Flag(result,"zero_counts_shape_unidentified");const double tail=.5*TMath::Erfc(std::sqrt(c.delta));result.zeroUpper=-std::log(tail)/(*std::min_element(ratio.begin(),ratio.end()));}
  for(size_t i=0;i<n.size();++i){double mu=0;for(int k=0;k<dim;++k)mu+=a[i][k]*result.coefficients[k];if(mu>0)result.deviance+=2*(mu-n[i]+(n[i]>0?n[i]*std::log(n[i]/mu):0));else if(n[i]>0){result.valid=false;Flag(result,"nonpositive_expected_bin");return result;}}
  // Profile B=sum(nonnegative window-yield coefficients). For each fixed B,
  // optimize only its shape fractions, including coefficient boundaries.
  Opt atMinimum=Optimize(a,n,result.b);
  if(!atMinimum.converged||std::fabs(atMinimum.nll-result.nll)>1.e-6*(1+result.side)){Flag(result,"profile_minimum_mismatch");return result;}
  const double cap=std::max(1.,result.b)*c.maxFactor,target=result.nll+c.delta;
  double low=result.b,high=std::min(cap,std::max(1.,2*result.b));Opt upper=Optimize(a,n,high);
  while(upper.converged&&upper.nll<target&&high<cap){high=std::min(cap,2*high);upper=Optimize(a,n,high);}
  if(!upper.converged){Flag(result,"profile_upper_fit_failed");return result;}
  if(upper.nll<target){Flag(result,"profile_upper_search_cap");return result;}
  for(int iteration=0;iteration<45;++iteration){double mid=(low+high)/2;Opt trial=Optimize(a,n,mid);if(!trial.converged){Flag(result,"profile_fit_failed");return result;}
    // Convex profile must not go below the fitted global minimum.
    if(trial.nll<result.nll-1.e-6*(1+result.side)){Flag(result,"profile_below_global_minimum");return result;}
    if(trial.nll<target)low=mid;else high=mid;
  }
  result.upper=(low+high)/2;result.profileValid=true;return result;
}
static void Overlay(const std::vector<Choice>&choices,const std::vector<int>&which,
  const std::vector<TH1D*>&hist,const std::map<std::string,Row>&metrics,const Config&c,
  int index,const std::string&stem,TFile&out){
  TCanvas*canvas=new TCanvas(Form("RefinementOverlay%d",index),"Full-data exploratory KF Lambda cuts",1800,1060);std::vector<TH1D*>legendHist;
  for(int panel=0;panel<2;++panel){canvas->cd();TPad*pad=new TPad(Form("RefinementPad%d_%d",index,panel),"",panel*.5,.06,(panel+1)*.5,.72);pad->SetLeftMargin(.15);pad->SetRightMargin(.04);pad->SetBottomMargin(.15);pad->SetTopMargin(.09);pad->SetTicks(1,1);pad->SetGridy();pad->Draw();pad->cd();double maximum=0;std::vector<TH1D*>drawn;
    for(size_t k=0;k<which.size();++k){int j=which[k];TH1D*h=(TH1D*)hist[j]->Clone(Form("RefinementCountsP%d_%d_%d",index,panel,j));h->SetDirectory(0);double norm=h->Integral(1,h->GetNbinsX());Require(norm>0,"empty normalization range");if(panel)h->Scale(1./norm);HistStyle(h,j);h->GetXaxis()->SetRange(Edge(h,c.displayLo),Edge(h,c.displayHi)-1);maximum=std::max(maximum,h->GetMaximum());drawn.push_back(h);if(!panel)legendHist.push_back(h);
      if(panel)Require(std::fabs(h->Integral(1,h->GetNbinsX())-1)<1.e-10,"normalization audit failure");out.cd();h->Write();pad->cd();}
    TH1D*frame=(TH1D*)drawn[0]->Clone(Form("RefinementFrame%d_%d",index,panel));frame->SetDirectory(0);frame->Reset("ICES");frame->SetMinimum(0);frame->SetMaximum(maximum*1.16);Axis(frame,Form(panel?"Fraction per %.3g MeV/c^{2} bin":"Candidates per %.3g MeV/c^{2} bin",1000*c.bw));frame->Draw("AXIS");
    for(size_t j=drawn.size();j>0;--j)drawn[j-1]->Draw("HIST SAME");
    for(size_t j=0;j<drawn.size();++j){TGraph*marks=new TGraph;for(int b=Edge(drawn[j],c.displayLo);b<Edge(drawn[j],c.displayHi);b+=2)if(drawn[j]->GetBinContent(b)>0)marks->SetPoint(marks->GetN(),drawn[j]->GetBinCenter(b),drawn[j]->GetBinContent(b));marks->SetMarkerColor(Color(which[j]));marks->SetMarkerStyle(Marker(which[j]));marks->SetMarkerSize(.9);marks->Draw("P SAME");}
    TLatex t;t.SetNDC();t.SetTextFont(42);t.SetTextSize(.041);t.DrawLatex(.15,.95,panel?"Unit-area normalized shape":"Raw candidate counts");pad->RedrawAxis();
  }
  canvas->cd();TLatex text;text.SetNDC();text.SetTextFont(42);text.SetTextSize(.029);text.DrawLatex(.045,.96,c.title.c_str());TLegend*legend=new TLegend(.045,.795,.98,.938);legend->SetBorderSize(0);legend->SetFillStyle(0);legend->SetTextFont(42);legend->SetTextSize(.024);legend->SetMargin(.05);
  for(size_t k=0;k<which.size();++k){int j=which[k];const Row&m=metrics.find(choices[j].id)->second;legend->AddEntry(legendHist[k],Form("%s: Nwin=%.0f; P=%.2f%%; S/S_{baseline}=%.1f%%",choices[j].label.c_str(),N(m,"Nwin"),100*N(m,"purity"),100*N(m,"signalRetention")),"lp");}legend->Draw();
  text.SetTextSize(.021);text.DrawLatex(.045,.75,"FULL-DATA EXPLORATORY refinement. All events were previously inspected; no new held-out validation.");text.SetTextSize(.019);text.DrawLatex(.07,.025,Form("P: default-sideband estimate. Fixed window [%.3f, %.3f); normalization [%.2f, %.2f). Original %.3g MeV bins; no mass constraint.",c.wl,c.wh,c.lo,c.hi,c.bw*1000));Save(canvas,stem,out);
}
static void Diagnostic(const Choice&choice,TH1D*h,const Config&c,
  const std::vector<Result>&results,const std::string&stem,int index,TFile&out){
  const int cols=int(c.degrees.size()),rows=int(c.regions.size());TCanvas*canvas=new TCanvas(Form("PositiveBackground%d",index),"Positive sideband-only background diagnostics",1800,1300);canvas->Divide(cols,rows);
  for(int r=0;r<rows;++r)for(int d=0;d<cols;++d){const Region&region=c.regions[r];const Result&fit=results[r*cols+d];canvas->cd(r*cols+d+1);gPad->SetLeftMargin(.14);gPad->SetBottomMargin(.14);gPad->SetTopMargin(.22);
    double ymax=1;TGraphErrors*points=new TGraphErrors;for(int b=Edge(h,region.ll);b<Edge(h,region.rh);++b){double lo=h->GetBinLowEdge(b),hi=h->GetBinLowEdge(b+1);if(!(hi<=region.lh+1.e-10||lo>=region.rl-1.e-10))continue;double count=h->GetBinContent(b);points->SetPoint(points->GetN(),h->GetBinCenter(b),count);points->SetPointError(points->GetN()-1,0,std::sqrt(count));ymax=std::max(ymax,count+std::sqrt(count));}
    TGraph*line=new TGraph;if(fit.valid)for(int k=0;k<=160;++k){double x=region.ll+(region.rh-region.ll)*k/160.;double y=0;for(int j=0;j<=fit.degree;++j)y+=fit.coefficients[j]*Basis(fit.degree,j,(x-region.ll)/(region.rh-region.ll))*c.bw/Integral(fit.degree,j,c.wl,c.wh,region);line->SetPoint(line->GetN(),x,y);ymax=std::max(ymax,y);}
    TH1D*frame=new TH1D(Form("BgFrame%d_%d_%d",index,r,d),"",1,region.ll,region.rh);frame->SetDirectory(0);frame->SetStats(false);frame->SetMinimum(0);frame->SetMaximum(ymax*1.15);Axis(frame,Form("Candidates / %.3g MeV/c^{2}",1000*c.bw));frame->GetYaxis()->SetTitleOffset(1.2);frame->Draw("AXIS");
    TBox*mask=new TBox(region.lh,0,region.rl,ymax*1.15);mask->SetFillColor(TColor::GetColor("#eeeeee"));mask->SetFillStyle(1001);mask->SetLineColor(TColor::GetColor("#eeeeee"));mask->Draw();
    line->SetLineColor(Color(d+1));line->SetLineWidth(2);if(fit.valid)line->Draw("L SAME");points->SetMarkerStyle(20);points->SetMarkerSize(.75);points->Draw("P SAME");gPad->RedrawAxis();
    TLatex text;text.SetNDC();text.SetTextFont(42);text.SetTextSize(.032);text.DrawLatex(.14,.95,Form("%s; %s; Bernstein %d",choice.label.c_str(),region.name.c_str(),fit.degree));
    text.SetTextSize(.03);if(fit.valid){text.DrawLatex(.14,.89,Form("SB count=%.0f; B=%.2f; P=%.2f%%",fit.side,fit.b,100*fit.purity));text.DrawLatex(.14,.835,fit.profileValid?Form("Profile upper B=%.2f (nominal, asymptotic)",fit.upper):"Profile upper B unavailable");}else{text.SetTextColor(kRed+1);text.DrawLatex(.14,.89,"NUMERICALLY INVALID; no estimate");}
    text.SetTextSize(.026);text.DrawLatex(.15,.725,"Gray region excluded from fit; peak counts hidden.");
    out.cd();points->Write(Form("SidebandPoints_%d_%d_%d",index,r,d));line->Write(Form("PositiveBackgroundCurve_%d_%d_%d",index,r,d));
  }Save(canvas,stem,out);
}
static void Checks(const std::vector<Choice>&choices,const std::vector<TH1D*>&hist,
  const Config&c,const std::string&stem,TFile&out){
  std::ofstream csv((stem+"_background.csv").c_str());Require(bool(csv),"cannot create background CSV");csv<<std::setprecision(17);
  csv<<"choice_id,split,sideband,degree,numerical_valid,iterations,kkt_residual,sideband_count,Nwin,B,purity,boundary,profile_valid,profile_upper_B,profile_purity_lower,profile_delta_nll,zero_count_conservative_upper_B,nominal_upper_tail,sideband_deviance,flags,window_yield_coefficients\n";
  for(size_t j=0;j<choices.size();++j){std::vector<Result>results;
    for(size_t r=0;r<c.regions.size();++r)for(size_t d=0;d<c.degrees.size();++d){Result fit=Background(hist[j],c.regions[r],c.degrees[d],c);results.push_back(fit);
      csv<<CsvQuote(choices[j].id)<<",full,"<<c.regions[r].name<<','<<fit.degree<<','<<int(fit.valid)<<','<<fit.iterations<<','<<fit.kkt<<','<<fit.side<<','<<fit.nwin<<',';
      if(fit.valid)csv<<fit.b<<','<<fit.purity;else csv<<',';csv<<','<<int(fit.boundary)<<','<<int(fit.profileValid)<<',';
      if(fit.profileValid)csv<<fit.upper<<','<<1-fit.upper/fit.nwin;else csv<<',';csv<<','<<c.delta<<',';if(fit.valid&&fit.side==0)csv<<fit.zeroUpper;
      csv<<','<<.5*TMath::Erfc(std::sqrt(c.delta))<<','<<fit.deviance<<','<<CsvQuote(fit.flags.empty()?"ok":fit.flags)<<',';std::ostringstream coefficients;coefficients<<std::setprecision(17);for(size_t k=0;k<fit.coefficients.size();++k){if(k)coefficients<<';';coefficients<<fit.coefficients[k];}csv<<CsvQuote(coefficients.str())<<'\n';
    }Diagnostic(choices[j],hist[j],c,results,Form("%s_background_%02u",stem.c_str(),(unsigned)j),j,out);
  }csv.close();Require(bool(csv),"background CSV write failure");
}
} // namespace LambdaPurityRefinement

Int_t checkLambdaPurityRefinement(const char*scanDirectory,const char*outputPrefix){
  using namespace LambdaPurityRefinement;
  try{
    Require(scanDirectory&&*scanDirectory&&outputPrefix&&*outputPrefix,"missing arguments");std::string dir(scanDirectory),stem(outputPrefix);Config c=Manifest(dir+"/plot_manifest.tsv");
    std::vector<Row>choiceRows=ReadCsv(dir+"/chosen_selections.csv");std::vector<Choice>choices;std::map<std::string,int>id;
    for(size_t j=0;j<choiceRows.size();++j){Choice ch;ch.id=Get(choiceRows[j],"choice_id");ch.label=Get(choiceRows[j],"label");ch.role=Get(choiceRows[j],"role");Require(id.insert(std::make_pair(ch.id,int(j))).second,"duplicate choice id");choices.push_back(ch);}Require(!choices.empty()&&choices[0].role=="baseline","baseline must be first");
    const int nbins=int(std::floor((c.hi-c.lo)/c.bw+.5));Require(nbins>0&&std::fabs(c.lo+nbins*c.bw-c.hi)<1.e-9,"invalid binning");std::vector<TH1D*>hist;std::vector<std::vector<bool> >seen(choices.size(),std::vector<bool>(nbins,false));Vector norms(choices.size(),-1.);
    for(size_t j=0;j<choices.size();++j){TH1D*h=new TH1D(Form("RefinementFullCounts_%u",(unsigned)j),"",nbins,c.lo,c.hi);h->SetDirectory(0);h->Sumw2();hist.push_back(h);}
    std::vector<Row>rows=ReadCsv(dir+"/histograms.csv");for(size_t i=0;i<rows.size();++i){const Row&r=rows[i];if(Get(r,"split")!="full")continue;Require(id.count(Get(r,"choice_id")),"unknown histogram choice");int j=id[Get(r,"choice_id")],b=int(N(r,"bin"));Require(b==N(r,"bin")&&b>=1&&b<=nbins&&!seen[j][b-1],"invalid/duplicate bin");seen[j][b-1]=true;
      Require(std::fabs(N(r,"low")-hist[j]->GetBinLowEdge(b))<1.e-9&&std::fabs(N(r,"high")-hist[j]->GetBinLowEdge(b+1))<1.e-9,"different bin edges");double count=N(r,"count"),error=N(r,"error");Require(count>=0&&std::fabs(count-std::floor(count+.5))<1.e-8&&error>=0&&std::fabs(error*error-count)<1.e-7*std::max(1.,count),"invalid unweighted count/error");hist[j]->SetBinContent(b,count);hist[j]->SetBinError(b,error);double norm=N(r,"normalization_count");if(norms[j]<0)norms[j]=norm;else Require(norms[j]==norm,"normalization metadata differs");}
    for(size_t j=0;j<choices.size();++j){for(int b=0;b<nbins;++b)Require(seen[j][b],"missing bin");Require(hist[j]->Integral(1,nbins)==norms[j],"normalization count mismatch");Edge(hist[j],c.wl);Edge(hist[j],c.wh);for(size_t r=0;r<c.regions.size();++r){Edge(hist[j],c.regions[r].ll);Edge(hist[j],c.regions[r].lh);Edge(hist[j],c.regions[r].rl);Edge(hist[j],c.regions[r].rh);}}
    rows=ReadCsv(dir+"/chosen_metrics.csv");std::map<std::string,Row>metrics;
    for(size_t i=0;i<rows.size();++i)if(Get(rows[i],"split")=="full"&&Get(rows[i],"sideband")=="default"){std::string key=Get(rows[i],"choice_id");Require(id.count(key)&&metrics.insert(std::make_pair(key,rows[i])).second,"duplicate/unknown metric");TH1D*h=hist[id[key]];Require(h->Integral(Edge(h,c.wl),Edge(h,c.wh)-1)==N(rows[i],"Nwin"),"metric window mismatch");}
    Require(metrics.size()==choices.size(),"missing full default metrics");
    const int references=std::min(2,int(choices.size()));const int overlays=std::max(1,int((choices.size()-references+1)/2));std::vector<std::string>products;products.push_back(stem+".root");products.push_back(stem+"_background.csv");products.push_back(stem+"_audit.txt");const char*ext[]={".pdf",".png"};
    for(int e=0;e<2;++e){for(int k=0;k<overlays;++k)products.push_back(Form("%s_overlay_%02d%s",stem.c_str(),k,ext[e]));for(size_t j=0;j<choices.size();++j)products.push_back(Form("%s_background_%02u%s",stem.c_str(),(unsigned)j,ext[e]));}
    for(size_t k=0;k<products.size();++k)Require(gSystem->AccessPathName(products[k].c_str()),"refusing overwrite "+products[k]);
    gROOT->SetBatch(true);gStyle->SetOptStat(0);gStyle->SetTextFont(42);TFile out((stem+".root").c_str(),"CREATE");Require(!out.IsZombie(),"cannot create ROOT output");for(size_t j=0;j<hist.size();++j)hist[j]->Write();
    std::vector<int>extras;const char*preferred[]={"robust70","nominal70"};
    for(int p=0;p<2;++p)if(id.count(preferred[p])&&id[preferred[p]]>=references)extras.push_back(id[preferred[p]]);
    for(size_t j=references;j<choices.size();++j)if(std::find(extras.begin(),extras.end(),int(j))==extras.end())extras.push_back(j);
    for(int k=0;k<overlays;++k){std::vector<int>which;for(int j=0;j<references;++j)which.push_back(j);for(size_t j=2*k;j<extras.size()&&j<size_t(2*k+2);++j)which.push_back(extras[j]);Overlay(choices,which,hist,metrics,c,k,Form("%s_overlay_%02d",stem.c_str(),k),out);}
    Checks(choices,hist,c,stem,out);
    std::ostringstream audit;audit<<"Full-data exploratory refinement; saved CSV only; no fresh held-out validation.\nsource="<<dir<<"\n"<<c.raw;
    audit<<"All bins checked; mass window fixed; normalized by complete original histogram integral before display zoom. No Gaussian signal model, new mass constraint, smoothing or rebinning.\n"
      <<"Background: sideband-only integrated-bin Poisson likelihood; nonnegative Bernstein coefficients with domains from outer sideband edges. Degree0/1/2 are configured sensitivity models, not alternative cut-selection scores.\n"
      <<"Convex mixture optimum checked by KKT residual. Boundary coefficients allowed; low sideband count is a warning, not rejection. Covariance/Hessian errors are deliberately not used at boundaries.\n"
      <<"Profile B=sum(window-yield coefficients), with shape reoptimized. Upper DeltaNLL threshold is NOMINAL ASYMPTOTIC and does not guarantee sparse/boundary coverage. No uncertainty from cut selection or background-family choice is included.\n"
      <<"For zero sideband counts, shape is unidentified: Bhat=0 is not BG-free evidence. Separate conservative upper B=-log(nominalTail)/min(SB/window basis ratio), valid only within the chosen positive-basis family and without selection correction.\n"
      <<"Nonnegative Bernstein coefficients are a restricted background family. Fits do not exclude signal leakage, peaking reflections, cut sculpting, or feeddown. No fresh validation or detector efficiency is claimed.\n";
    out.cd();TObjString(audit.str().c_str()).Write("RefinementBackgroundAudit");out.Flush();Require(!out.TestBit(TFile::kWriteError),"ROOT write failure");out.Close();std::ofstream log((stem+"_audit.txt").c_str());log<<audit.str();log.close();Require(bool(log),"audit write failure");std::cout<<"[checkLambdaPurityRefinement] PASS: "<<choices.size()<<" full-data choices; sideband models are cross-checks only.\n";return 0;
  }catch(const std::exception&e){std::cerr<<"[checkLambdaPurityRefinement] ERROR: "<<e.what()<<std::endl;return 1;}
}
#endif // ROOT5 CINT sees public entry points only

