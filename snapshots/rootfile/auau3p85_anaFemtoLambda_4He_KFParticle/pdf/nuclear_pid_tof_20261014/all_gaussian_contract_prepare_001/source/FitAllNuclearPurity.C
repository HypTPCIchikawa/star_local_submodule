// Study-only fixed-shape fits and exact TPC-selection purity integrals.
// TAG spectra are never added to ALL. No production calibration is written.
void FitAllNuclearPurity(const char*,const char*,const char*,const char*,const char*);
#ifndef __CINT__
#include "NuclearFitCommon.h"
#include "TFile.h"
#include "TH2.h"
#include <set>
namespace NuclearAll {
using namespace NuclearFit;
const char* const Species[4]={"d","t","he3","he4"};
struct Choice {double mean,sigma;bool active;Choice():mean(0),sigma(0),active(false){} };
struct ShapeRow {
  int index,rBin,node,species;double lo,hi,weight,r,mean,sigma,meanError,sigmaError;
  bool supported,candidateActive;Choice current;
};
struct Slice {int index,rBin;double lo,hi;std::vector<ShapeRow> rows;};
struct Window {int index,rBin;double lo,hi;};
struct Interval {double lo,hi;Interval(double a=0,double b=0):lo(a),hi(b){} };
typedef std::vector<Interval> Regions;
struct Estimate {
  Vec selectedFractions,selectedCounts;double signal,total,purity,error;
  bool covered,axisCovered,available,errorAvailable,boundaryYieldPresent,unsupportedNuisancePresent,droppedNuisancePresent;
  Estimate():signal(0),total(0),purity(0),error(0),covered(true),axisCovered(true),available(false),errorAvailable(false),boundaryYieldPresent(false),unsupportedNuisancePresent(false),droppedNuisancePresent(false){}
};
inline int SpeciesIndex(const std::string& name){for(int s=0;s<4;++s)if(name==Species[s])return s;return -1;}
inline void Scalar(std::ostream& out,double value){if(Finite(value))out<<value;else out<<"null";}
inline const char* Boolean(bool value){return value?"true":"false";}
inline bool Flag(int value){Require(value==0||value==1,"Boolean TSV value must be 0 or 1");return value!=0;}
std::vector<Slice> ReadShapes(const char* path,double tolerance){
  std::ifstream input(path);Require(bool(input),"Missing explicit shapes TSV");std::map<int,Slice> byIndex;std::set<std::string> seen;std::string line;
  const std::string header="index\trBin\trLo\trHi\tnode\tweight\trValue\tspecies\tmean\tsigma\tmean_error\tsigma_error\tsupported\tcurrent_mean\tcurrent_sigma\tcurrent_active\tcandidate_active";
  while(std::getline(input,line)){
    if(line==header)continue;std::istringstream in(line);ShapeRow row;std::string species,extra;int supported,currentActive,candidateActive;
    Require(bool(in>>row.index>>row.rBin>>row.lo>>row.hi>>row.node>>row.weight>>row.r>>species>>row.mean>>row.sigma>>row.meanError>>row.sigmaError>>supported>>row.current.mean>>row.current.sigma>>currentActive>>candidateActive)&&!(in>>extra),"Malformed 17-field shapes row");
    row.species=SpeciesIndex(species);row.supported=Flag(supported);row.current.active=Flag(currentActive);row.candidateActive=Flag(candidateActive);
    Require(row.index>=0&&row.rBin>=1&&row.node>=0&&row.node<3&&row.species>=0&&Finite(row.lo)&&Finite(row.hi)&&row.hi>row.lo&&row.lo>0&&Finite(row.r)&&row.r>row.lo&&row.r<row.hi&&Finite(row.weight)&&row.weight>0&&Finite(row.mean)&&Finite(row.sigma)&&row.sigma>0&&Finite(row.meanError)&&Finite(row.sigmaError),"Invalid shape geometry/kernel");
    Require(!row.current.active||(Finite(row.current.mean)&&Finite(row.current.sigma)&&row.current.sigma>0),"Active current selector requires finite positive sigma");
    std::ostringstream key;key<<row.index<<':'<<row.node<<':'<<row.species;Require(seen.insert(key.str()).second,"Duplicate shape index/node/species");
    Slice& slice=byIndex[row.index];if(slice.rows.empty()){slice.index=row.index;slice.rBin=row.rBin;slice.lo=row.lo;slice.hi=row.hi;}
    Require(slice.rBin==row.rBin&&std::fabs(slice.lo-row.lo)<tolerance&&std::fabs(slice.hi-row.hi)<tolerance,"Inconsistent slice geometry");slice.rows.push_back(row);
  }
  Require(input.eof()&&!byIndex.empty(),"Unreadable or empty shapes TSV");std::vector<Slice> result;
  for(std::map<int,Slice>::iterator i=byIndex.begin();i!=byIndex.end();++i){Slice& slice=i->second;Require(slice.rows.size()==12,"Exactly 3 GL nodes x 4 species are required");
    std::sort(slice.rows.begin(),slice.rows.end(),[](const ShapeRow& a,const ShapeRow& b){return a.node*4+a.species<b.node*4+b.species;});double weight=0;
    for(int q=0;q<3;++q){const ShapeRow& reference=slice.rows[q*4];Require(reference.node==q&&reference.species==0,"Missing canonical kernel row");weight+=reference.weight;
      for(int s=0;s<4;++s){const ShapeRow& row=slice.rows[q*4+s];Require(row.node==q&&row.species==s&&std::fabs(row.weight-reference.weight)<tolerance&&std::fabs(row.r-reference.r)<tolerance,"All species must use identical GL nodes/weights");}}
    Require(std::fabs(weight-1)<tolerance,"GL weights do not sum to one");result.push_back(slice);
  }return result;
}
std::map<int,Window> ReadWindows(const char* path){
  std::ifstream input(path);Require(bool(input),"Missing fit_windows TSV");std::map<int,Window> result;std::string line;
  while(std::getline(input,line)){if(line=="index\trBin\tfitLo\tfitHi")continue;std::istringstream in(line);Window row;std::string extra;
    Require(bool(in>>row.index>>row.rBin>>row.lo>>row.hi)&&!(in>>extra),"Malformed 4-field fit window row");Require(row.index>=0&&row.rBin>=1&&Finite(row.lo)&&Finite(row.hi)&&row.hi>row.lo,"Invalid fit window");Require(result.insert(std::make_pair(row.index,row)).second,"Duplicate fit window");}
  Require(input.eof()&&!result.empty(),"Unreadable/empty fit window table");return result;
}
inline int Select(const std::vector<Choice>& choices,double x,double maxNSigma){
  int best=-1;double distance=maxNSigma;
  for(unsigned s=0;s<choices.size();++s)if(choices[s].active){double pull=std::fabs((x-choices[s].mean)/choices[s].sigma);
    // Strict PID gate and original d<t<he3<he4 tie order.
    if(pull<distance){distance=pull;best=s;}}
  return best;
}
inline void AddEdge(Vec& edges,double edge,double lo,double hi){if(Finite(edge)&&edge>lo&&edge<hi)edges.push_back(edge);}
Regions SelectionRegions(const std::vector<Choice>& choices,int target,double maxNSigma,double lo,double hi){
  Require(choices.size()==4&&target>=0&&target<4&&Finite(maxNSigma)&&maxNSigma>0&&hi>lo,"Invalid selector integration request");Vec edges;edges.push_back(lo);edges.push_back(hi);
  for(unsigned a=0;a<choices.size();++a)if(choices[a].active){const Choice& x=choices[a];Require(Finite(x.mean)&&Finite(x.sigma)&&x.sigma>0,"Invalid active selector shape");
    AddEdge(edges,x.mean-maxNSigma*x.sigma,lo,hi);AddEdge(edges,x.mean+maxNSigma*x.sigma,lo,hi);
    for(unsigned b=0;b<a;++b)if(choices[b].active){const Choice& y=choices[b];
      AddEdge(edges,(y.sigma*x.mean+x.sigma*y.mean)/(y.sigma+x.sigma),lo,hi);
      if(y.sigma!=x.sigma)AddEdge(edges,(y.sigma*x.mean-x.sigma*y.mean)/(y.sigma-x.sigma),lo,hi);}}
  std::sort(edges.begin(),edges.end());edges.erase(std::unique(edges.begin(),edges.end()),edges.end());Regions result;
  for(unsigned i=1;i<edges.size();++i)if(edges[i]>edges[i-1]){double mid=edges[i-1]+.5*(edges[i]-edges[i-1]);
    if(Select(choices,mid,maxNSigma)==target){if(!result.empty()&&result.back().hi==edges[i-1])result.back().hi=edges[i];else result.push_back(Interval(edges[i-1],edges[i]));}}
  return result;
}
Estimate IntegrateSelection(const Result& fitted,const Request& request,const std::vector<Regions>& regions,const Vec& weights,int target,double minimumDenominator,double tolerance){
  Estimate out;const unsigned n=fitted.components.size(),columns=fitted.columnWindowYields.size();Require(n==4&&regions.size()==3&&weights.size()==3&&columns>=n,"Incomplete GL purity inputs");
  const double fitLo=request.bins.front().lo,fitHi=request.bins.back().hi;out.selectedFractions.assign(columns,0);out.selectedCounts.assign(columns,0);
  for(unsigned q=0;q<regions.size();++q)for(unsigned k=0;k<regions[q].size();++k){const Interval& region=regions[q][k];
    if(region.lo<fitLo-tolerance||region.hi>fitHi+tolerance)out.covered=false;double lo=std::max(fitLo,region.lo),hi=std::min(fitHi,region.hi);if(!(hi>lo))continue;
    for(unsigned s=0;s<n;++s)if(fitted.components[s].active&&fitted.windowAcceptance[s]>0){Require(fitted.components[s].nodes.size()==3,"Missing final GL Gaussian nodes");const Node& node=fitted.components[s].nodes[q];out.selectedFractions[s]+=weights[q]*GaussianIntegral(lo,hi,node.mean,node.sigma)/fitted.windowAcceptance[s];}
    for(unsigned j=n;j<columns;++j)out.selectedFractions[j]+=weights[q]*BackgroundIntegral(lo,hi,fitLo,fitHi,request.backgroundModel,fitted.backgroundSlope,j-n);
  }
  for(unsigned j=0;j<columns;++j){Require(Finite(out.selectedFractions[j])&&out.selectedFractions[j]>=0&&out.selectedFractions[j]<=1+tolerance,"Invalid selected component fraction");out.selectedCounts[j]=out.selectedFractions[j]*fitted.columnWindowYields[j];out.total+=out.selectedCounts[j];
    if(out.selectedFractions[j]>0&&fitted.components.size()<=j)out.boundaryYieldPresent=out.boundaryYieldPresent||!fitted.yieldCovarianceColumnActive[j];
    if(out.selectedFractions[j]>0&&j<n&&fitted.components[j].active)out.boundaryYieldPresent=out.boundaryYieldPresent||!fitted.yieldCovarianceColumnActive[j];}
  out.signal=out.selectedCounts[target];out.available=Finite(out.total)&&Finite(out.signal)&&out.total>minimumDenominator;
  if(out.available){out.purity=out.signal/out.total;Require(out.purity>=0&&out.purity<=1+tolerance,"Purity outside nonnegative model range");}
  if(out.available&&fitted.yieldCovarianceValid&&fitted.yieldCovarianceColumnActive[target]&&fitted.fixedShapeYieldCovariance.size()==columns){
    Vec gradient(columns,0);for(unsigned j=0;j<columns;++j)gradient[j]=((j==unsigned(target)?out.selectedFractions[target]*out.total:0)-out.signal*out.selectedFractions[j])/(out.total*out.total);
    double variance=0;for(unsigned j=0;j<columns;++j){Require(fitted.fixedShapeYieldCovariance[j].size()==columns,"Wrong yield covariance dimension");for(unsigned k=0;k<columns;++k)variance+=gradient[j]*fitted.fixedShapeYieldCovariance[j][k]*gradient[k];}
    out.errorAvailable=Finite(variance)&&variance>=0;if(out.errorAvailable)out.error=std::sqrt(variance);
  }return out;
}
std::vector<Choice> Choices(const Slice& slice,int q,bool current){std::vector<Choice> result(4);
  for(int s=0;s<4;++s){const ShapeRow& row=slice.rows[q*4+s];if(current)result[s]=row.current;else{result[s].mean=row.mean;result[s].sigma=row.sigma;result[s].active=row.candidateActive;}}return result;}
bool ShapeSupported(const Slice& slice,int target){for(int q=0;q<3;++q)if(!slice.rows[q*4+target].supported)return false;return true;}
bool SelectorAxisCovered(const Slice& slice,int target,bool current,double maxNSigma,double lo,double hi){
  bool active=false;for(int q=0;q<3;++q){std::vector<Choice> choice=Choices(slice,q,current);if(!choice[target].active)continue;active=true;
    if(choice[target].mean-maxNSigma*choice[target].sigma<lo||choice[target].mean+maxNSigma*choice[target].sigma>hi)return false;}return active;
}
void WriteEstimate(std::ostream& out,const Estimate& estimate,bool eligible,bool supported,bool numerical,const std::string& policy){
  out<<"{\"selectionPolicy\":"<<Quote(policy)<<",\"signalSelected\":"<<estimate.signal<<",\"totalSelected\":"<<estimate.total<<",\"backgroundSelected\":"<<estimate.total-estimate.signal<<",\"purity\":";if(estimate.available)Scalar(out,estimate.purity);else out<<"null";
  out<<",\"conditionalFixedShapeError\":";if(estimate.errorAvailable)Scalar(out,estimate.error);else out<<"null";
  out<<",\"selectedColumnFractions\":";WriteVec(out,estimate.selectedFractions);out<<",\"selectedColumnCounts\":";WriteVec(out,estimate.selectedCounts);
  out<<",\"available\":"<<Boolean(estimate.available)<<",\"selectionCoveredByFitWindow\":"<<Boolean(estimate.covered)<<",\"selectionWindowInsideDetectorAxis\":"<<Boolean(estimate.axisCovered)<<",\"boundaryYieldPresent\":"<<Boolean(estimate.boundaryYieldPresent)<<",\"unsupportedContributingNuisancePresent\":"<<Boolean(estimate.unsupportedNuisancePresent)<<",\"droppedNuisancePresent\":"<<Boolean(estimate.droppedNuisancePresent)<<",\"shapeSupported\":"<<Boolean(supported)<<",\"numericalValid\":"<<Boolean(numerical)<<",\"conditionalEligible\":"<<Boolean(eligible)<<",\"referenceOnly\":"<<Boolean(!eligible)<<",\"productionEnabled\":false,\"fullUncertaintyValidated\":false}";
}
void SavePredictions(TFile& output,TH1D& native,int first,const Request& request,const Result& result,int index,int rebin){
  output.cd();Require(native.Write()>0,"Cannot write native ALL projection");TH1D* display=static_cast<TH1D*>(native.Clone(Form("raw_display_%03d",index)));display->SetDirectory(0);double before=native.Integral(0,native.GetNbinsX()+1);display->Rebin(rebin);Require(display->Integral(0,display->GetNbinsX()+1)==before,"Display rebin changed original count");Require(display->Write()>0,"Cannot write display projection");delete display;
  const unsigned n=request.components.size();std::vector<TH1D*> hist;
  for(unsigned j=0;j<n+3;++j){std::string name=j==0?"total":(j==1?"smooth":(j==2?"fit_mask":request.components[j-3].name));TH1D* h=static_cast<TH1D*>(native.Clone(Form("prediction_%s_%03d",name.c_str(),index)));h->SetDirectory(0);h->Reset();hist.push_back(h);}
  for(unsigned b=0;b<request.bins.size();++b){hist[2]->SetBinContent(first+b,1);if(b<result.totalCounts.size()){hist[0]->SetBinContent(first+b,result.totalCounts[b]);hist[1]->SetBinContent(first+b,result.smoothCounts[b]);for(unsigned j=0;j<n;++j)hist[j+3]->SetBinContent(first+b,result.componentCounts[b][j]);}}
  for(unsigned j=0;j<hist.size();++j){Require(hist[j]->Write()>0,"Cannot write model prediction");delete hist[j];}
}
void WriteFit(std::ostream& json,const Result& fitted,const Request& request,const std::string& failure){
  json<<"\"fitFailure\":"<<Quote(failure)<<",\"dataDeviance\":";Scalar(json,fitted.dataDeviance);json<<",\"chi2\":";Scalar(json,fitted.pearson);
  json<<",\"conditionalNdf\":"<<fitted.conditionalNdf<<",\"status\":"<<fitted.status<<",\"edm\":";Scalar(json,fitted.edm);
  json<<",\"pointAccepted\":"<<Boolean(fitted.pointAccepted)<<",\"finite\":"<<Boolean(fitted.finite)<<",\"amplitudeConverged\":"<<Boolean(fitted.amplitudeConverged)<<",\"amplitudeKktResidual\":";Scalar(json,fitted.amplitudeKkt);
  json<<",\"amplitudeIterations\":"<<fitted.amplitudeIterations<<",\"atBoundary\":"<<Boolean(fitted.atBoundary)<<",\"backgroundModel\":"<<Quote(request.backgroundModel)<<",\"backgroundSlope\":"<<fitted.backgroundSlope<<",\"columnNames\":[";
  for(unsigned j=0;j<fitted.columnNames.size();++j){if(j)json<<',';json<<Quote(fitted.columnNames[j]);}json<<"],\"columnWindowYields\":";WriteVec(json,fitted.columnWindowYields);json<<",\"windowAcceptance\":";WriteVec(json,fitted.windowAcceptance);
  json<<",\"fixedShapeYieldCovarianceValid\":"<<Boolean(fitted.yieldCovarianceValid)<<",\"fixedShapeYieldCovariance\":";WriteMatrix(json,fitted.fixedShapeYieldCovariance);
  json<<",\"yieldCovarianceColumnActive\":[";for(unsigned j=0;j<fitted.yieldCovarianceColumnActive.size();++j){if(j)json<<',';json<<Boolean(fitted.yieldCovarianceColumnActive[j]);}json<<']';
  json<<",\"uncertaintyScope\":\"Interior nonnegative yields only, retained Gaussian functions and BG slope fixed; excludes TAG uncertainty, TAG-ALL correlation, transfer bias and background-model choice\",\"fullUncertaintyValidated\":false";
}
}
#endif
void FitAllNuclearPurity(const char* cacheRoot,const char* output,const char* settings,const char* shapesTsv,const char* windowsTsv){
#ifndef __CINT__
  using namespace NuclearAll;
  try{
    Settings cfg=ReadSettings(settings);ValidateSettings(cfg);double tolerance=Number(cfg,"axisTolerance"),maxNSigma=Number(cfg,"maxNSigma"),minimumDenominator=Number(cfg,"minimumPurityDenominator");
    Require(tolerance>0&&maxNSigma>0&&minimumDenominator>=0&&Number(cfg,"maximumConditionalChi2Ndf")>0,"Invalid explicit purity criteria");
    const std::string sample=Text(cfg,"sampleHistogram");Require(sample=="ALL"||sample=="REMAINING"||sample=="TAG_D"||sample=="TAG_T"||sample=="TAG_HE3","Unknown sample; spectra must never be added");
    std::vector<Slice> slices=ReadShapes(shapesTsv,tolerance);std::map<int,Window> windows=ReadWindows(windowsTsv);Require(slices.size()==windows.size(),"Shapes and fit window coverage differ");
    const TString folder(output);const char* products[]={"all_fits.root","all_fits.json","purity.tsv","completed.json"};
    for(unsigned j=0;j<sizeof(products)/sizeof(products[0]);++j)Require(gSystem->AccessPathName(folder+"/"+products[j]),"Existing ALL purity products will not be replaced");
    Require(!gSystem->AccessPathName(folder+"/context.json")||!gSystem->AccessPathName(folder+"/../context.json"),"ALL output must be provenance-prepared");
    TFile cache(cacheRoot,"READ");Require(!cache.IsZombie(),"Cannot read immutable cache");TH2* source=dynamic_cast<TH2*>(cache.Get(sample.c_str()));Require(source!=0,"Missing selected cache histogram");
    TFile product(folder+"/all_fits.root","CREATE");Require(!product.IsZombie(),"Cannot CREATE new ALL product");std::ofstream json((folder+"/all_fits.json").Data()),table((folder+"/purity.tsv").Data());Require(bool(json)&&bool(table),"Cannot create new ALL JSON/TSV");json<<std::setprecision(17);table<<std::setprecision(17);
    json<<"{\"schemaVersion\":1,\"fitStage\":\"conditional_fixed_GL3_nuclear_shapes\",\"sampleHistogram\":"<<Quote(sample)<<",\"productionEnabled\":false,\"fullUncertaintyValidated\":false,\"maxNSigma\":"<<maxNSigma<<",\"requireBestSpecies\":true,\"tieOrder\":[\"d\",\"t\",\"he3\",\"he4\"],\"RBinPolicy\":\"Three uniform-R GL nodes, same weights for each component and selector; no sub-native statistical independence\",\"likelihoodPolicy\":\"Only the explicitly selected sample; TAG never added to ALL\",\"purityPolicy\":\"S/(sum all modeled Gaussian contributions + positive smooth BG), exact interval partition and bin-window normalized integrals; no observed coverage extrapolation\",\"slices\":[";
    table<<"index\trBin\trLo\trHi\tspecies\tselection\tsignal\ttotal\tpurity\tconditionalError\tcovered\tshapeSupported\tnumericalValid\tconditionalEligible\n";
    std::set<int> used;bool allNumerical=true;double totalD=0;unsigned eligibleCount=0;
    for(unsigned i=0;i<slices.size();++i){const Slice& slice=slices[i];std::map<int,Window>::const_iterator wi=windows.find(slice.index);Require(wi!=windows.end()&&wi->second.rBin==slice.rBin,"Missing/mismatched fit window");const Window& window=wi->second;
      Require(used.insert(slice.rBin).second&&slice.rBin<=source->GetNbinsX()&&std::fabs(source->GetXaxis()->GetBinLowEdge(slice.rBin)-slice.lo)<tolerance&&std::fabs(source->GetXaxis()->GetBinUpEdge(slice.rBin)-slice.hi)<tolerance&&std::fabs(slice.hi-slice.lo-Number(cfg,"nativeRigidityBinWidth"))<tolerance,"R slice is not one unique native bin");
      TH1D* native=source->ProjectionY(Form("raw_native_%03d",slice.index),slice.rBin,slice.rBin,"e");native->SetDirectory(0);int first=native->FindBin(window.lo+tolerance),last=native->FindBin(window.hi-tolerance),rebin=Integer(cfg,"displayDedxRebinFactor");
      Require(rebin>0&&native->GetNbinsX()%rebin==0&&first>=1&&last>=first&&last<=native->GetNbinsX()&&std::fabs(native->GetXaxis()->GetBinLowEdge(first)-window.lo)<tolerance&&std::fabs(native->GetXaxis()->GetBinUpEdge(last)-window.hi)<tolerance&&(first-1)%rebin==0&&last%rebin==0,"Fit window is not exact native/draw-rebin bounds");
      Request request;request.nativeAxisMin=native->GetXaxis()->GetXmin();request.nativeAxisMax=native->GetXaxis()->GetXmax();request.backgroundModel=Text(cfg,"backgroundModel");request.backgroundSlope=Number(cfg,"backgroundStartSlope");request.freeBackgroundSlope=Flag(Integer(cfg,"freeBackgroundSlope"));
      for(int b=0;b<=native->GetNbinsX()+1;++b){double count=native->GetBinContent(b),variance=std::pow(native->GetBinError(b),2);Require(Finite(count)&&count>=0&&Finite(variance)&&variance>=0&&std::fabs(variance-count)<=Number(cfg,"countVarianceRelativeTolerance")*std::max(1.,count),"Cached native spectrum is not Poisson-type nonnegative counts");
        if(b>=1&&b<=native->GetNbinsX())Require(std::fabs(native->GetBinWidth(b)-Number(cfg,"nativeDedxBinWidth"))<tolerance,"Unexpected native dEdx width");if(b>=first&&b<=last)request.bins.push_back(Bin(native->GetBinLowEdge(b),native->GetXaxis()->GetBinUpEdge(b),count));}
      std::vector<std::string> dropped(4);int active=0;Vec weights(3);for(int q=0;q<3;++q)weights[q]=slice.rows[q*4].weight;
      for(int s=0;s<4;++s){Component c;c.name=Species[s];c.mean=slice.rows[4+s].mean;c.sigma=slice.rows[4+s].sigma;c.freeMean=false;c.freeSigma=false;
        for(int q=0;q<3;++q){const ShapeRow& row=slice.rows[q*4+s];c.nodes.push_back(Node(row.weight,row.mean,row.sigma));
          if(row.mean<request.nativeAxisMin-Number(cfg,"meanAxisPadding")||row.mean>request.nativeAxisMax+Number(cfg,"meanAxisPadding")||row.sigma<Number(cfg,"minimumSigma")||row.sigma>Number(cfg,"maximumSigma"))dropped[s]="outside_declared_mathematical_support";}
        if(dropped[s].empty()&&!(KernelIntegral(c,window.lo,window.hi)>Number(cfg,"minimumGaussianAcceptance")))dropped[s]="fit_window_gaussian_acceptance_below_floor";
        c.active=dropped[s].empty();if(c.active)++active;request.components.push_back(c);}
      const double count=TotalCount(request.bins);const int backgroundColumns=request.backgroundModel=="linear"?2:1;double perColumn=count/std::max(1,active+backgroundColumns);
      for(int s=0;s<4;++s)request.components[s].startWindowYield=request.components[s].active?perColumn:0;request.backgroundWindowYield=perColumn*backgroundColumns;
      Result fitted;std::string failure;try{fitted=Fit(request,cfg);}catch(const std::exception& e){failure=e.what();}
      bool numerical=failure.empty()&&fitted.finite&&fitted.amplitudeConverged&&fitted.pointAccepted&&fitted.status==0;allNumerical=allNumerical&&numerical;if(fitted.finite)totalD+=fitted.dataDeviance;
      bool quality=numerical&&fitted.conditionalNdf>0&&fitted.pearson/fitted.conditionalNdf<=Number(cfg,"maximumConditionalChi2Ndf")&&!fitted.atBoundary;
      SavePredictions(product,*native,first,request,fitted,slice.index,rebin);
      if(i)json<<',';json<<"{\"index\":"<<slice.index<<",\"rBin\":"<<slice.rBin<<",\"lo\":"<<slice.lo<<",\"hi\":"<<slice.hi<<",\"r\":"<<(slice.lo+slice.hi)/2<<",\"fitLo\":"<<window.lo<<",\"fitHi\":"<<window.hi<<",\"firstFitDedxBin\":"<<first<<",\"lastFitDedxBin\":"<<last<<",\"fitCount\":"<<count<<",\"observedNativeCountWithFlow\":"<<native->Integral(0,native->GetNbinsX()+1)<<",\"numericalValid\":"<<Boolean(numerical)<<",\"conditionalFitQuality\":"<<Boolean(quality)<<',';WriteFit(json,fitted,request,failure);
      json<<",\"droppedComponents\":{";for(int s=0;s<4;++s){if(s)json<<',';json<<Quote(Species[s])<<':'<<Quote(dropped[s]);}json<<"},\"components\":{";
      for(int s=0;s<4;++s){if(s)json<<',';const ShapeRow& center=slice.rows[4+s];json<<Quote(Species[s])<<":{\"mean\":"<<center.mean<<",\"sigma\":"<<center.sigma<<",\"meanError\":"<<center.meanError<<",\"sigmaError\":"<<center.sigmaError<<",\"shapeConstraint\":\"fixed conditional GL3 function\",\"shapeSupported\":"<<Boolean(ShapeSupported(slice,s))<<",\"identityValidated\":false,\"windowYield\":";
        if(fitted.columnWindowYields.size()>unsigned(s))Scalar(json,fitted.columnWindowYields[s]);else json<<"null";json<<",\"windowAcceptance\":";if(fitted.windowAcceptance.size()>unsigned(s))Scalar(json,fitted.windowAcceptance[s]);else json<<"null";
        json<<",\"effectiveCount\":";if(fitted.effectiveCount.size()>unsigned(s))Scalar(json,fitted.effectiveCount[s]);else json<<"null";json<<",\"conditionalDropDeltaD\":";if(fitted.conditionalDropDeltaD.size()>unsigned(s))Scalar(json,fitted.conditionalDropDeltaD[s]);else json<<"null";json<<",\"nodes\":[";
        for(int q=0;q<3;++q){if(q)json<<',';const ShapeRow& row=slice.rows[q*4+s];json<<"{\"r\":"<<row.r<<",\"weight\":"<<row.weight<<",\"mean\":"<<row.mean<<",\"sigma\":"<<row.sigma<<",\"meanError\":"<<row.meanError<<",\"sigmaError\":"<<row.sigmaError<<",\"supported\":"<<Boolean(row.supported)<<",\"candidateActive\":"<<Boolean(row.candidateActive)<<",\"currentActive\":"<<Boolean(row.current.active)<<",\"currentMean\":";Scalar(json,row.current.mean);json<<",\"currentSigma\":";Scalar(json,row.current.sigma);json<<'}';}
        json<<"],\"selections\":{";
        for(int method=0;method<3;++method){const char* name=method==0?"current_tpc_best_species":(method==1?"candidate_tpc_best_species":"single_gaussian_2sigma_reference");if(method)json<<',';json<<Quote(name)<<':';
          Estimate estimate;bool supported=ShapeSupported(slice,s),eligible=false;
          if(fitted.components.size()==4&&fitted.columnWindowYields.size()>=4){std::vector<Regions> regions(3);
            for(int q=0;q<3;++q){if(method<2)regions[q]=SelectionRegions(Choices(slice,q,method==0),s,maxNSigma,request.nativeAxisMin,request.nativeAxisMax);
              else {const ShapeRow& row=slice.rows[q*4+s];double lo=std::max(request.nativeAxisMin,row.mean-maxNSigma*row.sigma),hi=std::min(request.nativeAxisMax,row.mean+maxNSigma*row.sigma);if(hi>lo)regions[q].push_back(Interval(lo,hi));}}
            estimate=IntegrateSelection(fitted,request,regions,weights,s,minimumDenominator,tolerance);estimate.axisCovered=method<2?SelectorAxisCovered(slice,s,method==0,maxNSigma,request.nativeAxisMin,request.nativeAxisMax):true;
            if(method==2)for(int q=0;q<3;++q){const ShapeRow& row=slice.rows[q*4+s];estimate.axisCovered=estimate.axisCovered&&row.mean-maxNSigma*row.sigma>=request.nativeAxisMin&&row.mean+maxNSigma*row.sigma<=request.nativeAxisMax;}
            for(int other=0;other<4;++other)if(other!=s){estimate.unsupportedNuisancePresent=estimate.unsupportedNuisancePresent||(!ShapeSupported(slice,other)&&estimate.selectedCounts[other]>0);estimate.droppedNuisancePresent=estimate.droppedNuisancePresent||!dropped[other].empty();}
            eligible=method<2&&sample=="ALL"&&quality&&supported&&dropped[s].empty()&&estimate.available&&estimate.covered&&estimate.axisCovered&&fitted.columnWindowYields[s]>0&&!estimate.boundaryYieldPresent;
          }
          WriteEstimate(json,estimate,eligible,supported,numerical,name);eligibleCount+=eligible;
          table<<slice.index<<'\t'<<slice.rBin<<'\t'<<slice.lo<<'\t'<<slice.hi<<'\t'<<Species[s]<<'\t'<<name<<'\t'<<estimate.signal<<'\t'<<estimate.total<<'\t';if(estimate.available)table<<estimate.purity;else table<<"nan";table<<'\t';if(estimate.errorAvailable)table<<estimate.error;else table<<"nan";table<<'\t'<<estimate.covered<<'\t'<<supported<<'\t'<<numerical<<'\t'<<eligible<<'\n';
        }json<<"}}";
      }json<<"}}";std::cout<<"NUCLEAR_ALL_SLICE index="<<slice.index<<" sample="<<sample<<" R="<<(slice.lo+slice.hi)/2<<" D="<<fitted.dataDeviance<<" numerical="<<numerical<<std::endl;delete native;
    }
    json<<"],\"totalFiniteDataDeviance\":"<<totalD<<",\"numericalConverged\":"<<Boolean(allNumerical)<<",\"conditionalEligibleEstimates\":"<<eligibleCount<<",\"executionIsNotPhysicalValidation\":true}\n";json.close();table.close();product.Close();cache.Close();std::cout<<"NUCLEAR_ALL_PURITY_COMPLETE rows="<<slices.size()<<std::endl;
  }catch(const std::exception& e){std::cerr<<"FitAllNuclearPurity: "<<e.what()<<std::endl;gSystem->Exit(2);}
#endif
}
