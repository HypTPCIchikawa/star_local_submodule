// Study only: observed t + 3He + 4He peaks jointly anchor the nonnegative yields.
// Every named isotope mean/width is FIXED to its independent original PID polynomial.
void FitPeakAnchored(const char*,const char*,const char*,const char*,const char*);
#ifndef __CINT__
#include "../FitLocalNuclearCutoff.C"
#include <iostream>

namespace NuclearPeakAnchored {
using namespace NuclearLocalCutoff;

void Numbers(std::ostream& out,const Vec& values) {
  out<<'[';for(unsigned j=0;j<values.size();++j){if(j)out<<',';Scalar(out,values[j]);}out<<']';
}
void Strings(std::ostream& out,const std::vector<std::string>& values) {
  out<<'[';for(unsigned j=0;j<values.size();++j){if(j)out<<',';out<<Quote(values[j]);}out<<']';
}
void Covariance(std::ostream& out,const std::vector<Vec>& values) {
  out<<'[';for(unsigned j=0;j<values.size();++j){if(j)out<<',';Numbers(out,values[j]);}out<<']';
}
void AddFlag(std::vector<std::string>& flags,const std::string& value) {
  if(std::find(flags.begin(),flags.end(),value)==flags.end())flags.push_back(value);
}
template<class E> double Integral(const typename E::Component& c,double lo,double hi) {
  // Deliberately independent of c.active: omitted d must retain its nonzero bound.
  if(c.nodes.empty())return E::Integral(lo,hi,c.mean,c.sigma);
  double value=0;for(unsigned q=0;q<c.nodes.size();++q)
    value+=c.nodes[q].weight*E::Integral(lo,hi,c.nodes[q].mean,c.nodes[q].sigma);
  return value;
}
struct Region {
  double lo,hi,observed,predicted,pearson,maxPull;int first,last,bins;bool available;
  Region():lo(0),hi(0),observed(0),predicted(0),pearson(0),maxPull(0),first(0),last(-1),bins(0),available(false){}
};
Region Bounds(TH1D& raw,double lo,double hi,double alignment,double tol) {
  Require(alignment>0&&Finite(lo)&&Finite(hi)&&hi>lo,"Invalid explicit region/alignment");
  Region r;r.lo=std::floor((lo+tol)/alignment)*alignment;
  r.hi=std::ceil((hi-tol)/alignment)*alignment;
  r.first=raw.FindBin(r.lo+tol);r.last=raw.FindBin(r.hi-tol);r.bins=r.last-r.first+1;
  Require(r.first>=1&&r.last<=raw.GetNbinsX()&&r.last>=r.first&&
          std::fabs(raw.GetBinLowEdge(r.first)-r.lo)<tol&&
          std::fabs(raw.GetXaxis()->GetBinUpEdge(r.last)-r.hi)<tol,
          "Diagnostic region must contain exact whole native bins within the axis");
  return r;
}
Region PeakRegion(TH1D& raw,const Slice& slice,int left,int right,double leftSigma,
                  double rightSigma,double alignment,double tol) {
  double lo=std::numeric_limits<double>::infinity(),hi=-lo;
  for(int q=0;q<3;++q) {
    const Row& a=slice.rows[4*q+left];const Row& b=slice.rows[4*q+right];
    lo=std::min(lo,a.mean-leftSigma*a.sigma);hi=std::max(hi,b.mean+rightSigma*b.sigma);
  }
  return Bounds(raw,lo,hi,alignment,tol);
}
template<class E> void Diagnose(Region& region,const TH1D& raw,const typename E::Result& fit,int first) {
  region.available=!fit.totalCounts.empty();
  for(int b=region.first;b<=region.last;++b) {
    region.observed+=raw.GetBinContent(b);const int local=b-first;
    if(local<0||local>=int(fit.totalCounts.size())){region.available=false;continue;}
    const double mu=fit.totalCounts[local],n=raw.GetBinContent(b);
    if(!Finite(mu)||mu<0){region.available=false;continue;}
    region.predicted+=mu;
    if(mu>0){const double pull=(n-mu)/std::sqrt(mu);region.pearson+=pull*pull;region.maxPull=std::max(region.maxPull,std::fabs(pull));}
    else if(n>0)region.available=false;
  }
}
void WriteRegion(std::ostream& out,const Region& r) {
  out<<"{\"lo\":"<<r.lo<<",\"hi\":"<<r.hi<<",\"firstNativeBin\":"<<r.first
     <<",\"lastNativeBin\":"<<r.last<<",\"binCount\":"<<r.bins<<",\"observedCount\":"<<r.observed
     <<",\"predictionAvailable\":"<<Boolean(r.available)<<",\"predictedCount\":";
  if(r.available)Scalar(out,r.predicted);else out<<"null";
  out<<",\"pearson\":";if(r.available)Scalar(out,r.pearson);else out<<"null";
  out<<",\"maxAbsPearsonPull\":";if(r.available)Scalar(out,r.maxPull);else out<<"null";
  out<<",\"ndf\":null,\"ndfPolicy\":\"Regional diagnostic only; fitted-parameter subtraction belongs to the joint likelihood, not each overlapping region\"}";
}
template<class E> Estimate RawGate(const typename E::Request& req,const Slice& slice,const TH1D& raw,
                                 const typename E::Result& fit,int first,int target,double ns,
                                 double tol,std::vector<Regions>& regions) {
  Estimate out;regions.clear();
  for(int q=0;q<3;++q) {
    std::vector<Choice> choices;for(int s=0;s<4;++s)choices.push_back(slice.rows[4*q+s].gate);
    regions.push_back(SelectionRegions(choices,target,ns,req.nativeAxisMin,req.nativeAxisMax));
    const Row& row=slice.rows[4*q+target];out.identitySupport=out.identitySupport&&row.supported;
    out.gateActive=out.gateActive||row.gate.active;
    if(row.gate.active)out.axisCovered=out.axisCovered&&row.gate.mean-ns*row.gate.sigma>=req.nativeAxisMin&&
      row.gate.mean+ns*row.gate.sigma<=req.nativeAxisMax;
    for(unsigned k=0;k<regions[q].size();++k)out.covered=out.covered&&
      regions[q][k].lo>=slice.fitLo-tol&&regions[q][k].hi<=slice.fitHi+tol;
  }
  const bool prediction=fit.totalCounts.size()==req.bins.size();
  for(int b=1;b<=raw.GetNbinsX();++b) {
    double fraction=0,full=0,touch=0;
    for(int q=0;q<3;++q) {
      const double overlap=Overlap(regions[q],raw.GetBinLowEdge(b),raw.GetXaxis()->GetBinUpEdge(b));
      const double w=slice.rows[4*q].weight;fraction+=w*overlap;
      if(overlap>=1-tol)full+=w;if(overlap>tol)touch+=w;
    }
    const double n=raw.GetBinContent(b);out.rawCount+=fraction*n;out.rawLower+=full*n;out.rawUpper+=touch*n;
    if(prediction&&b>=first&&b<first+int(fit.totalCounts.size()))out.absoluteResidual+=fraction*std::fabs(n-fit.totalCounts[b-first]);
  }
  out.localResidual=out.rawCount>0?out.absoluteResidual/out.rawCount:0;return out;
}
void WriteRawSelection(std::ostream& out,const Estimate& e,const std::vector<Regions>& regions,
                       bool numerical,bool prediction,double omitted) {
  out<<"{\"available\":false,\"purity\":null,\"signalSelected\":null,\"totalSelected\":null,\"backgroundSelected\":null,\"conditionalError\":null,\"reason\":\"peak_anchored_fixed_response_quality_pilot_not_purity_certification\",\"rawGateCount\":";
  Scalar(out,e.rawCount);out<<",\"rawGateCountLowerFullBins\":";Scalar(out,e.rawLower);
  out<<",\"rawGateCountUpperTouchingBins\":";Scalar(out,e.rawUpper);
  out<<",\"gateAbsoluteResidualCount\":";if(prediction)Scalar(out,e.absoluteResidual);else out<<"null";
  out<<",\"localResidualFrac\":";if(prediction&&e.rawCount>0)Scalar(out,e.localResidual);else out<<"null";
  out<<",\"omittedSelectedDDescriptiveUpperBound\":";Scalar(out,omitted);
  out<<",\"numericalValid\":"<<Boolean(numerical)<<",\"identitySupport\":"<<Boolean(e.identitySupport)
     <<",\"identitySupportMeaning\":\"Original function reference, not independent calibration certification\",\"selectionCovered\":"<<Boolean(e.covered)
     <<",\"axisCovered\":"<<Boolean(e.axisCovered)<<",\"gateActive\":"<<Boolean(e.gateActive)<<",\"gateIntervalsAtNodes\":";
  WriteRegions(out,regions);out<<'}';
}
template<class E> void SaveDiagnostics(TFile& file,const TH1D& raw,const typename E::Result& fit,
                                     int first,int index,const Region& control,const Region& helium) {
  file.cd();const char* names[]={"pearson_pull","t_control_mask","helium_mask"};
  for(int part=0;part<3;++part) {
    TH1D* h=static_cast<TH1D*>(raw.Clone(Form("diagnostic_%s_%03d",names[part],index)));
    h->SetDirectory(0);h->Reset();
    for(int b=1;b<=raw.GetNbinsX();++b) {
      if(part==1)h->SetBinContent(b,b>=control.first&&b<=control.last?1:0);
      else if(part==2)h->SetBinContent(b,b>=helium.first&&b<=helium.last?1:0);
      else {const int local=b-first;if(local>=0&&local<int(fit.totalCounts.size())&&fit.totalCounts[local]>0)
        h->SetBinContent(b,(raw.GetBinContent(b)-fit.totalCounts[local])/std::sqrt(fit.totalCounts[local]));}
    }
    Require(h->Write()>0,"Cannot save anchored-region diagnostic histogram");delete h;
  }
}
void ReadDiagnosticBack(const char* path,const std::vector<Slice>& slices,double tol) {
  TFile file(path,"READ");Require(!file.IsZombie()&&!file.TestBit(TFile::kRecovered),"Cannot reopen anchored diagnostics");
  for(unsigned i=0;i<slices.size();++i) {
    const int index=slices[i].index;
    TH1* raw=dynamic_cast<TH1*>(file.Get(Form("raw_native_%03d",index)));
    TH1* prediction=dynamic_cast<TH1*>(file.Get(Form("prediction_total_%03d",index)));
    TH1* fitMask=dynamic_cast<TH1*>(file.Get(Form("prediction_fit_mask_%03d",index)));
    TH1* pull=dynamic_cast<TH1*>(file.Get(Form("diagnostic_pearson_pull_%03d",index)));
    TH1* control=dynamic_cast<TH1*>(file.Get(Form("diagnostic_t_control_mask_%03d",index)));
    TH1* helium=dynamic_cast<TH1*>(file.Get(Form("diagnostic_helium_mask_%03d",index)));
    Require(raw&&prediction&&fitMask&&pull&&control&&helium,"Missing anchored-region histogram");
    for(int b=0;b<=raw->GetNbinsX()+1;++b) {
      Require(Finite(pull->GetBinContent(b)),"Nonfinite saved Pearson pull");
      Require((control->GetBinContent(b)==0||control->GetBinContent(b)==1)&&
              (helium->GetBinContent(b)==0||helium->GetBinContent(b)==1),"Invalid region mask");
      Require((!control->GetBinContent(b)&&!helium->GetBinContent(b))||fitMask->GetBinContent(b)==1,
              "Saved control/helium diagnostic lies outside joint fit");
      if(prediction->GetBinContent(b)>0) {
        const double value=(raw->GetBinContent(b)-prediction->GetBinContent(b))/std::sqrt(prediction->GetBinContent(b));
        Require(std::fabs(value-pull->GetBinContent(b))<=tol*std::max(1.,std::fabs(value)),"Closed ROOT Pearson pull differs from saved data/model");
      }
    }
  }
  file.Close();std::cout<<"PEAK_ANCHORED_DIAGNOSTIC_READBACK_PASS rows="<<slices.size()<<std::endl;
}

template<class E> void Run(const char* cachePath,const char* output,const Settings& cfg,
                          const std::vector<Slice>& slices) {
  Require(Text(cfg,"sampleHistogram")=="ALL"&&Text(cfg,"targetSpecies")=="he4"&&
          Text(cfg,"fitSpecies")=="t,he3,he4"&&Integer(cfg,"fixedShape")==1,
          "Require explicit ALL/t+he3+he4 fixed-shape anchored pilot");
  Require(Text(cfg,"shapePolicy")=="fixed_original_pid_functions"&&Text(cfg,"targetShapeMode")=="fixed",
          "Every isotope must retain its original independent polynomial mean/width");
  Require(Text(cfg,"responseStudyMode")=="peak_anchored_fixed_pid","Unknown explicit anchored study mode");
  Settings::const_iterator free=cfg.find("freeSpecies");
  Require(free!=cfg.end()&&free->second.empty(),"No free isotope center/width is allowed in anchored pilot");
  Require(Text(cfg,"axisCountBoundMethod")=="observed_plus_nsigma_sqrt","Unknown explicit descriptive count-bound method");
  const double tol=Number(cfg,"axisTolerance"),ns=Number(cfg,"maxNSigma"),floor=E::AcceptanceFloor(cfg);
  const double alignment=Number(cfg,"fitWindowAlignmentKeV"),left=Number(cfg,"roiLeftSigma"),right=Number(cfg,"roiRightSigma");
  const double heLeft=Number(cfg,"heRegionLeftSigma"),heRight=Number(cfg,"heRegionRightSigma"),tSigma=Number(cfg,"tControlSigma");
  const double material=Number(cfg,"materialOmissionFraction"),countSigma=Number(cfg,"axisCountBoundNSigma");
  Require(tol>0&&ns>0&&floor>0&&alignment>0&&left>0&&right>0&&heLeft>0&&heRight>0&&tSigma>0&&material>0&&countSigma>=0,"Invalid explicit pilot geometry/diagnostic settings");
  TString folder(output);Require(!gSystem->AccessPathName(folder),"Output folder must be provenance-prepared");
  Require(gSystem->AccessPathName(folder+"/cutoff_fits.json")&&gSystem->AccessPathName(folder+"/cutoff_fits.root"),"Existing anchored outputs will not be replaced");
  TFile cache(cachePath,"READ");Require(!cache.IsZombie()&&!cache.TestBit(TFile::kRecovered),"Cannot open intact ALL cache");
  TH2* source=dynamic_cast<TH2*>(cache.Get("ALL"));Require(source!=0,"Missing ALL histogram");
  TFile product(folder+"/cutoff_fits.root","CREATE");Require(!product.IsZombie(),"Cannot CREATE anchored ROOT");
  std::ofstream json((folder+"/cutoff_fits.json").Data());Require(bool(json),"Cannot CREATE anchored JSON");json<<std::setprecision(17);
  json<<"{\"schemaVersion\":1,\"responseModel\":"<<Quote(Text(cfg,"responseModel"))
      <<",\"responseStudyMode\":"<<Quote(Text(cfg,"responseStudyMode"))<<",\"backgroundModel\":"<<Quote(Text(cfg,"backgroundModel"))
      <<",\"targetSpecies\":\"he4\",\"sampleHistogram\":\"ALL\",\"shapeSource\":\"legacy_log10_polynomial\",\"shapePolicy\":\"fixed_original_pid_functions\",\"fixedShape\":true,\"backgroundPeakShapesFree\":false,\"fixedBackgroundNuclearMeans\":true,\"fixedBackgroundNuclearWidths\":true,\"freeSpecies\":[],\"fitSpecies\":[\"t\",\"he3\",\"he4\"],\"productionEnabled\":false,\"purityAvailable\":false,\"fullUncertaintyValidated\":false,\"calibrationUncertaintyPropagated\":false,\"normalizationPolicy\":\"One contiguous joint ALL-data likelihood contains the t peak and the 3He left flank through 4He right flank. Fixed independent nuclear functions, nonnegative window yields. Yfull=Yfit/Afit; Yaxis=Yfull*Aaxis. TOF-tag absolute yields are not transferred to ALL.\",\"omissionBoundPolicy\":\"Descriptive worst-case omitted-d yield from ALL in-axis observed count plus Nsigma sqrt(N), divided by d axis acceptance. Not an exact confidence bound, not a yield constraint or purity certification.\",\"axisCountBoundNSigma\":"<<countSigma
      <<",\"materialOmissionFraction\":"<<material<<",\"regionalDiagnosticsNotIndependent\":true,\"regionalNdfPolicy\":\"Only the joint likelihood has conditionalNdf. Overlapping regional Pearson contributions do not define independent reduced chi-squared values.\",\"slices\":[";
  std::set<int> used;
  for(unsigned i=0;i<slices.size();++i) {
    const Slice& slice=slices[i];Require(used.insert(slice.rBin).second&&slice.rBin<=source->GetNbinsX()&&
      std::fabs(source->GetXaxis()->GetBinLowEdge(slice.rBin)-slice.lo)<tol&&
      std::fabs(source->GetXaxis()->GetBinUpEdge(slice.rBin)-slice.hi)<tol,"Anchored R bin differs from cache");
    TH1D* raw=source->ProjectionY(Form("raw_native_%03d",slice.index),slice.rBin,slice.rBin,"e");raw->SetDirectory(0);
    Region joint=Bounds(*raw,slice.fitLo,slice.fitHi,alignment,tol);
    Region control=PeakRegion(*raw,slice,1,1,tSigma,tSigma,alignment,tol);
    Region helium=PeakRegion(*raw,slice,2,3,heLeft,heRight,alignment,tol);
    Region required=PeakRegion(*raw,slice,1,3,left,right,alignment,tol);
    for(int q=0;q<3;++q)required.lo=std::min(required.lo,std::floor((slice.rows[4*q+2].mean-left*slice.rows[4*q+2].sigma+tol)/alignment)*alignment);
    Require(std::fabs(joint.lo-slice.fitLo)<tol&&std::fabs(joint.hi-slice.fitHi)<tol&&
      joint.lo<=required.lo+tol&&joint.hi>=required.hi-tol&&joint.first<=control.first&&joint.last>=control.last&&
      joint.first<=helium.first&&joint.last>=helium.last,"Joint window must include observed t control and full 3He/4He diagnostic regions");
    const int first=joint.first,last=joint.last,rebin=Integer(cfg,"displayDedxRebinFactor");
    Require(rebin>0&&raw->GetNbinsX()%rebin==0,"Invalid display-only rebin factor");
    typename E::Request req;req.nativeAxisMin=raw->GetXaxis()->GetXmin();req.nativeAxisMax=raw->GetXaxis()->GetXmax();
    req.backgroundModel=Text(cfg,"backgroundModel");req.backgroundSlope=Number(cfg,"backgroundStartSlope");req.freeBackgroundSlope=Flag(Integer(cfg,"freeBackgroundSlope"));
    for(int b=0;b<=raw->GetNbinsX()+1;++b) {
      const double n=raw->GetBinContent(b),variance=std::pow(raw->GetBinError(b),2);
      Require(Finite(n)&&n>=0&&Finite(variance)&&std::fabs(variance-n)<=Number(cfg,"countVarianceRelativeTolerance")*std::max(1.,n),"Native ALL cache must retain Poisson count/error/flow");
      if(b>=1&&b<=raw->GetNbinsX())Require(std::fabs(raw->GetBinWidth(b)-Number(cfg,"nativeDedxBinWidth"))<tol,"Native dEdx width changed");
      if(b>=first&&b<=last)req.bins.push_back(typename E::Bin(raw->GetBinLowEdge(b),raw->GetXaxis()->GetBinUpEdge(b),n));
    }
    Vec acceptance(4,0),axisAcceptance(4,0);
    for(int s=0;s<4;++s) {
      typename E::Component c;c.name=Species[s];c.mean=slice.rows[4+s].mean;c.sigma=slice.rows[4+s].sigma;
      c.freeMean=false;c.freeSigma=false;
      for(int q=0;q<3;++q) {const Row& row=slice.rows[4*q+s];c.nodes.push_back(typename E::Node(row.weight,row.mean,row.sigma));}
      acceptance[s]=Integral<E>(c,joint.lo,joint.hi);axisAcceptance[s]=Integral<E>(c,req.nativeAxisMin,req.nativeAxisMax);
      c.active=s!=0; // d is explicitly omitted, never activated via a tiny normalized tail.
      if(c.active)Require(Finite(acceptance[s])&&acceptance[s]>floor,"Named observed peak has insufficient joint-window acceptance");
      req.components.push_back(c);
    }
    double count=0;for(unsigned b=0;b<req.bins.size();++b)count+=req.bins[b].count;
    const double per=count/(3+(req.backgroundModel=="linear"?2:1));
    for(int s=0;s<4;++s)req.components[s].startWindowYield=s?per:0;
    req.backgroundWindowYield=per*(req.backgroundModel=="linear"?2:1);
    typename E::Result fit;std::string failure;
    try {fit=E::Fit(req,cfg);}catch(const std::exception& error){failure=error.what();}
    const bool prediction=fit.totalCounts.size()==req.bins.size()&&fit.components.size()==4;
    const bool numerical=failure.empty()&&prediction&&fit.finite&&fit.amplitudeConverged&&fit.pointAccepted&&fit.status==0&&!fit.atBoundary;
    if(prediction)for(int s=0;s<4;++s) {
      Require(fit.components[s].mean==req.components[s].mean&&fit.components[s].sigma==req.components[s].sigma&&
        !fit.components[s].freeMean&&!fit.components[s].freeSigma,"Fixed original isotope functions changed in fit");
      for(int q=0;q<3;++q)Require(
        std::fabs(fit.components[s].nodes[q].mean-req.components[s].nodes[q].mean)<=tol*std::max(1.,std::fabs(req.components[s].nodes[q].mean))&&
        std::fabs(fit.components[s].nodes[q].sigma-req.components[s].nodes[q].sigma)<=tol*std::max(1.,std::fabs(req.components[s].nodes[q].sigma)),
        "Fixed original R-node function changed beyond numerical arithmetic tolerance");
    }
    Diagnose<E>(joint,*raw,fit,first);Diagnose<E>(control,*raw,fit,first);Diagnose<E>(helium,*raw,fit,first);
    SaveSpectra<E>(product,*raw,req,fit,first,slice.index,rebin);SaveDiagnostics<E>(product,*raw,fit,first,slice.index,control,helium);
    unsigned activeYields=0;for(unsigned j=0;j<fit.yieldCovarianceColumnActive.size();++j)activeYields+=fit.yieldCovarianceColumnActive[j];
    if(prediction)Require(fit.conditionalNdf==int(req.bins.size())-int(fit.shapeNames.size())-int(activeYields),"Joint NDF must retain true free background/yield parameter counts");
    const double observed=raw->Integral(1,raw->GetNbinsX()),axisBound=observed+countSigma*std::sqrt(observed);
    const bool dRatio=Finite(axisAcceptance[0])&&axisAcceptance[0]>0;
    const double dFitBound=dRatio?axisBound*acceptance[0]/axisAcceptance[0]:axisBound;
    const double dHeAcceptance=Integral<E>(req.components[0],helium.lo,helium.hi);
    const double dHeBound=dRatio?axisBound*dHeAcceptance/axisAcceptance[0]:axisBound;
    std::vector<std::vector<Regions> > gateRegions(4);std::vector<Estimate> gates;
    for(int s=0;s<4;++s)gates.push_back(RawGate<E>(req,slice,*raw,fit,first,s,ns,tol,gateRegions[s]));
    double dSelectedAcceptance=0;
    for(int q=0;q<3;++q)for(unsigned k=0;k<gateRegions[3][q].size();++k) {
      const Interval& interval=gateRegions[3][q][k];const Row& d=slice.rows[4*q];
      dSelectedAcceptance+=d.weight*E::Integral(interval.lo,interval.hi,d.mean,d.sigma);
    }
    const double dGateBound=dRatio?axisBound*dSelectedAcceptance/axisAcceptance[0]:axisBound;
    const bool dMaterial=!dRatio||!(helium.observed>0)||!(gates[3].rawCount>0)||
      dHeBound/helium.observed>material||dGateBound/gates[3].rawCount>material;
    std::vector<std::string> flags;
    if(!numerical)AddFlag(flags,"joint_fit_numerically_unavailable");
    if(!fit.yieldCovarianceValid)AddFlag(flags,"conditional_yield_covariance_invalid");
    if(dMaterial)AddFlag(flags,"omitted_d_bound_material_or_unavailable");
    if(!gates[3].covered||!gates[3].axisCovered)AddFlag(flags,"he4_gate_not_covered_by_joint_window_or_axis");
    bool tAxisBoundAvailable=false;double tFull=0,tAxis=0;
    if(prediction&&fit.windowAcceptance.size()>1&&fit.windowAcceptance[1]>0&&fit.columnWindowYields.size()>1) {
      tFull=fit.columnWindowYields[1]/fit.windowAcceptance[1];tAxis=tFull*axisAcceptance[1];tAxisBoundAvailable=Finite(tAxis);
      if(tAxis>axisBound)AddFlag(flags,"fitted_t_axis_yield_exceeds_descriptive_ALL_count_margin");
    }
    if(i)json<<',';
    json<<"{\"index\":"<<slice.index<<",\"rBin\":"<<slice.rBin<<",\"lo\":"<<slice.lo<<",\"hi\":"<<slice.hi
        <<",\"r\":"<<.5*(slice.lo+slice.hi)<<",\"fitLo\":"<<joint.lo<<",\"fitHi\":"<<joint.hi
        <<",\"fitCount\":"<<count<<",\"shapeSource\":\"legacy_log10_polynomial\",\"allNuclearShapeParametersFixed\":true,\"fixedShapeAlignmentIndependentlyValidated\":false,\"purityAvailable\":false,\"fitFailure\":"<<Quote(failure)
        <<",\"predictionAvailable\":"<<Boolean(prediction)<<",\"numericalValid\":"<<Boolean(numerical)<<",\"status\":"<<fit.status
        <<",\"covarianceStatus\":"<<fit.covarianceStatus<<",\"pointAccepted\":"<<Boolean(fit.pointAccepted)<<",\"atBoundary\":"<<Boolean(fit.atBoundary)
        <<",\"amplitudeConverged\":"<<Boolean(fit.amplitudeConverged)<<",\"amplitudeKktResidual\":";Scalar(json,fit.amplitudeKkt);
    json<<",\"edm\":";Scalar(json,fit.edm);json<<",\"pearson\":";if(prediction)Scalar(json,fit.pearson);else json<<"null";
    json<<",\"dataDeviance\":";if(prediction)Scalar(json,fit.dataDeviance);else json<<"null";
    json<<",\"conditionalNdf\":"<<fit.conditionalNdf<<",\"freeShapeParameterCount\":"<<fit.shapeNames.size()
        <<",\"interiorYieldParameterCount\":"<<activeYields<<",\"shapeNames\":";Strings(json,fit.shapeNames);
    json<<",\"columnNames\":";Strings(json,fit.columnNames);json<<",\"columnWindowYields\":";Numbers(json,fit.columnWindowYields);
    json<<",\"yieldCovarianceValid\":"<<Boolean(fit.yieldCovarianceValid)<<",\"fixedShapeYieldCovariance\":";Covariance(json,fit.fixedShapeYieldCovariance);
    json<<",\"observedNativeInAxisCount\":"<<observed<<",\"descriptiveAxisCountUpperMargin\":"<<axisBound
        <<",\"omittedD\":{\"active\":false,\"axisAcceptance\":"<<axisAcceptance[0]<<",\"fitAcceptance\":"<<acceptance[0]
        <<",\"heliumAcceptance\":"<<dHeAcceptance<<",\"selectedHe4Acceptance\":"<<dSelectedAcceptance
        <<",\"ratioAvailable\":"<<Boolean(dRatio)<<",\"fitYieldDescriptiveUpperBound\":"<<dFitBound
        <<",\"heliumYieldDescriptiveUpperBound\":"<<dHeBound<<",\"selectedHe4YieldDescriptiveUpperBound\":"<<dGateBound
        <<",\"materialOrUnassessed\":"<<Boolean(dMaterial)<<",\"exactConfidenceBound\":false}"
        <<",\"tNormalization\":{\"sameALLObservedPeakInLikelihood\":true,\"fullLineYield\":";
    if(tAxisBoundAvailable)Scalar(json,tFull);else json<<"null";
    json<<",\"axisYield\":";if(tAxisBoundAvailable)Scalar(json,tAxis);else json<<"null";
    json<<",\"withinDescriptiveALLCountMargin\":";if(tAxisBoundAvailable)json<<Boolean(tAxis<=axisBound);else json<<"null";
    json<<"},\"postFitFlags\":";Strings(json,flags);
    json<<",\"regions\":{\"joint\":";WriteRegion(json,joint);json<<",\"tControl\":";WriteRegion(json,control);
    json<<",\"helium\":";WriteRegion(json,helium);json<<"},\"components\":{";
    for(int s=0;s<4;++s) {
      const typename E::Component& c=prediction?fit.components[s]:req.components[s];
      const double a=prediction&&fit.windowAcceptance.size()>unsigned(s)?fit.windowAcceptance[s]:acceptance[s];
      const bool yAvailable=prediction&&fit.columnWindowYields.size()>unsigned(s);
      const double y=yAvailable?fit.columnWindowYields[s]:0;
      if(s)json<<',';json<<Quote(Species[s])<<":{\"mean\":"<<c.mean<<",\"sigma\":"<<c.sigma
          <<",\"initialMean\":"<<req.components[s].mean<<",\"initialSigma\":"<<req.components[s].sigma
          <<",\"meanError\":null,\"sigmaError\":null,\"fixedFromFunction\":true,\"shapeSource\":\"legacy_log10_polynomial\",\"calibrationUncertaintyPropagated\":false,\"active\":"<<Boolean(c.active)
          <<",\"freeMean\":false,\"freeSigma\":false,\"windowAcceptance\":"<<acceptance[s]<<",\"modelReturnedWindowAcceptance\":"<<a<<",\"axisAcceptance\":"<<axisAcceptance[s]
          <<",\"windowYield\":";if(yAvailable)Scalar(json,y);else json<<"null";
      json<<",\"fullLineYield\":";if(yAvailable&&c.active&&a>0)Scalar(json,y/a);else json<<"null";
      json<<",\"axisYield\":";if(yAvailable&&c.active&&a>0)Scalar(json,y*axisAcceptance[s]/a);else json<<"null";
      json<<",\"rNodes\":[";
      for(unsigned q=0;q<c.nodes.size();++q){if(q)json<<',';json<<"{\"weight\":"<<c.nodes[q].weight<<",\"mean\":"<<c.nodes[q].mean<<",\"sigma\":"<<c.nodes[q].sigma<<'}';}
      json<<"],\"selection\":";WriteRawSelection(json,gates[s],gateRegions[s],numerical,prediction,s==3?dGateBound:0);json<<'}';
    }
    json<<"}}";delete raw;
    std::cout<<"PEAK_ANCHORED_SLICE index="<<slice.index<<" numerical="<<numerical<<" chi2="<<fit.pearson<<" ndf="<<fit.conditionalNdf<<" dMaterial="<<dMaterial<<std::endl;
  }
  json<<"],\"executionIsNotPhysicalValidation\":true}\n";json.close();Require(bool(json),"Cannot finish anchored JSON");product.Close();
  Readback((folder+"/cutoff_fits.root").Data(),*source,slices,Number(cfg,"countVarianceRelativeTolerance"));
  ReadDiagnosticBack((folder+"/cutoff_fits.root").Data(),slices,Number(cfg,"countVarianceRelativeTolerance"));cache.Close();
  std::cout<<"PEAK_ANCHORED_COMPLETE rows="<<slices.size()<<std::endl;
}
}
#endif
void FitPeakAnchored(const char* cachePath,const char* output,const char* settings,const char* shapes,const char* gates) {
#ifndef __CINT__
  try {
    const NuclearFit::Settings cfg=NuclearFit::ReadSettings(settings);
    const std::vector<NuclearLocalCutoff::Slice> slices=NuclearLocalCutoff::ReadInputs(shapes,gates,NuclearFit::Number(cfg,"axisTolerance"));
    const std::string response=NuclearFit::Text(cfg,"responseModel");
    if(response=="gaussian")NuclearPeakAnchored::Run<NuclearLocalCutoff::GaussianEngine>(cachePath,output,cfg,slices);
    else {NuclearFit::Require(response=="lognormal","Unknown explicit fixed response model");
      NuclearPeakAnchored::Run<NuclearLocalCutoff::LognormalEngine>(cachePath,output,cfg,slices);}
  }catch(const std::exception& error){std::cerr<<"FitPeakAnchored: "<<error.what()<<std::endl;gSystem->Exit(2);}
#endif
}
