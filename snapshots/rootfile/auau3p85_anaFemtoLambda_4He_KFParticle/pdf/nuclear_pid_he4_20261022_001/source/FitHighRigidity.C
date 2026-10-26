// Study-only TOF-conditional adapter. The low-R implementation stays unchanged.
void FitHighRigidity(const char*,const char*,const char*,const char*,const char*);
#ifndef __CINT__
#include "FitLocalNuclearCutoff.C"
namespace NuclearHighRigidity {
using namespace NuclearLocalCutoff;
inline bool AllowedSample(const std::string& sample) {
  return sample=="TAG_D" || sample=="TAG_T" || sample=="TAG_HE3";
}
struct RawGate {double count,lower,upper;RawGate():count(0),lower(0),upper(0){} };
RawGate Measure(const Slice& slice,const TH1D& raw,int target,double ns,double tol) {
  std::vector<Regions> regions(3);
  for(int q=0;q<3;++q) {
    std::vector<Choice> choices(4);
    for(int s=0;s<4;++s) choices[s]=slice.rows[4*q+s].gate;
    regions[q]=SelectionRegions(choices,target,ns,raw.GetXaxis()->GetXmin(),raw.GetXaxis()->GetXmax());
  }
  RawGate value;
  for(int b=1;b<=raw.GetNbinsX();++b) for(int q=0;q<3;++q) {
    double fraction=Overlap(regions[q],raw.GetBinLowEdge(b),raw.GetXaxis()->GetBinUpEdge(b));
    double weighted=raw.GetBinContent(b)*slice.rows[4*q].weight;
    value.count+=weighted*fraction;
    if(fraction>=1-tol) value.lower+=weighted;
    if(fraction>tol) value.upper+=weighted;
  }
  Require(Finite(value.count)&&value.count>=0&&value.lower<=value.count+tol*std::max(1.,value.count)&&value.upper+tol*std::max(1.,value.count)>=value.count,"Invalid fit-independent TAG gate count");
  return value;
}
template<class E> void Run(const char* cachePath,const char* output,const Settings& cfg,const std::vector<Slice>& slices) {
  const std::string sample=Text(cfg,"sampleHistogram"),targetName=Text(cfg,"targetSpecies");
  const int target=SpeciesIndex(targetName);
  Require(AllowedSample(sample)&&target>=0,"Explicit supported TOF sample/target required");
  Require(TargetMode(cfg)=="free_target_mean_sigma","High-R adapter explicitly requires free target response");
  Require(Number(cfg,"mass2Max")>Number(cfg,"mass2Min")&&Number(cfg,"rigidityMin")>=Number(cfg,"tofThreshold"),"Invalid explicit TOF window/domain");
  const double tol=Number(cfg,"axisTolerance"),ns=Number(cfg,"maxNSigma"),minimum=Number(cfg,"minimumPurityDenominator"),floor=E::AcceptanceFloor(cfg);
  Require(tol>0&&ns>0&&minimum>=0&&floor>0,"Invalid explicit integration criteria");
  TString folder(output);
  Require(!gSystem->AccessPathName(folder),"Output directory must be provenance-prepared");
  const char* products[]={"cutoff_fits.json","cutoff_fits.root","purity.tsv","raw_gate_counts.tsv"};
  for(int j=0;j<4;++j) Require(gSystem->AccessPathName(folder+"/"+products[j]),"Existing high-R product will not be replaced");
  TFile cache(cachePath,"READ");
  Require(!cache.IsZombie()&&!cache.TestBit(TFile::kRecovered),"Cannot read intact native TAG cache");
  TH2* source=dynamic_cast<TH2*>(cache.Get(sample.c_str()));
  Require(source!=0,"Missing explicitly selected TAG histogram");
  TFile product(folder+"/cutoff_fits.root","CREATE");
  Require(!product.IsZombie(),"Cannot CREATE fitted spectra");
  std::ofstream json((folder+"/cutoff_fits.json").Data()),table((folder+"/purity.tsv").Data()),counts((folder+"/raw_gate_counts.tsv").Data());
  Require(bool(json)&&bool(table)&&bool(counts),"Cannot create high-R tables");
  json<<std::setprecision(17);table<<std::setprecision(17);counts<<std::setprecision(17);
  json<<"{\"schemaVersion\":1,\"fitStage\":\"high_r_tof_conditional_target_free\",\"responseModel\":"<<Quote(Text(cfg,"responseModel"))
      <<",\"backgroundModel\":"<<Quote(Text(cfg,"backgroundModel"))<<",\"targetSpecies\":"<<Quote(targetName)
      <<",\"targetShapeMode\":\"free_target_mean_sigma\",\"sampleHistogram\":"<<Quote(sample)
      <<",\"mass2Window\":["<<Number(cfg,"mass2Min")<<','<<Number(cfg,"mass2Max")
      <<"],\"mass2Bounds\":\"strict\",\"exactProductionTofWindow\":"<<Boolean(Flag(Integer(cfg,"exactProductionTofWindow")))
      <<",\"productionEnabled\":false,\"fullUncertaintyValidated\":false"
      <<",\"likelihoodPolicy\":\"Only explicitly selected pre-TPC TOF sample; no ALL or TAG addition; response seeds from the same sample are not independent priors\""
      <<",\"gatePolicy\":\"All four original TPC hypotheses compete first; strict nSigma, d<t<he3<he4 ties; target TOF/domain applied only after the original winner; no reassignment\""
      <<",\"responsePolicy\":\"He4 Gaussian mean/width free plus linear only; d/t/He3 amplitudes disabled by explicit user model; zeros are model assumptions, not measured absence\""
      <<",\"rawGateCountPolicy\":\"Fit-independent fractional native-dEdx overlap with full/touching bounds and uniform-R GL3\""
      <<",\"errorPolicy\":\"Conditional fixed-shape yield error and target shape covariance are separate, not complete purity uncertainty\",\"slices\":[";
  table<<"index\trBin\trLo\trHi\tspecies\tsignalSelected\ttotalSelected\tpurity\tconditionalError\trawGateCount\tlocalResidualFrac\tidentitySupport\tnumericalValid\tboundaryTarget\tselectionCovered\taxisCovered\n";
  counts<<"index\trBin\tlo\thi\tspecies\trawGateCount\trawGateCountLowerFullBins\trawGateCountUpperTouchingBins\n";
  std::set<int> used;
  for(unsigned i=0;i<slices.size();++i) {
    const Slice& slice=slices[i];
    Require(used.insert(slice.rBin).second&&slice.rBin<=source->GetNbinsX()&&std::fabs(source->GetXaxis()->GetBinLowEdge(slice.rBin)-slice.lo)<tol&&std::fabs(source->GetXaxis()->GetBinUpEdge(slice.rBin)-slice.hi)<tol&&std::fabs(slice.hi-slice.lo-Number(cfg,"nativeRigidityBinWidth"))<tol,"R slices must be exact unique native bins");
    Require(slice.lo>=Number(cfg,"rigidityMin")-tol&&slice.hi<=Number(cfg,"rigidityMax")+tol,"Slice outside explicit target domain");
    TH1D* raw=source->ProjectionY(Form("raw_native_%03d",slice.index),slice.rBin,slice.rBin,"e");raw->SetDirectory(0);
    const int first=raw->FindBin(slice.fitLo+tol),last=raw->FindBin(slice.fitHi-tol),rebin=Integer(cfg,"displayDedxRebinFactor");
    Require(first>=1&&last>=first&&last<=raw->GetNbinsX()&&rebin>0&&raw->GetNbinsX()%rebin==0&&(first-1)%rebin==0&&last%rebin==0&&std::fabs(raw->GetBinLowEdge(first)-slice.fitLo)<tol&&std::fabs(raw->GetXaxis()->GetBinUpEdge(last)-slice.fitHi)<tol,"Fit window not exact native/display bounds");
    typename E::Request request;request.nativeAxisMin=raw->GetXaxis()->GetXmin();request.nativeAxisMax=raw->GetXaxis()->GetXmax();request.backgroundModel=Text(cfg,"backgroundModel");request.backgroundSlope=Number(cfg,"backgroundStartSlope");request.freeBackgroundSlope=Flag(Integer(cfg,"freeBackgroundSlope"));
    for(int b=0;b<=raw->GetNbinsX()+1;++b) {
      double n=raw->GetBinContent(b),v=std::pow(raw->GetBinError(b),2);
      Require(Finite(n)&&n>=0&&Finite(v)&&std::fabs(v-n)<=Number(cfg,"countVarianceRelativeTolerance")*std::max(1.,n),"TAG cache count/error/flow is not Poisson-type");
      if(b>=1&&b<=raw->GetNbinsX()) Require(std::fabs(raw->GetBinWidth(b)-Number(cfg,"nativeDedxBinWidth"))<tol,"Unexpected native dEdx width");
      if(b>=first&&b<=last) request.bins.push_back(typename E::Bin(raw->GetBinLowEdge(b),raw->GetXaxis()->GetBinUpEdge(b),n));
    }
    RawGate measured[4];
    for(int s=0;s<4;++s) {
      measured[s]=Measure(slice,*raw,s,ns,tol);
      counts<<slice.index<<'\t'<<slice.rBin<<'\t'<<slice.lo<<'\t'<<slice.hi<<'\t'<<Species[s]<<'\t'<<measured[s].count<<'\t'<<measured[s].lower<<'\t'<<measured[s].upper<<'\n';
    }
    Vec acceptance(4,0),axisAcceptance(4,0),omittedUpper(4,0);std::vector<bool> omissionRatioAvailable(4,true);
    int active=0;double inAxisCount=raw->Integral(1,raw->GetNbinsX()),omittedSum=0;
    for(int s=0;s<4;++s) {
      typename E::Component c;c.name=Species[s];c.mean=slice.rows[4+s].mean;c.sigma=slice.rows[4+s].sigma;c.freeMean=s==target;c.freeSigma=s==target;
      for(int q=0;q<3;++q) {const Row& row=slice.rows[4*q+s];c.nodes.push_back(typename E::Node(row.weight,row.mean,row.sigma));acceptance[s]+=row.weight*E::Integral(slice.fitLo,slice.fitHi,row.mean,row.sigma);axisAcceptance[s]+=row.weight*E::Integral(request.nativeAxisMin,request.nativeAxisMax,row.mean,row.sigma);}
      c.active=s==target&&Finite(acceptance[s])&&acceptance[s]>floor;
      if(c.active) ++active;else if(s==target){bool ratio;omittedUpper[s]=OmittedBound(inAxisCount,acceptance[s],axisAcceptance[s],ratio);omissionRatioAvailable[s]=ratio;omittedSum+=omittedUpper[s];}
      request.components.push_back(c);
    }
    double count=0;for(unsigned b=0;b<request.bins.size();++b) count+=request.bins[b].count;
    const int bg=request.backgroundModel=="linear"?2:1;double per=count/(active+bg);
    for(int s=0;s<4;++s) request.components[s].startWindowYield=request.components[s].active?per:0;
    request.backgroundWindowYield=per*bg;
    typename E::Request identityReference=request;
    if(cfg.count("targetWarmStartFile")){std::ifstream seeds(Text(cfg,"targetWarmStartFile").c_str());Require(bool(seeds),"Missing high-R seeds");int index;double mean,sigma;while(seeds>>index>>mean>>sigma)if(index==slice.index){typename E::Component& c=request.components[target];double shift=mean-c.mean,ratio=sigma/c.sigma;for(unsigned q=0;q<c.nodes.size();++q){c.nodes[q].mean+=shift;c.nodes[q].sigma*=ratio;}c.mean=mean;c.sigma=sigma;}Require(seeds.eof(),"Malformed seeds");}
    typename E::Result fit;std::string failure;
    try {fit=E::Fit(request,cfg);}catch(const std::exception& e) {failure=e.what();}
    bool numerical=failure.empty()&&request.components[target].active&&fit.finite&&fit.amplitudeConverged&&fit.pointAccepted&&fit.status==0&&!fit.atBoundary;
    const TargetIdentity identity=DiagnoseTarget<E>(identityReference,fit,slice,target,cfg,true,numerical);
    if(!fit.shapeCovarianceValid) fit.shapeCovariance.clear();
    SaveSpectra<E>(product,*raw,request,fit,first,slice.index,rebin);
    if(i) json<<',';
    json<<"{\"index\":"<<slice.index<<",\"rBin\":"<<slice.rBin<<",\"lo\":"<<slice.lo<<",\"hi\":"<<slice.hi<<",\"r\":"<<.5*(slice.lo+slice.hi)
        <<",\"fitLo\":"<<slice.fitLo<<",\"fitHi\":"<<slice.fitHi<<",\"fitCount\":"<<count<<",\"observedNativeInAxisCount\":"<<inAxisCount
        <<",\"omittedLocalYieldUpperBoundSum\":"<<omittedSum<<",\"targetCenterShiftSigma\":"<<identity.shift<<",\"targetWidthRatio\":"<<identity.widthRatio
        <<",\"targetLabelCrossing\":"<<Boolean(identity.labelCrossing)<<",\"targetIdentityScreenPass\":"<<Boolean(identity.screenPass)<<",\"targetFreeResponseSupported\":"<<Boolean(identity.freeSupported)
        <<",\"shapeCovarianceValid\":"<<Boolean(fit.shapeCovarianceValid)<<",\"shapeCovariance\":";NuclearFit::WriteMatrix(json,fit.shapeCovariance);
    json<<",\"shapePoint\":";NuclearFit::WriteVec(json,fit.shapePoint);
    json<<",\"fitFailure\":"<<Quote(failure)<<",\"numericalValid\":"<<Boolean(numerical)<<",\"status\":"<<fit.status<<",\"amplitudeConverged\":"<<Boolean(fit.amplitudeConverged)
        <<",\"amplitudeKktResidual\":";Scalar(json,fit.amplitudeKkt);
    json<<",\"backgroundSlope\":";Scalar(json,fit.backgroundSlope);
    json<<",\"atBoundary\":"<<Boolean(fit.atBoundary)<<",\"dataDeviance\":";Scalar(json,fit.dataDeviance);
    json<<",\"pearson\":";Scalar(json,fit.pearson);
    json<<",\"conditionalNdf\":"<<fit.conditionalNdf<<",\"columnWindowYields\":";NuclearFit::WriteVec(json,fit.columnWindowYields);
    json<<",\"components\":{";
    for(int s=0;s<4;++s) {
      if(s) json<<',';
      double mean=slice.rows[4+s].mean,sigma=slice.rows[4+s].sigma;
      if(fit.components.size()==4) {mean=fit.components[s].mean;sigma=fit.components[s].sigma;}
      json<<Quote(Species[s])<<":{\"mean\":"<<mean<<",\"sigma\":"<<sigma<<",\"initialMean\":"<<slice.rows[4+s].mean<<",\"initialSigma\":"<<slice.rows[4+s].sigma<<",\"meanError\":";
      if(fit.meanError.size()>unsigned(s)&&fit.meanError[s]>=0) Scalar(json,fit.meanError[s]);else json<<"null";
      json<<",\"sigmaError\":";if(fit.sigmaError.size()>unsigned(s)&&fit.sigmaError[s]>=0) Scalar(json,fit.sigmaError[s]);else json<<"null";
      json<<",\"effectiveCount\":";if(fit.effectiveCount.size()>unsigned(s)) Scalar(json,fit.effectiveCount[s]);else json<<"null";
      json<<",\"conditionalDropDeltaD\":";if(fit.conditionalDropDeltaD.size()>unsigned(s)) Scalar(json,fit.conditionalDropDeltaD[s]);else json<<"null";
      json<<",\"isFreeTarget\":"<<Boolean(s==target)<<",\"active\":"<<Boolean(request.components[s].active)<<",\"initialResponseWindowAcceptance\":"<<acceptance[s]
          <<",\"responseAxisAcceptance\":"<<axisAcceptance[s]<<",\"omittedLocalYieldUpperBound\":"<<omittedUpper[s]<<",\"selection\":";
      Estimate estimate;std::vector<Regions> regions;
      if(fit.components.size()==4&&fit.columnWindowYields.size()>=4&&fit.totalCounts.size()==request.bins.size()) {
        estimate=Integrate<E>(fit,request,slice,*raw,first,s,ns,minimum,tol,&regions);
        Require(std::fabs(estimate.rawCount-measured[s].count)<=tol*std::max(1.,measured[s].count),"Fitted and independent raw gate counts disagree");
      } else {estimate.identitySupport=false;estimate.covered=false;estimate.axisCovered=false;}
      estimate.rawCount=measured[s].count;estimate.rawLower=measured[s].lower;estimate.rawUpper=measured[s].upper;
      WriteEstimate(json,estimate,numerical,regions,omittedSum,s==target&&identity.screenPass,s==target&&identity.freeSupported);json<<'}';
      table<<slice.index<<'\t'<<slice.rBin<<'\t'<<slice.lo<<'\t'<<slice.hi<<'\t'<<Species[s]<<'\t'<<estimate.signal<<'\t'<<estimate.total<<'\t';
      if(estimate.available) table<<estimate.purity;else table<<"nan";
      table<<'\t';if(estimate.errorAvailable) table<<estimate.error;else table<<"nan";
      table<<'\t'<<estimate.rawCount<<'\t'<<estimate.localResidual<<'\t'<<estimate.identitySupport<<'\t'<<numerical<<'\t'<<estimate.boundaryTarget<<'\t'<<estimate.covered<<'\t'<<estimate.axisCovered<<'\n';
    }
    json<<"}}";std::cout<<"NUCLEAR_HIGH_R_SLICE sample="<<sample<<" target="<<targetName<<" R="<<.5*(slice.lo+slice.hi)<<" numerical="<<numerical<<std::endl;delete raw;
  }
  json<<"],\"executionIsNotPhysicalValidation\":true}\n";json.close();table.close();counts.close();
  Require(bool(json)&&bool(table)&&bool(counts),"Cannot close output tables");
  product.Close();Readback((folder+"/cutoff_fits.root").Data(),*source,slices,Number(cfg,"countVarianceRelativeTolerance"));cache.Close();
  std::cout<<"NUCLEAR_HIGH_R_COMPLETE rows="<<slices.size()<<std::endl;
}
}
#endif
void FitHighRigidity(const char* cachePath,const char* output,const char* settings,const char* shapes,const char* gates) {
#ifndef __CINT__
  try {
    NuclearFit::Settings cfg=NuclearFit::ReadSettings(settings);
    std::vector<NuclearLocalCutoff::Slice> slices=NuclearLocalCutoff::ReadInputs(shapes,gates,NuclearFit::Number(cfg,"axisTolerance"));
    if(NuclearFit::Text(cfg,"responseModel")=="gaussian") NuclearHighRigidity::Run<NuclearLocalCutoff::GaussianEngine>(cachePath,output,cfg,slices);
    else if(NuclearFit::Text(cfg,"responseModel")=="lognormal") NuclearHighRigidity::Run<NuclearLocalCutoff::LognormalEngine>(cachePath,output,cfg,slices);
    else NuclearFit::Require(false,"Unknown response model");
  }catch(const std::exception& e) {std::cerr<<"FitHighRigidity: "<<e.what()<<std::endl;gSystem->Exit(2);}
#endif
}
