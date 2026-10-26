// Study-only independent Gaussian core/tail columns. No PID update or purity estimate.
void FitFreeCoreTailPilot(const char*,const char*,const char*,const char*,const char*);
#ifndef __CINT__
#include "../FitLocalNuclearCutoff.C"
#include "FreeGaussianCoreTail.h"

namespace FreeCoreTailPilot {
using namespace NuclearLocalCutoff;
namespace Free=NuclearFreeCoreTail;
namespace Fixed=NuclearCoreTailFit;
typedef NuclearFit::Request Request;
typedef NuclearFit::Result Result;

std::set<int> Names(const Settings& cfg,const char* key){
  std::istringstream input(Text(cfg,key));std::string name;std::set<int> out;
  while(std::getline(input,name,',')){int s=SpeciesIndex(name);Require(s>=0&&out.insert(s).second,"Invalid/duplicate named free-tail species");}
  Require(!out.empty(),"Empty explicit free-tail species policy");return out;
}
void Strings(std::ostream& out,const std::vector<std::string>& values){
  out<<'[';for(unsigned i=0;i<values.size();++i){if(i)out<<',';out<<Quote(values[i]);}out<<']';
}
void Numbers(std::ostream& out,const Vec& values){
  out<<'[';for(unsigned i=0;i<values.size();++i){if(i)out<<',';Scalar(out,values[i]);}out<<']';
}
void Matrix(std::ostream& out,const NuclearFit::Matrix& values){
  out<<'[';for(unsigned i=0;i<values.size();++i){if(i)out<<',';Numbers(out,values[i]);}out<<']';
}
void Booleans(std::ostream& out,const std::vector<bool>& values){
  out<<'[';for(unsigned i=0;i<values.size();++i){if(i)out<<',';out<<Boolean(values[i]);}out<<']';
}
std::vector<Fixed::Shape> Shapes(const Settings& cfg,const std::set<int>& tail){
  const Fixed::Shape seed(Number(cfg,"tailFraction"),Number(cfg,"tailOffsetSeedSigma"),Number(cfg,"tailWidthRatio"));
  Fixed::ValidateShape(seed);std::vector<Fixed::Shape> out;
  for(int s=0;s<4;++s)out.push_back(tail.count(s)?seed:Fixed::Shape(0,0,1));return out;
}

// This is an observed frozen-gate diagnostic, deliberately independent of fitted yield.
Estimate RawGate(const Request& req,const Slice& slice,const TH1D& raw,const Result& fit,
                 int first,int target,double ns,double tol,std::vector<Regions>& regions){
  Estimate out;regions.clear();
  for(int q=0;q<3;++q){
    std::vector<Choice> choices;for(int s=0;s<4;++s)choices.push_back(slice.rows[4*q+s].gate);
    regions.push_back(SelectionRegions(choices,target,ns,req.nativeAxisMin,req.nativeAxisMax));
    const Row& row=slice.rows[4*q+target];out.identitySupport=out.identitySupport&&row.supported;
    out.gateActive=out.gateActive||row.gate.active;
    if(row.gate.active)out.axisCovered=out.axisCovered&&row.gate.mean-ns*row.gate.sigma>=req.nativeAxisMin&&row.gate.mean+ns*row.gate.sigma<=req.nativeAxisMax;
    for(unsigned k=0;k<regions[q].size();++k)out.covered=out.covered&&regions[q][k].lo>=slice.fitLo-tol&&regions[q][k].hi<=slice.fitHi+tol;
  }
  const bool prediction=fit.totalCounts.size()==req.bins.size();
  for(int b=1;b<=raw.GetNbinsX();++b){
    double fraction=0,full=0,touch=0;
    for(int q=0;q<3;++q){const double overlap=Overlap(regions[q],raw.GetBinLowEdge(b),raw.GetXaxis()->GetBinUpEdge(b)),w=slice.rows[4*q].weight;fraction+=w*overlap;if(overlap>=1-tol)full+=w;if(overlap>tol)touch+=w;}
    const double n=raw.GetBinContent(b);out.rawCount+=fraction*n;out.rawLower+=full*n;out.rawUpper+=touch*n;
    if(prediction&&b>=first&&b<first+int(fit.totalCounts.size()))out.absoluteResidual+=fraction*std::fabs(n-fit.totalCounts[b-first]);
  }
  out.localResidual=out.rawCount>0?out.absoluteResidual/out.rawCount:0;return out;
}
void UnavailableSelection(std::ostream& out,const Estimate& e,const std::vector<Regions>& regions,
                          bool numerical,bool prediction,double omitted){
  out<<"{\"available\":false,\"purity\":null,\"signalSelected\":null,\"totalSelected\":null,\"backgroundSelected\":null,\"conditionalError\":null,\"purityLowerIncludingOmittedBound\":null,\"selectedColumnFractions\":null,\"selectedColumnCounts\":null,\"reason\":\"group_purity_not_evaluated_in_response_quality_pilot\",\"rawGateCount\":";Scalar(out,e.rawCount);
  out<<",\"rawGateCountLowerFullBins\":";Scalar(out,e.rawLower);out<<",\"rawGateCountUpperTouchingBins\":";Scalar(out,e.rawUpper);
  out<<",\"gateAbsoluteResidualCount\":";if(prediction)Scalar(out,e.absoluteResidual);else out<<"null";
  out<<",\"localResidualFrac\":";if(prediction&&e.rawCount>0)Scalar(out,e.localResidual);else out<<"null";
  out<<",\"omittedLocalYieldUpperBound\":";Scalar(out,omitted);
  out<<",\"omissionBoundOverRawGateCount\":";if(e.rawCount>0)Scalar(out,omitted/e.rawCount);else out<<"null";
  out<<",\"numericalValid\":"<<Boolean(numerical)<<",\"identitySupport\":"<<Boolean(e.identitySupport)<<",\"identitySupportMeaning\":\"Frozen reference support, not independent isotope validation\",\"selectionCovered\":"<<Boolean(e.covered)<<",\"axisCovered\":"<<Boolean(e.axisCovered)<<",\"gateActive\":"<<Boolean(e.gateActive)<<",\"gateIntervalsAtNodes\":";WriteRegions(out,regions);out<<'}';
}
void SaveDynamicColumns(TFile& product,const TH1D& raw,const Result& fit,int first,int index){
  product.cd();
  for(unsigned c=0;c<fit.components.size();++c){
    TH1D* h=static_cast<TH1D*>(raw.Clone(Form("prediction_freecolumn_%02d_%03d",c,index)));h->SetDirectory(0);h->Reset();h->SetTitle(fit.components[c].name.c_str());
    for(unsigned b=0;b<fit.componentCounts.size();++b){Require(fit.componentCounts[b].size()==fit.components.size(),"Bad dynamic prediction dimensions");h->SetBinContent(first+b,fit.componentCounts[b][c]);}
    Require(h->Write()>0,"Cannot save independent Gaussian column");delete h;
  }
}
void DynamicReadback(const char* path,const std::vector<Slice>& slices,
                     const std::vector<Free::Result>& saved,const std::vector<int>& first,double tol){
  Require(slices.size()==saved.size()&&saved.size()==first.size(),"Missing dynamic readback reference");
  TFile input(path,"READ");Require(!input.IsZombie()&&!input.TestBit(TFile::kRecovered),"Cannot reopen dynamic prediction ROOT");
  for(unsigned i=0;i<slices.size();++i){
    const Free::Result& result=saved[i];const Result& model=result.model;
    std::vector<TH1*> columns;
    for(unsigned c=0;c<model.components.size();++c){
      TH1* h=dynamic_cast<TH1*>(input.Get(Form("prediction_freecolumn_%02d_%03d",c,slices[i].index)));
      Require(h!=0,"Missing independent Gaussian prediction in closed ROOT");columns.push_back(h);
      for(int b=0;b<=h->GetNbinsX()+1;++b){
        const int local=b-first[i];const double expected=local>=0&&local<int(model.componentCounts.size())?model.componentCounts[local][c]:0;
        Require(Finite(h->GetBinContent(b))&&std::fabs(h->GetBinContent(b)-expected)<=tol*std::max(1.,std::fabs(expected)),"Closed ROOT changed independent Gaussian count/flow");
      }
    }
    for(unsigned s=0;s<result.isotopeMap.size();++s){
      const Free::IsotopeColumns& map=result.isotopeMap[s];
      Require(map.core<columns.size()&&(map.tail<0||unsigned(map.tail)<columns.size()),"Invalid dynamic readback isotope map");
      TH1* grouped=dynamic_cast<TH1*>(input.Get(Form("prediction_%s_%03d",map.name.c_str(),slices[i].index)));
      Require(grouped!=0,"Missing folded isotope prediction in closed ROOT");
      for(int b=0;b<=grouped->GetNbinsX()+1;++b){
        double sum=columns[map.core]->GetBinContent(b);if(map.tail>=0)sum+=columns[map.tail]->GetBinContent(b);
        Require(Finite(sum)&&std::fabs(sum-grouped->GetBinContent(b))<=tol*std::max(1.,std::fabs(sum)),"Closed ROOT dynamic core/tail sum differs from folded isotope");
      }
    }
  }
  input.Close();std::cout<<"FREE_CORETAIL_DYNAMIC_READBACK_PASS rows="<<slices.size()<<std::endl;
}
void Run(const char* cachePath,const char* output,const Settings& cfg,const std::vector<Slice>& slices){
  const double tol=Number(cfg,"axisTolerance"),floor=Number(cfg,"minimumGaussianAcceptance"),ns=Number(cfg,"maxNSigma");
  const int target=SpeciesIndex(Text(cfg,"targetSpecies"));
  Require(target>=0&&Text(cfg,"sampleHistogram")=="ALL"&&Text(cfg,"responseModel")=="gaussian_coretail_free_components"&&Text(cfg,"responseStudyMode")=="free_coretail","Explicit free-coretail ALL response declaration required");
  const std::set<int> free=Names(cfg,"freeSpecies"),tail=Names(cfg,"tailSpecies");
  Require(free.count(target)&&free==tail,"Tail hypotheses must match the explicitly free core hypotheses");
  const std::vector<Fixed::Shape> shapes=Shapes(cfg,tail);
  TString directory(output);Require(!gSystem->AccessPathName(directory),"Pilot directory must be provenance-prepared");
  const char* products[]={"cutoff_fits.json","cutoff_fits.root"};for(int k=0;k<2;++k)Require(gSystem->AccessPathName(directory+"/"+products[k]),"Existing free-pilot product will not be replaced");
  TFile cache(cachePath,"READ");Require(!cache.IsZombie()&&!cache.TestBit(TFile::kRecovered),"Cannot read native ALL cache");TH2* source=dynamic_cast<TH2*>(cache.Get("ALL"));Require(source!=0,"Missing native ALL");
  TFile product(directory+"/cutoff_fits.root","CREATE");Require(!product.IsZombie(),"Cannot CREATE free-pilot ROOT");
  std::ofstream json((directory+"/cutoff_fits.json").Data());Require(bool(json),"Cannot create free-pilot JSON");json<<std::setprecision(17);
  json<<"{\"schemaVersion\":1,\"responseModel\":\"gaussian_coretail_free_components\",\"responseStudyMode\":\"free_coretail\",\"backgroundModel\":"<<Quote(Text(cfg,"backgroundModel"))<<",\"targetSpecies\":"<<Quote(Text(cfg,"targetSpecies"))<<",\"sampleHistogram\":\"ALL\",\"productionEnabled\":false,\"fullUncertaintyValidated\":false,\"purityPolicy\":\"Not estimated in this response-quality pilot; null is not zero\",\"freeSpecies\":"<<Quote(Text(cfg,"freeSpecies"))<<",\"tailSpecies\":"<<Quote(Text(cfg,"tailSpecies"))<<",\"profileName\":"<<Quote(Text(cfg,"profileName"))<<",\"tailFraction\":"<<Number(cfg,"tailFraction")<<",\"tailOffsetSeedSigma\":"<<Number(cfg,"tailOffsetSeedSigma")<<",\"tailWidthRatio\":"<<Number(cfg,"tailWidthRatio")<<",\"responsePolicy\":\"Independent core/tail Gaussian means, sigmas and nonnegative yields; seeds are not fitted constraints; isotope grouping is only a post-fit screen\",\"gatePolicy\":\"Exact frozen old TPC gates and fit windows; no re-selection, purity evaluation or coefficient adoption\",\"qualityPolicy\":\"Seed/model selection is not included in an uncorrected goodness-of-fit p-value; raw dynamic NDF and covariance are retained\",\"slices\":[";
  std::set<int> used;std::vector<Free::Result> savedResults;std::vector<int> savedFirst;
  for(unsigned i=0;i<slices.size();++i){
    const Slice& slice=slices[i];Require(used.insert(slice.rBin).second&&slice.rBin<=source->GetNbinsX()&&std::fabs(source->GetXaxis()->GetBinLowEdge(slice.rBin)-slice.lo)<tol&&std::fabs(source->GetXaxis()->GetBinUpEdge(slice.rBin)-slice.hi)<tol,"Free-pilot R bin differs from source");
    TH1D* raw=source->ProjectionY(Form("raw_native_%03d",slice.index),slice.rBin,slice.rBin,"e");raw->SetDirectory(0);
    const int first=raw->FindBin(slice.fitLo+tol),last=raw->FindBin(slice.fitHi-tol),rebin=Integer(cfg,"displayDedxRebinFactor");
    Require(first>=1&&last>=first&&last<=raw->GetNbinsX()&&rebin>0&&raw->GetNbinsX()%rebin==0&&std::fabs(raw->GetBinLowEdge(first)-slice.fitLo)<tol&&std::fabs(raw->GetXaxis()->GetBinUpEdge(last)-slice.fitHi)<tol,"Free-pilot window changed native bounds");
    Request req;req.nativeAxisMin=raw->GetXaxis()->GetXmin();req.nativeAxisMax=raw->GetXaxis()->GetXmax();req.backgroundModel=Text(cfg,"backgroundModel");req.backgroundSlope=Number(cfg,"backgroundStartSlope");req.freeBackgroundSlope=Flag(Integer(cfg,"freeBackgroundSlope"));
    for(int b=0;b<=raw->GetNbinsX()+1;++b){
      const double n=raw->GetBinContent(b),v=std::pow(raw->GetBinError(b),2);
      Require(Finite(n)&&n>=0&&Finite(v)&&std::fabs(v-n)<=Number(cfg,"countVarianceRelativeTolerance")*std::max(1.,n),"Free-pilot cache is not Poisson native count/error/flow");
      if(b>=1&&b<=raw->GetNbinsX())Require(std::fabs(raw->GetBinWidth(b)-Number(cfg,"nativeDedxBinWidth"))<tol,"Free-pilot dEdx width changed");
      if(b>=first&&b<=last)req.bins.push_back(NuclearFit::Bin(raw->GetBinLowEdge(b),raw->GetXaxis()->GetBinUpEdge(b),n));
    }
    Vec acceptance(4,0),axisAcceptance(4,0),omitted(4,0);std::vector<bool> omissionRatio(4,true);double omittedSum=0;int active=0;
    for(int s=0;s<4;++s){
      NuclearFit::Component c;c.name=Species[s];c.mean=slice.rows[4+s].mean;c.sigma=slice.rows[4+s].sigma;c.freeMean=free.count(s);c.freeSigma=free.count(s);
      for(int q=0;q<3;++q){const Row& r=slice.rows[4*q+s];c.nodes.push_back(NuclearFit::Node(r.weight,r.mean,r.sigma));}
      acceptance[s]=NuclearFit::KernelIntegral(c,slice.fitLo,slice.fitHi);axisAcceptance[s]=NuclearFit::KernelIntegral(c,req.nativeAxisMin,req.nativeAxisMax);c.active=Finite(acceptance[s])&&acceptance[s]>floor;
      if(c.active)++active;else{bool ratio;omitted[s]=OmittedBound(raw->Integral(1,raw->GetNbinsX()),acceptance[s],axisAcceptance[s],ratio);omissionRatio[s]=ratio;omittedSum+=omitted[s];}
      req.components.push_back(c);
    }
    const double count=NuclearFit::TotalCount(req.bins),per=count/(active+(req.backgroundModel=="linear"?2:1));
    for(int s=0;s<4;++s)req.components[s].startWindowYield=req.components[s].active?per:0;req.backgroundWindowYield=per*(req.backgroundModel=="linear"?2:1);
    Free::Result result;Result folded;std::string failure;
    try{result=Free::Fit(req,cfg,shapes);}catch(const std::exception& e){failure=e.what();}
    // A folding defect is an output-integrity failure, not a zero-valued fitted model.
    if(!result.model.components.empty())folded=Free::FoldIsotopes(result);
    const Result& fit=result.model;const bool numerical=failure.empty()&&req.components[target].active&&fit.finite&&fit.amplitudeConverged&&fit.pointAccepted&&fit.status==0&&!fit.atBoundary;
    const bool prediction=fit.totalCounts.size()==req.bins.size();
    SaveSpectra<GaussianEngine>(product,*raw,req,folded,first,slice.index,rebin);SaveDynamicColumns(product,*raw,fit,first,slice.index);
    savedResults.push_back(result);savedFirst.push_back(first);
    unsigned activeYields=0;for(unsigned j=0;j<fit.yieldCovarianceColumnActive.size();++j)activeYields+=fit.yieldCovarianceColumnActive[j];
    if(prediction)Require(fit.conditionalNdf==int(req.bins.size())-int(fit.shapeNames.size())-int(activeYields),"Free-pilot NDF must count all dynamic free shapes and interior yields");
    if(i)json<<',';
    json<<"{\"index\":"<<slice.index<<",\"rBin\":"<<slice.rBin<<",\"lo\":"<<slice.lo<<",\"hi\":"<<slice.hi<<",\"r\":"<<.5*(slice.lo+slice.hi)<<",\"fitLo\":"<<slice.fitLo<<",\"fitHi\":"<<slice.fitHi<<",\"fitCount\":"<<count<<",\"observedNativeInAxisCount\":"<<raw->Integral(1,raw->GetNbinsX())<<",\"omittedLocalYieldUpperBoundSum\":"<<omittedSum<<",\"fitFailure\":"<<Quote(failure)<<",\"numericalValid\":"<<Boolean(numerical)<<",\"isotopeGroupingSupported\":"<<Boolean(numerical&&result.groupingSupported)<<",\"status\":"<<fit.status<<",\"covarianceStatus\":"<<fit.covarianceStatus<<",\"atBoundary\":"<<Boolean(fit.atBoundary)<<",\"pointAccepted\":"<<Boolean(fit.pointAccepted)<<",\"amplitudeConverged\":"<<Boolean(fit.amplitudeConverged)<<",\"amplitudeKktResidual\":";Scalar(json,fit.amplitudeKkt);
    json<<",\"edm\":";Scalar(json,fit.edm);json<<",\"shapeCovarianceValid\":"<<Boolean(fit.shapeCovarianceValid)<<",\"yieldCovarianceValid\":"<<Boolean(fit.yieldCovarianceValid)<<",\"pearson\":";Scalar(json,fit.pearson);json<<",\"dataDeviance\":";Scalar(json,fit.dataDeviance);
    json<<",\"conditionalNdf\":"<<fit.conditionalNdf<<",\"freeShapeParameterCount\":"<<fit.shapeNames.size()<<",\"interiorYieldParameterCount\":"<<activeYields<<",\"shapeNames\":";Strings(json,fit.shapeNames);json<<",\"shapePoint\":";Numbers(json,fit.shapePoint);json<<",\"shapeCovariance\":";Matrix(json,fit.shapeCovariance);
    json<<",\"columnNames\":";Strings(json,fit.columnNames);json<<",\"columnWindowYields\":";Numbers(json,fit.columnWindowYields);json<<",\"yieldCovarianceColumnActive\":";Booleans(json,fit.yieldCovarianceColumnActive);json<<",\"fixedShapeYieldCovariance\":";Matrix(json,fit.fixedShapeYieldCovariance);
    json<<",\"backgroundSlope\":";Scalar(json,fit.backgroundSlope);json<<",\"postFitFlags\":";Strings(json,result.flags);json<<",\"isotopeMap\":[";
    for(unsigned s=0;s<result.isotopeMap.size();++s){const Free::IsotopeColumns& map=result.isotopeMap[s];if(s)json<<',';json<<"{\"name\":"<<Quote(map.name)<<",\"coreColumn\":"<<map.core<<",\"tailColumn\":";if(map.tail>=0)json<<map.tail;else json<<"null";json<<",\"active\":"<<Boolean(map.active)<<",\"referenceMean\":";Scalar(json,map.referenceMean);json<<",\"referenceSigma\":";Scalar(json,map.referenceSigma);json<<'}';}
    json<<"],\"dynamicComponents\":[";
    for(unsigned c=0;c<fit.components.size();++c){
      const NuclearFit::Component& model=fit.components[c];if(c)json<<',';json<<"{\"name\":"<<Quote(model.name)<<",\"active\":"<<Boolean(model.active)<<",\"freeMean\":"<<Boolean(model.freeMean)<<",\"freeSigma\":"<<Boolean(model.freeSigma)<<",\"mean\":";Scalar(json,model.mean);json<<",\"sigma\":";Scalar(json,model.sigma);
      json<<",\"meanError\":";if(c<fit.meanError.size()&&fit.meanError[c]>=0)Scalar(json,fit.meanError[c]);else json<<"null";json<<",\"sigmaError\":";if(c<fit.sigmaError.size()&&fit.sigmaError[c]>=0)Scalar(json,fit.sigmaError[c]);else json<<"null";
      const double a=c<fit.windowAcceptance.size()?fit.windowAcceptance[c]:0,y=c<fit.columnWindowYields.size()?fit.columnWindowYields[c]:0;
      json<<",\"windowAcceptance\":";if(c<fit.windowAcceptance.size())Scalar(json,a);else json<<"null";json<<",\"axisAcceptance\":";Scalar(json,NuclearFit::KernelIntegral(model,req.nativeAxisMin,req.nativeAxisMax));json<<",\"windowYield\":";if(c<fit.columnWindowYields.size())Scalar(json,y);else json<<"null";json<<",\"fullLineYield\":";if(a>0&&c<fit.columnWindowYields.size())Scalar(json,y/a);else json<<"null";
      json<<",\"rNodes\":[";for(unsigned q=0;q<model.nodes.size();++q){if(q)json<<',';json<<"{\"weight\":";Scalar(json,model.nodes[q].weight);json<<",\"mean\":";Scalar(json,model.nodes[q].mean);json<<",\"sigma\":";Scalar(json,model.nodes[q].sigma);json<<'}';}json<<"]}";
    }
    json<<"],\"components\":{";
    for(int s=0;s<4;++s){
      const bool fitted=folded.components.size()==4;const NuclearFit::Component& core=fitted?folded.components[s]:req.components[s];
      std::vector<Regions> regions;const Estimate gate=RawGate(req,slice,*raw,fit,first,s,ns,tol,regions);
      if(s)json<<',';json<<Quote(Species[s])<<":{\"mean\":";Scalar(json,core.mean);json<<",\"sigma\":";Scalar(json,core.sigma);json<<",\"initialMean\":"<<req.components[s].mean<<",\"initialSigma\":"<<req.components[s].sigma<<",\"coreMean\":";Scalar(json,core.mean);json<<",\"coreSigma\":";Scalar(json,core.sigma);
      json<<",\"active\":"<<Boolean(req.components[s].active)<<",\"freeMean\":"<<Boolean(req.components[s].freeMean)<<",\"freeSigma\":"<<Boolean(req.components[s].freeSigma)<<",\"initialResponseWindowAcceptance\":"<<acceptance[s]<<",\"responseAxisAcceptance\":"<<axisAcceptance[s]<<",\"omittedLocalYieldUpperBound\":"<<omitted[s]<<",\"omissionRatioAvailable\":"<<Boolean(omissionRatio[s]);
      json<<",\"fullLineTailFraction\":";if(result.fullLineTailFraction.size()>unsigned(s)&&result.fullLineTailFraction[s]>=0)Scalar(json,result.fullLineTailFraction[s]);else json<<"null";
      json<<",\"windowTailFraction\":";if(result.windowTailFraction.size()>unsigned(s)&&result.windowTailFraction[s]>=0)Scalar(json,result.windowTailFraction[s]);else json<<"null";
      json<<",\"postFitFlags\":";if(result.isotopeFlags.size()>unsigned(s))Strings(json,result.isotopeFlags[s]);else Strings(json,std::vector<std::string>(1,"fit_unavailable"));
      json<<",\"selection\":";UnavailableSelection(json,gate,regions,numerical,prediction,omittedSum);json<<'}';
    }
    json<<"}}";delete raw;
    std::cout<<"FREE_CORETAIL_PILOT_SLICE index="<<slice.index<<" numerical="<<numerical<<" chi2="<<fit.pearson<<" ndf="<<fit.conditionalNdf<<std::endl;
  }
  json<<"],\"executionIsNotPhysicalValidation\":true}\n";json.close();Require(bool(json),"Cannot finish free-pilot JSON");product.Close();
  Readback((directory+"/cutoff_fits.root").Data(),*source,slices,Number(cfg,"countVarianceRelativeTolerance"));
  DynamicReadback((directory+"/cutoff_fits.root").Data(),slices,savedResults,savedFirst,Number(cfg,"countVarianceRelativeTolerance"));cache.Close();
  std::cout<<"CORETAIL_PILOT_COMPLETE rows="<<slices.size()<<std::endl;
}
}
#endif
void FitFreeCoreTailPilot(const char* cachePath,const char* output,const char* settings,const char* shapes,const char* gates){
#ifndef __CINT__
  try{NuclearFit::Settings cfg=NuclearFit::ReadSettings(settings);std::vector<NuclearLocalCutoff::Slice> slices=NuclearLocalCutoff::ReadInputs(shapes,gates,NuclearFit::Number(cfg,"axisTolerance"));FreeCoreTailPilot::Run(cachePath,output,cfg,slices);}catch(const std::exception& error){std::cerr<<"FitFreeCoreTailPilot: "<<error.what()<<std::endl;gSystem->Exit(2);}
#endif
}
