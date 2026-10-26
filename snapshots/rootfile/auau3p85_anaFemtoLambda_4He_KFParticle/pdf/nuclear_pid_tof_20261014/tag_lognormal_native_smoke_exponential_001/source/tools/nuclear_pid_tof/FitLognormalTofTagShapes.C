void FitLognormalTofTagShapes(const char*,const char*,const char*,const char*,const char*);
#ifndef __CINT__
#include "NuclearLognormalFitCommon.h"
#include "TFile.h"
#include "TH2.h"
#include <set>
namespace NuclearLognormalTag {
using namespace NuclearLognormalFit;
struct Slice {int index;std::string sample;double lo,hi,fitLo,fitHi;std::vector<std::string> names;};
struct Seed {Component component;double backgroundSlope,backgroundYield;};
std::vector<std::string> Split(const std::string& value){std::vector<std::string> out;std::istringstream in(value);std::string item;while(std::getline(in,item,',')){Require(!item.empty(),"Empty component name");out.push_back(item);}return out;}
std::vector<Slice> ReadSlices(const char* path){std::ifstream file(path);Require(bool(file),"Missing R slices TSV");std::vector<Slice> rows;std::set<int> seen;std::string line;
  while(std::getline(file,line)){std::istringstream in(line);Slice s;std::string components,extra;
    Require(bool(in>>s.index>>s.sample>>s.lo>>s.hi>>s.fitLo>>s.fitHi>>components)&&!(in>>extra),"Malformed 7-field R slice");
    Require(s.index>=0&&seen.insert(s.index).second&&(s.sample=="TAG_D"||s.sample=="TAG_T"||s.sample=="TAG_HE3")&&Finite(s.lo)&&Finite(s.hi)&&Finite(s.fitLo)&&Finite(s.fitHi)&&s.lo>0&&s.hi>s.lo&&s.fitHi>s.fitLo,"Invalid R slice");
    s.names=Split(components);std::set<std::string> names;for(unsigned j=0;j<s.names.size();++j)Require((s.names[j]=="d"||s.names[j]=="t"||s.names[j]=="he3"||s.names[j]=="he4")&&names.insert(s.names[j]).second,"Unknown or duplicate component");rows.push_back(s);}
  Require(file.eof()&&!rows.empty(),"Empty or unreadable R slices");return rows;}
std::map<std::string,Seed> ReadSeeds(const char* path){std::ifstream file(path);Require(bool(file),"Missing seeds TSV");std::map<std::string,Seed> out;std::string line;
  while(std::getline(file,line)){std::istringstream in(line);int index,active,fm,fs;Seed seed;std::string extra;
    Require(bool(in>>index>>seed.component.name>>seed.component.mean>>seed.component.sigma>>seed.component.startWindowYield>>active>>fm>>fs)&&!(in>>extra),"Malformed 8-field seed row");
    Require((active==0||active==1)&&(fm==0||fm==1)&&(fs==0||fs==1)&&Finite(seed.component.mean)&&Finite(seed.component.sigma)&&Finite(seed.component.startWindowYield)&&seed.component.startWindowYield>=0,"Invalid component seed");
    seed.component.active=active;seed.component.freeMean=fm;seed.component.freeSigma=fs;
    std::ostringstream key;key<<index<<':'<<seed.component.name;Require(out.insert(std::make_pair(key.str(),seed)).second,"Duplicate seed");}
  Require(file.eof(),"Unparsed seed data");return out;}
Seed NeedSeed(const std::map<std::string,Seed>& seeds,int index,const std::string& name){std::ostringstream key;key<<index<<':'<<name;std::map<std::string,Seed>::const_iterator i=seeds.find(key.str());Require(i!=seeds.end(),"Missing component seed "+key.str());return i->second;}
void WriteResult(std::ostream& out,const Result& r,const Request& request){
  out<<"\"dataDeviance\":"<<r.dataDeviance<<",\"chi2\":"<<r.pearson<<",\"conditionalNdf\":"<<r.conditionalNdf
    <<",\"status\":"<<r.status<<",\"covarianceStatus\":"<<r.covarianceStatus<<",\"edm\":"<<r.edm
    <<",\"finite\":"<<(r.finite?"true":"false")<<",\"pointAccepted\":"<<(r.pointAccepted?"true":"false")
    <<",\"amplitudeConverged\":"<<(r.amplitudeConverged?"true":"false")<<",\"amplitudeKkt\":"<<r.amplitudeKkt<<",\"amplitudeIterations\":"<<r.amplitudeIterations
    <<",\"atBoundary\":"<<(r.atBoundary?"true":"false")<<",\"shapeCovarianceValid\":"<<(r.shapeCovarianceValid?"true":"false")
    <<",\"fixedShapeYieldCovarianceValid\":"<<(r.yieldCovarianceValid?"true":"false")<<",\"backgroundSlope\":"<<r.backgroundSlope
    <<",\"columnWindowYields\":";WriteVec(out,r.columnWindowYields);
  out<<",\"columnNames\":[";for(unsigned i=0;i<r.columnNames.size();++i){if(i)out<<',';out<<Quote(r.columnNames[i]);}out<<"]";
  out<<",\"yieldCovarianceColumnActive\":[";for(unsigned i=0;i<r.yieldCovarianceColumnActive.size();++i){if(i)out<<',';out<<(r.yieldCovarianceColumnActive[i]?"true":"false");}out<<"]";
  out<<",\"shapePoint\":";WriteVec(out,r.shapePoint);out<<",\"shapeParameterNames\":[";for(unsigned i=0;i<r.shapeNames.size();++i){if(i)out<<',';out<<Quote(r.shapeNames[i]);}out<<"]";
  out<<",\"shapeCovariance\":";WriteMatrix(out,r.shapeCovariance);out<<",\"fixedShapeYieldCovariance\":";WriteMatrix(out,r.fixedShapeYieldCovariance);
  out<<",\"covarianceScope\":\"shape: local mode/logCurvatureWidth/expSlope, amplitudes profiled at retained point; yield: final shapes fixed, interior active columns only; NOT joint TAG/ALL or transfer uncertainty\""
    <<",\"dropDeltaDScope\":\"Remove one component and reprofile amplitudes with all retained shapes fixed; NOT a significance or fully refitted null test\""
    <<",\"components\":{\"";
  for(unsigned j=0;j<r.components.size();++j){const Component& c=r.components[j];if(j)out<<",\"";out<<c.name<<"\":{\"mean\":"<<c.mean<<",\"sigma\":"<<c.sigma<<",\"mode\":"<<c.mean<<",\"localCurvatureWidth\":"<<c.sigma
      <<",\"windowYield\":"<<r.columnWindowYields[j]<<",\"windowAcceptance\":"<<r.windowAcceptance[j]<<",\"effectiveCount\":"<<r.effectiveCount[j]
      <<",\"meanError\":"<<r.meanError[j]<<",\"sigmaError\":"<<r.sigmaError[j]<<",\"active\":"<<(c.active?"true":"false")
      <<",\"freeMean\":"<<(c.freeMean?"true":"false")<<",\"freeSigma\":"<<(c.freeSigma?"true":"false")
      <<",\"detectorAcceptance\":"<<KernelIntegral(c,request.nativeAxisMin,request.nativeAxisMax)<<",\"centerOutsideDetector\":"<<((c.mean<request.nativeAxisMin||c.mean>request.nativeAxisMax)?"true":"false")
      <<",\"zeroYield\":"<<(r.columnWindowYields[j]<=TotalCount(request.bins)*1.e-12?"true":"false")<<",\"conditionalDropDeltaD\":"<<r.conditionalDropDeltaD[j]
      <<",\"identityValidated\":false,\"responseModel\":\"lognormal\",\"parameterMeaning\":\"mean alias is mode; sigma alias is local curvature width, NOT arithmetic mean or standard deviation\"}";}
  out<<"}";
}
}
#endif
void FitLognormalTofTagShapes(const char* cacheRoot,const char* output,const char* settings,const char* rSlices,const char* seeds){
#ifndef __CINT__
  using namespace NuclearLognormalTag;
  try{
    Settings cfg=ReadSettings(settings);ValidateSettings(cfg);std::vector<Slice> slices=ReadSlices(rSlices);std::map<std::string,Seed> starts=ReadSeeds(seeds);
    unsigned expectedSeeds=0;for(unsigned i=0;i<slices.size();++i){expectedSeeds+=slices[i].names.size()+1;NeedSeed(starts,slices[i].index,"BG");for(unsigned j=0;j<slices[i].names.size();++j)NeedSeed(starts,slices[i].index,slices[i].names[j]);}
    Require(starts.size()==expectedSeeds,"Unused component/BG seed rows are forbidden");
    TString folder(output);const char* outputs[]={"tag_fits.root","tag_fits.json","slices.tsv","completed.json"};
    for(unsigned i=0;i<sizeof(outputs)/sizeof(outputs[0]);++i)Require(gSystem->AccessPathName(folder+"/"+outputs[i]),"Existing TAG output will not be replaced");
    Require(!gSystem->AccessPathName(folder+"/context.json")||!gSystem->AccessPathName(folder+"/../context.json"),"TAG output must be provenance-prepared");
    TFile cache(cacheRoot,"READ");Require(!cache.IsZombie(),"Cannot read immutable TAG cache");TFile product(folder+"/tag_fits.root","CREATE");Require(!product.IsZombie(),"Cannot CREATE TAG output");
    std::ofstream json((folder+"/tag_fits.json").Data()),manifest((folder+"/slices.tsv").Data());Require(bool(json)&&bool(manifest),"Cannot create TAG JSON/TSV");json<<std::setprecision(17);manifest<<std::setprecision(17);
    json<<"{\"schemaVersion\":1,\"fitStage\":\"independent_lognormal_tag_shapes\",\"responseModel\":\"lognormal\",\"responseParameterization\":\"mode_and_local_curvature\",\"curveAdoptionAllowed\":false,\"gaussianVarianceUnfoldAllowed\":false,\"productionEnabled\":false,\"widthMeaning\":\"Component mean=lognormal mode; sigma=local curvature width, NOT arithmetic mean or standard deviation; projected response prototype, no Gaussian variance-unfold adoption\",\"likelihoodPolicy\":\"Only each TAG spectrum once, no ALL+TAG addition; same native R cell cannot be repeated within a TAG sample\",\"labelPolicy\":\"Nominal named components, no hard order or legacy-relative bounds; identities must be validated independently\",\"slices\":[";
    manifest<<"index\tsample\tRlo\tRhi\tfirstNativeRBin\tlastNativeRBin\tfitLo\tfitHi\tfirstFitDedxBin\tlastFitDedxBin\tfitCount\tstatus\tampConverged\tdataDeviance\tchi2\tndf\tatBoundary\n";
    std::map<std::string,std::set<int> > used;bool allNumerical=true;double totalD=0;
    for(unsigned i=0;i<slices.size();++i){const Slice& slice=slices[i];TH2* source=dynamic_cast<TH2*>(cache.Get(slice.sample.c_str()));Require(source!=0,"Missing TAG cache histogram "+slice.sample);
      const double tol=Number(cfg,"axisTolerance");int first=source->GetXaxis()->FindBin(slice.lo+tol),last=source->GetXaxis()->FindBin(slice.hi-tol);
      Require(first>=1&&last>=first&&last<=source->GetNbinsX()&&std::fabs(source->GetXaxis()->GetBinLowEdge(first)-slice.lo)<tol&&std::fabs(source->GetXaxis()->GetBinUpEdge(last)-slice.hi)<tol,"Requested R interval is not exact native bins");
      for(int b=first;b<=last;++b){Require(used[slice.sample].insert(b).second,"Repeated native R bin within TAG sample");Require(std::fabs(source->GetXaxis()->GetBinWidth(b)-Number(cfg,"nativeRigidityBinWidth"))<tol,"Unexpected native R width");}
      TH1D* native=source->ProjectionY(Form("raw_native_%03d",slice.index),first,last,"e");native->SetDirectory(0);
      Request request;request.nativeAxisMin=native->GetXaxis()->GetXmin();request.nativeAxisMax=native->GetXaxis()->GetXmax();request.backgroundModel=Text(cfg,"backgroundModel");request.freeBackgroundSlope=Integer(cfg,"freeBackgroundSlope")!=0;
      Seed background=NeedSeed(starts,slice.index,"BG");request.backgroundSlope=background.component.mean;request.backgroundWindowYield=background.component.startWindowYield;
      int yfirst=native->FindBin(slice.fitLo+tol),ylast=native->FindBin(slice.fitHi-tol);int rebin=Integer(cfg,"displayDedxRebinFactor");
      Require(rebin>0&&native->GetNbinsX()%rebin==0&&yfirst>=1&&ylast>=yfirst&&ylast<=native->GetNbinsX()&&std::fabs(native->GetXaxis()->GetBinLowEdge(yfirst)-slice.fitLo)<tol&&std::fabs(native->GetXaxis()->GetBinUpEdge(ylast)-slice.fitHi)<tol&&(yfirst-1)%rebin==0&&ylast%rebin==0,"Fit window must use exact native bins and display rebin boundaries");
      for(int b=0;b<=native->GetNbinsX()+1;++b){double count=native->GetBinContent(b),variance=std::pow(native->GetBinError(b),2);
        Require(Finite(count)&&count>=0&&Finite(variance)&&variance>=0&&std::fabs(count-variance)<=Number(cfg,"countVarianceRelativeTolerance")*std::max(1.,count),"Cached spectrum is not finite Poisson counts");
        if(b>=1&&b<=native->GetNbinsX())Require(std::fabs(native->GetBinWidth(b)-Number(cfg,"nativeDedxBinWidth"))<tol,"Unexpected native dEdx width");
        if(b>=yfirst&&b<=ylast)request.bins.push_back(Bin(native->GetXaxis()->GetBinLowEdge(b),native->GetXaxis()->GetBinUpEdge(b),count));}
      for(unsigned j=0;j<slice.names.size();++j)request.components.push_back(NeedSeed(starts,slice.index,slice.names[j]).component);
      const int normalize=Integer(cfg,"normalizeWarmStartToObserved");Require(normalize==0||normalize==1,"normalizeWarmStartToObserved requires 0/1");double seedNormalization=1;
      if(normalize&&TotalCount(request.bins)>0){double seedTotal=request.backgroundWindowYield;for(unsigned j=0;j<request.components.size();++j)if(request.components[j].active)seedTotal+=request.components[j].startWindowYield;
        Require(seedTotal>0,"Relative warm-start weights sum to zero");seedNormalization=TotalCount(request.bins)/seedTotal;request.backgroundWindowYield*=seedNormalization;
        for(unsigned j=0;j<request.components.size();++j)if(request.components[j].active)request.components[j].startWindowYield*=seedNormalization;}
      Result result=Fit(request,cfg);bool numerical=result.finite&&result.amplitudeConverged&&result.pointAccepted&&result.status==0&&(result.shapePoint.empty()||result.shapeCovarianceValid);
      allNumerical=allNumerical&&numerical;if(result.finite)totalD+=result.dataDeviance;
      product.cd();native->Write();TH1D* display=static_cast<TH1D*>(native->Clone(Form("raw_display_%03d",slice.index)));display->SetDirectory(0);double count=native->Integral(0,native->GetNbinsX()+1);display->Rebin(rebin);Require(display->Integral(0,display->GetNbinsX()+1)==count,"Display rebin changed count sum");display->Write();
      const unsigned n=request.components.size();std::vector<TH1D*> predicted;
      for(unsigned j=0;j<n+3;++j){std::string key=j==0?"total":(j==1?"smooth":(j==2?"fit_mask":request.components[j-3].name));TH1D* h=static_cast<TH1D*>(native->Clone(Form("prediction_%s_%03d",key.c_str(),slice.index)));h->Reset();h->SetDirectory(0);predicted.push_back(h);}
      for(unsigned b=0;b<request.bins.size();++b){predicted[2]->SetBinContent(yfirst+b,1);if(b>=result.totalCounts.size())continue;predicted[0]->SetBinContent(yfirst+b,result.totalCounts[b]);predicted[1]->SetBinContent(yfirst+b,result.smoothCounts[b]);for(unsigned j=0;j<n;++j)predicted[j+3]->SetBinContent(yfirst+b,result.componentCounts[b][j]);}
      for(unsigned j=0;j<predicted.size();++j){predicted[j]->Write();delete predicted[j];}
      if(i)json<<',';json<<"{\"index\":"<<slice.index<<",\"sample\":"<<Quote(slice.sample)<<",\"lo\":"<<slice.lo<<",\"hi\":"<<slice.hi<<",\"r\":"<<(slice.lo+slice.hi)/2<<",\"firstNativeRBin\":"<<first<<",\"lastNativeRBin\":"<<last<<",\"fitLo\":"<<slice.fitLo<<",\"fitHi\":"<<slice.fitHi<<",\"firstFitDedxBin\":"<<yfirst<<",\"lastFitDedxBin\":"<<ylast<<",\"fitCount\":"<<TotalCount(request.bins)<<",\"warmStartToObservedNormalized\":"<<(normalize?"true":"false")<<",\"warmStartNormalization\":"<<seedNormalization<<",\"backgroundModel\":"<<Quote(request.backgroundModel)<<",\"responseModel\":\"lognormal\",\"widthMeaning\":\"mode and local curvature width, not Gaussian mean/std; no variance-unfold adoption\",\"numericalConverged\":"<<(numerical?"true":"false")<<',';
      if(result.components.empty())json<<"\"status\":-2,\"emptyFitWindow\":true,\"components\":{}";else WriteResult(json,result,request);json<<'}';
      manifest<<slice.index<<'\t'<<slice.sample<<'\t'<<slice.lo<<'\t'<<slice.hi<<'\t'<<first<<'\t'<<last<<'\t'<<slice.fitLo<<'\t'<<slice.fitHi<<'\t'<<yfirst<<'\t'<<ylast<<'\t'<<TotalCount(request.bins)<<'\t'<<result.status<<'\t'<<result.amplitudeConverged<<'\t'<<result.dataDeviance<<'\t'<<result.pearson<<'\t'<<result.conditionalNdf<<'\t'<<result.atBoundary<<'\n';
      std::cout<<"LOGNORMAL_TAG_SLICE index="<<slice.index<<" sample="<<slice.sample<<" R="<<(slice.lo+slice.hi)/2<<" D="<<result.dataDeviance<<" status="<<result.status<<" AMP="<<result.amplitudeConverged<<std::endl;delete native;delete display;
    }
    json<<"],\"totalDataDeviance\":"<<totalD<<",\"numericalConverged\":"<<(allNumerical?"true":"false")<<",\"executionIsNotPhysicalValidation\":true}\n";json.close();manifest.close();product.Close();cache.Close();std::cout<<"LOGNORMAL_TAG_FIT_COMPLETE rows="<<slices.size()<<std::endl;
  }catch(const std::exception& e){std::cerr<<"FitLognormalTofTagShapes: "<<e.what()<<std::endl;gSystem->Exit(2);}
#endif
}
