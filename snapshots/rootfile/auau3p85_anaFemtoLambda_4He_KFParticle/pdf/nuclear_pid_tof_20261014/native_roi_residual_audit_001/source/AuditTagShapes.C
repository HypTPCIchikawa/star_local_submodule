// Diagnostics of retained TAG predictions only. No refit or gate changes.
void AuditTagShapes(const char*,const char*,const char*,const char*);
#ifndef __CINT__
#include "NuclearFitCommon.h"
#include "TFile.h"
#include "TH1.h"
#include <limits>
#include <set>
namespace NuclearTagResidual {
using namespace NuclearFit;
struct Row {int index;std::string sample,species;double lo,hi,mean,sigma,fitLo,fitHi;bool numerical;};
struct Metrics {
  double observed,predicted,primary,l1,pearson,deviance;int groups,activeGroups,observedGroups,zeroModelPositiveCountGroups,partialEdgeGroups;
  Metrics():observed(0),predicted(0),primary(0),l1(0),pearson(0),deviance(0),groups(0),activeGroups(0),observedGroups(0),zeroModelPositiveCountGroups(0),partialEdgeGroups(0){}
};
inline void Scalar(std::ostream& out,double x){if(Finite(x))out<<x;else out<<"null";}
inline void Cell(std::ostream& out,double x){if(Finite(x))out<<x;else out<<"nan";}
std::vector<double> Values(const std::string& text){std::istringstream in(text);std::vector<double> values;std::string token;
  while(std::getline(in,token,',')){std::istringstream field(token);double value;std::string extra;Require(bool(field>>value)&&!(field>>extra)&&Finite(value)&&value>0,"Invalid explicit residual diagnostic list");values.push_back(value);}Require(!values.empty(),"Empty residual diagnostic list");return values;}
std::vector<Row> Rows(const char* path){std::ifstream in(path);Require(bool(in),"Missing residual audit rows TSV");std::vector<Row> rows;std::string line;std::set<std::string> seen;
  while(std::getline(in,line)){if(line=="index\tsample\trLo\trHi\tspecies\tmean\tsigma\tfitLo\tfitHi\tnumericalValid")continue;Row r;int numerical;std::string extra;std::istringstream field(line);
    Require(bool(field>>r.index>>r.sample>>r.lo>>r.hi>>r.species>>r.mean>>r.sigma>>r.fitLo>>r.fitHi>>numerical)&&!(field>>extra),"Malformed 10-field diagnostic row");
    Require(r.index>=0&&(r.sample=="TAG_D"||r.sample=="TAG_T"||r.sample=="TAG_HE3")&&Finite(r.lo)&&Finite(r.hi)&&r.hi>r.lo&&Finite(r.mean)&&Finite(r.sigma)&&Finite(r.fitLo)&&Finite(r.fitHi)&&r.fitHi>r.fitLo&&(numerical==0||numerical==1),"Invalid diagnostic row");
    Require(r.species=="full_fit"||(r.sigma>0&&(r.species=="d"||r.species=="t"||r.species=="he3"||r.species=="he4")),"Invalid primary component");r.numerical=numerical;std::ostringstream key;key<<r.index<<':'<<r.species;Require(seen.insert(key.str()).second,"Duplicate diagnostic primary");rows.push_back(r);}
  Require(in.eof()&&!rows.empty(),"Unreadable/empty residual audit rows");return rows;}
Metrics Aggregate(const TH1& raw,const TH1& model,const TH1* primary,int first,int last,int factor){
  Metrics result;if(last<first)return result;int bin=first;
  while(bin<=last){int group=(bin-1)/factor,end=std::min(last,(group+1)*factor),begin=bin;double n=0,mu=0,s=0;
    for(;bin<=end;++bin){n+=raw.GetBinContent(bin);mu+=model.GetBinContent(bin);if(primary)s+=primary->GetBinContent(bin);}
    Require(Finite(n)&&n>=0&&Finite(mu)&&mu>=0&&Finite(s)&&s>=0,"Nonfinite/negative observed or retained prediction");++result.groups;result.partialEdgeGroups+=end-begin+1!=factor;result.observed+=n;result.predicted+=mu;result.primary+=s;result.l1+=std::fabs(n-mu);result.observedGroups+=n>0;
    if(mu>0){++result.activeGroups;result.pearson+=(n-mu)*(n-mu)/mu;if(n>0){double x=(mu-n)/n;result.deviance+=x>-1?2*n*LogRemainder(x):invalid;}else result.deviance+=2*mu;}
    else if(n>0){++result.zeroModelPositiveCountGroups;result.pearson=result.deviance=invalid;}
  }return result;
}
double Ratio(double a,double b){return b>0?a/b:std::numeric_limits<double>::quiet_NaN();}
}
#endif
void AuditTagShapes(const char* fitRoot,const char* output,const char* settings,const char* rowsTsv){
#ifndef __CINT__
  using namespace NuclearTagResidual;
  try{
    Settings cfg=ReadSettings(settings);std::vector<Row> rows=Rows(rowsTsv);Vec sigmas=Values(Text(cfg,"roiNSigma")),widths=Values(Text(cfg,"aggregateWidthsKeV"));const double tolerance=Number(cfg,"axisTolerance"),nativeWidth=Number(cfg,"nativeDedxBinWidth");Require(tolerance>0&&nativeWidth>0,"Invalid numerical diagnostic geometry");
    const TString folder(output);Require(gSystem->AccessPathName(folder+"/residuals.tsv")&&gSystem->AccessPathName(folder+"/residual_audit.json"),"Diagnostic products already exist");Require(!gSystem->AccessPathName(folder+"/context.json"),"Diagnostic provenance must be prepared");
    TFile source(fitRoot,"READ");Require(!source.IsZombie()&&!source.TestBit(TFile::kRecovered),"Cannot read retained TAG prediction ROOT");std::ofstream csv((folder+"/residuals.tsv").Data()),json((folder+"/residual_audit.json").Data());Require(bool(csv)&&bool(json),"Cannot CREATE residual diagnostics");csv<<std::setprecision(17);json<<std::setprecision(17);
    csv<<"index\tsample\trLo\trHi\tspecies\tnSigma\taggregateKeV\tmean\tprojectionSigma\tphysicalGaussianCoverage\tfitWindowGaussianCoverage\trequestedLo\trequestedHi\tactualLo\tactualHi\tfirstNativeDedxBin\tlastNativeDedxBin\trequestedGaussianRoiProbability\tactualGaussianRoiProbability\tobserved\tpredicted\tprimaryModelCount\tl1\tl1OverObserved\tpearson\tdeviance\tgroups\tactiveGroups\tobservedGroups\tpearsonPerActiveGroup\tsqrtPearsonOverObserved\tzeroModelPositiveCountGroups\tpartialEdgeGroups\tnumericalValid\n";
    json<<"{\"schemaVersion\":1,\"noFit\":true,\"productionEnabled\":false,\"widthMeaning\":\"Retained free TAG native-R projection Gaussian width, NOT an intrinsic shared width\",\"roiPolicy\":\"Per-model retained mean +/- declared sigma; clipped to physical axis and actual fitted window, outward snapped to whole native dEdx bins\",\"aggregationPolicy\":\"Sum the same observed and retained expected native cells, detector-axis-aligned groups, partial edge groups retained; NO refit\",\"chi2Policy\":\"Pearson / number of aggregate groups with positive expected count; NOT fitted chi2/ndf\",\"identityValidated\":false,\"records\":[";unsigned records=0;
    for(unsigned r=0;r<rows.size();++r){const Row& row=rows[r];TH1* raw=dynamic_cast<TH1*>(source.Get(Form("raw_native_%03d",row.index)));TH1* model=dynamic_cast<TH1*>(source.Get(Form("prediction_total_%03d",row.index)));TH1* mask=dynamic_cast<TH1*>(source.Get(Form("prediction_fit_mask_%03d",row.index)));TH1* primary=row.species=="full_fit"?0:dynamic_cast<TH1*>(source.Get(Form("prediction_%s_%03d",row.species.c_str(),row.index)));
      Require(raw&&model&&mask&&(row.species=="full_fit"||primary)&&raw->GetDimension()==1&&raw->GetNbinsX()==model->GetNbinsX(),"Missing/invalid retained native or component prediction");const double physicalLo=raw->GetXaxis()->GetXmin(),physicalHi=raw->GetXaxis()->GetXmax();
      for(int b=1;b<=raw->GetNbinsX();++b)Require(std::fabs(raw->GetBinWidth(b)-nativeWidth)<tolerance&&raw->GetBinLowEdge(b)==model->GetBinLowEdge(b),"Unexpected native axis");
      std::vector<double> roi=row.species=="full_fit"?Vec(1,0):sigmas;
      for(unsigned k=0;k<roi.size();++k){double requestedLo=row.species=="full_fit"?row.fitLo:row.mean-roi[k]*row.sigma,requestedHi=row.species=="full_fit"?row.fitHi:row.mean+roi[k]*row.sigma;
        double lo=std::max(std::max(requestedLo,physicalLo),row.fitLo),hi=std::min(std::min(requestedHi,physicalHi),row.fitHi);int first=1,last=0;
        if(hi>lo){first=raw->FindBin(lo+tolerance);last=raw->FindBin(hi-tolerance);first=std::max(1,first);last=std::min(raw->GetNbinsX(),last);}
        for(int b=first;b<=last;++b)Require(mask->GetBinContent(b)==1,"ROI extends beyond saved fit mask");double actualLo=last>=first?raw->GetBinLowEdge(first):lo,actualHi=last>=first?raw->GetXaxis()->GetBinUpEdge(last):lo;
        double physicalProbability=row.species=="full_fit"?1:GaussianIntegral(physicalLo,physicalHi,row.mean,row.sigma),fitProbability=row.species=="full_fit"?1:GaussianIntegral(row.fitLo,row.fitHi,row.mean,row.sigma),requestedProbability=row.species=="full_fit"?1:GaussianIntegral(requestedLo,requestedHi,row.mean,row.sigma),actualProbability=row.species=="full_fit"?1:GaussianIntegral(actualLo,actualHi,row.mean,row.sigma);
        for(unsigned w=0;w<widths.size();++w){double factorValue=widths[w]/nativeWidth;int factor=int(std::floor(factorValue+.5));Require(factor>0&&std::fabs(factorValue-factor)<tolerance,"Aggregate width is not an integer native-bin grouping");Metrics m=Aggregate(*raw,*model,primary,first,last,factor);
          double pearson=m.zeroModelPositiveCountGroups?std::numeric_limits<double>::quiet_NaN():m.pearson,deviance=m.zeroModelPositiveCountGroups?std::numeric_limits<double>::quiet_NaN():m.deviance;
          if(records++)json<<',';json<<"{\"index\":"<<row.index<<",\"sample\":"<<Quote(row.sample)<<",\"lo\":"<<row.lo<<",\"hi\":"<<row.hi<<",\"species\":"<<Quote(row.species)<<",\"nSigma\":"<<roi[k]<<",\"aggregateKeV\":"<<widths[w]<<",\"mean\":"<<row.mean<<",\"projectionSigma\":"<<row.sigma<<",\"physicalGaussianCoverage\":"<<physicalProbability<<",\"fitWindowGaussianCoverage\":"<<fitProbability<<",\"requestedLo\":"<<requestedLo<<",\"requestedHi\":"<<requestedHi<<",\"actualLo\":"<<actualLo<<",\"actualHi\":"<<actualHi<<",\"firstNativeDedxBin\":"<<first<<",\"lastNativeDedxBin\":"<<last<<",\"requestedGaussianRoiProbability\":"<<requestedProbability<<",\"actualGaussianRoiProbability\":"<<actualProbability<<",\"observed\":"<<m.observed<<",\"predicted\":"<<m.predicted<<",\"primaryModelCount\":"<<m.primary<<",\"l1\":"<<m.l1<<",\"l1OverObserved\":";Scalar(json,Ratio(m.l1,m.observed));json<<",\"pearson\":";Scalar(json,pearson);json<<",\"deviance\":";Scalar(json,deviance);json<<",\"groups\":"<<m.groups<<",\"activeGroups\":"<<m.activeGroups<<",\"observedGroups\":"<<m.observedGroups<<",\"pearsonPerActiveGroup\":";Scalar(json,Ratio(pearson,m.activeGroups));json<<",\"sqrtPearsonOverObserved\":";Scalar(json,std::sqrt(Ratio(pearson,m.observed)));json<<",\"zeroModelPositiveCountGroups\":"<<m.zeroModelPositiveCountGroups<<",\"partialEdgeGroups\":"<<m.partialEdgeGroups<<",\"numericalValid\":"<<(row.numerical?"true":"false")<<'}';
          csv<<row.index<<'\t'<<row.sample<<'\t'<<row.lo<<'\t'<<row.hi<<'\t'<<row.species<<'\t'<<roi[k]<<'\t'<<widths[w]<<'\t'<<row.mean<<'\t'<<row.sigma<<'\t'<<physicalProbability<<'\t'<<fitProbability<<'\t'<<requestedLo<<'\t'<<requestedHi<<'\t'<<actualLo<<'\t'<<actualHi<<'\t'<<first<<'\t'<<last<<'\t'<<requestedProbability<<'\t'<<actualProbability<<'\t'<<m.observed<<'\t'<<m.predicted<<'\t'<<m.primary<<'\t'<<m.l1<<'\t';Cell(csv,Ratio(m.l1,m.observed));csv<<'\t';Cell(csv,pearson);csv<<'\t';Cell(csv,deviance);csv<<'\t'<<m.groups<<'\t'<<m.activeGroups<<'\t'<<m.observedGroups<<'\t';Cell(csv,Ratio(pearson,m.activeGroups));csv<<'\t';Cell(csv,std::sqrt(Ratio(pearson,m.observed)));csv<<'\t'<<m.zeroModelPositiveCountGroups<<'\t'<<m.partialEdgeGroups<<'\t'<<row.numerical<<'\n';
        }
      }
    }
    json<<"],\"recordCount\":"<<records<<",\"residualMetricsAreSupportDiagnosticsNotNewFitNdf\":true}\n";json.close();csv.close();source.Close();std::cout<<"NUCLEAR_TAG_RESIDUAL_AUDIT_COMPLETE records="<<records<<std::endl;
  }catch(const std::exception& e){std::cerr<<"AuditTagShapes: "<<e.what()<<std::endl;gSystem->Exit(2);}
#endif
}
