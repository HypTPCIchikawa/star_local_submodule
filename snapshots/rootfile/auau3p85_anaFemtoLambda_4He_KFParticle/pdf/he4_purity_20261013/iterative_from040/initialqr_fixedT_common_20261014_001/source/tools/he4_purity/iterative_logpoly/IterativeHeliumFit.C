// Study-only continuation from the frozen from040 independent fits.
// Local widths are R-projection widths; COMMON widths are intrinsic GL kernels.
void IterativeHeliumFit(const char* cacheRoot,const char* output,const char* settings,
                        const char* seeds,const char* targets);
#ifndef __CINT__
#include "../global_logpoly/GlobalHeliumFit.C"
#include "TH1I.h"
#include <limits>
namespace HeliumIterative {
using GlobalHelium::Vec;using GlobalHelium::Matrix;using GlobalHelium::Settings;
using GlobalHelium::Value;using GlobalHelium::Text;using GlobalHelium::Poly;
using GlobalHelium::Extreme;using GlobalHelium::Gaussian;using GlobalHelium::Quote;
using GlobalHelium::ReadSettings;using GlobalHelium::WriteVec;
const double invalid=1.e20;
struct Seed {
  int index,oldChosen;bool oldValid,oldQuality,tPresent;std::string bg;
  double lo,hi,anchorMean,anchorSigma,mu[3],sigma[3],yield[3],bgYield,bgSlope;
};
struct Weight {
  double info[2],meanScale[2],widthScale[2];
  Weight(){for(int j=0;j<2;++j){info[j]=0;meanScale[j]=widthScale[j]=1;}}
};
struct Target {Vec coefficients;std::map<int,Weight> weights;bool present;Target():present(false){}};
struct Shape {
  double mu[3],sigma[3],slope;
  Shape(){for(int j=0;j<3;++j){mu[j]=0;sigma[j]=1;}slope=0;}
};
struct Data {
  Seed seed;GlobalHelium::Slice raw;TH1D *native,*h;Shape shape;Weight weight;
  bool included[3],shapeFree[3],accepted,ampConverged,finite,boundary,labelSwapped;
  int status,covStatus,ampIterations,freeCount;double edm,ampKkt,dataD,penalty,pearson;
  Vec a,labelPermutation,nodes,nodeWeights;Matrix covariance;
  Data():native(0),h(0),accepted(false),ampConverged(false),finite(false),boundary(false),
    labelSwapped(false),status(-1),covStatus(-1),ampIterations(0),freeCount(0),
    edm(invalid),ampKkt(invalid),dataD(invalid),penalty(0),pearson(0){
      for(int j=0;j<3;++j){included[j]=true;shapeFree[j]=true;labelPermutation.push_back(j);}
  }
};
Target ReadTarget(const char* path){
  std::ifstream f(path);if(!f)throw std::runtime_error("missing explicit targets TSV");
  Target target;std::string type;if(!(f>>type))throw std::runtime_error("empty target file");
  if(type=="none"){std::string extra;if(f>>extra)throw std::runtime_error("unexpected fields after none");return target;}
  if(type!="coefficients")throw std::runtime_error("targets must start coefficients or none");
  target.present=true;target.coefficients.resize(14);
  for(int j=0;j<14;++j)if(!(f>>target.coefficients[j])||!TMath::Finite(target.coefficients[j]))throw std::runtime_error("invalid target coefficients");
  while(f>>type){if(type!="slice")throw std::runtime_error("invalid target row type");int index;Weight w;
    if(!(f>>index>>w.info[0]>>w.meanScale[0]>>w.widthScale[0]>>w.info[1]>>w.meanScale[1]>>w.widthScale[1]))throw std::runtime_error("invalid target weight row");
    for(int j=0;j<2;++j)if(!TMath::Finite(w.info[j])||!TMath::Finite(w.meanScale[j])||!TMath::Finite(w.widthScale[j])||w.info[j]<0||w.meanScale[j]<=0||w.widthScale[j]<=0)throw std::runtime_error("nonfinite or negative target information/scale");
    if(!target.weights.insert(std::make_pair(index,w)).second)throw std::runtime_error("duplicate target slice");
  }return target;
}
std::vector<Seed> ReadSeeds(const char* path){
  std::ifstream f(path);if(!f)throw std::runtime_error("missing explicit from040 seeds");
  std::vector<Seed> rows;std::set<int> seen;std::string line;
  while(std::getline(f,line)){if(line.empty())throw std::runtime_error("empty seed row");
    Seed s;int valid,quality,t;std::istringstream r(line);std::string extra;
    if(!(r>>s.index>>s.lo>>s.hi>>s.anchorMean>>s.anchorSigma>>s.oldChosen>>valid>>quality>>s.bg
       >>t>>s.mu[0]>>s.sigma[0]>>s.yield[0]>>s.mu[1]>>s.sigma[1]>>s.yield[1]
       >>s.mu[2]>>s.sigma[2]>>s.yield[2]>>s.bgYield>>s.bgSlope)||(r>>extra))throw std::runtime_error("invalid 21-field from040 seed");
    if(!seen.insert(s.index).second||s.index<0||!(s.lo>0&&s.hi>s.lo&&s.anchorSigma>0)||
       (valid!=0&&valid!=1)||(quality!=0&&quality!=1)||(t!=0&&t!=1)||
       (s.bg!="linear"&&s.bg!="exponential"))throw std::runtime_error("invalid seed identifiers or flags");
    s.oldValid=valid;s.oldQuality=quality;s.tPresent=t;
    for(int j=0;j<3;++j)if(!TMath::Finite(s.mu[j])||!TMath::Finite(s.sigma[j])||!TMath::Finite(s.yield[j])||s.yield[j]<0||((j>0||s.tPresent)&&s.sigma[j]<=0))throw std::runtime_error("nonfinite from040 Gaussian seed");
    if(!TMath::Finite(s.bgYield)||s.bgYield<0||!TMath::Finite(s.bgSlope))throw std::runtime_error("invalid background seed");
    rows.push_back(s);
  }if(rows.empty())throw std::runtime_error("no from040 seeds");return rows;
}
double SharedMean(const Vec& g,int species,double r){return Poly(&g[7*species],4,std::log10(r));}
double SharedSigma(const Vec& g,int species,double r){return Poly(&g[7*species+4],3,std::log10(r));}
Vec ToNormal(const Vec& g,double mid,double half){
  Vec n(14);for(int sp=0;sp<2;++sp){int o=sp*7;const double* c=&g[o];
    n[o]=Poly(c,4,mid);n[o+1]=half*(c[1]+2*c[2]*mid+3*c[3]*mid*mid);
    n[o+2]=half*half*(c[2]+3*c[3]*mid);n[o+3]=half*half*half*c[3];
    n[o+4]=Poly(c+4,3,mid);n[o+5]=half*(c[5]+2*c[6]*mid);n[o+6]=half*half*c[6];
  }return n;
}
Vec FromNormal(const Vec& n,double mid,double half){
  Vec g(14);for(int sp=0;sp<2;++sp){int o=sp*7;double c3=n[o+3]/std::pow(half,3),b2=n[o+2]/(half*half),b1=n[o+1]/half;
    g[o+3]=c3;g[o+2]=b2-3*mid*c3;g[o+1]=b1-2*mid*b2+3*mid*mid*c3;
    g[o]=n[o]-mid*b1+mid*mid*b2-mid*mid*mid*c3;
    g[o+6]=n[o+6]/(half*half);g[o+5]=n[o+5]/half-2*mid*g[o+6];
    g[o+4]=n[o+4]-mid*n[o+5]/half+mid*mid*g[o+6];
  }return g;
}
double WidthShapeMinimum(double b,double a){const double c[]={0,b,a};return Extreme(c,2,-1,1,true);}
Vec ToOptim(const Vec& g,double mid,double half,double epsilon){
  Vec p=ToNormal(g,mid,half);for(int sp=0;sp<2;++sp){int o=sp*7;double positive=p[o+4]+WidthShapeMinimum(p[o+5],p[o+6])-epsilon;
    if(!(positive>0))throw std::runtime_error("shared seed width must exceed positive numerical floor");p[o+4]=std::log(positive);
  }return p;
}
Vec FromOptim(const Vec& p,double mid,double half,double epsilon){
  Vec n=p;for(int sp=0;sp<2;++sp){int o=sp*7;n[o+4]=epsilon+std::exp(p[o+4])-WidthShapeMinimum(p[o+5],p[o+6]);}
  return FromNormal(n,mid,half);
}
struct Profile {Vec a;double d,kkt;int iterations;bool converged;Profile():d(invalid),kkt(invalid),iterations(0),converged(false){}};
Profile WindowProfile(const GlobalHelium::Slice& data,const Matrix& t,const Vec& initial,const Settings& cfg){
  const int n=t[0].size();Profile p;p.a=initial.size()==unsigned(n)?initial:Vec(n,1./n);
  for(int j=0;j<n;++j)p.a[j]=std::max(0.,p.a[j]);
  double previous=GlobalHelium::Deviance(data,t,p.a);
  for(int iteration=0;iteration<int(Value(cfg,"amplitudeMaxIterations"));++iteration){
    Vec grad(n,0);Matrix h(n,Vec(n,0));
    for(unsigned b=0;b<t.size();++b){double v=0;for(int j=0;j<n;++j)v+=p.a[j]*t[b][j];v=std::max(1.e-20,v);
      double observed=data.bins[b].n/data.count;
      for(int j=0;j<n;++j){grad[j]+=t[b][j]*(1-observed/v);for(int k=0;k<n;++k)h[j][k]+=observed*t[b][j]*t[b][k]/(v*v);}
    }
    std::vector<int> active;double kkt=0;
    for(int j=0;j<n;++j)if(!(p.a[j]<1.e-12&&grad[j]>0)){active.push_back(j);kkt=std::max(kkt,std::fabs(grad[j]));}
    p.iterations=iteration+1;if(kkt<Value(cfg,"amplitudeGradientTolerance"))break;
    Matrix ha(active.size(),Vec(active.size()));Vec rhs(active.size());
    for(unsigned j=0;j<active.size();++j){rhs[j]=-grad[active[j]];for(unsigned k=0;k<active.size();++k)ha[j][k]=h[active[j]][active[k]];ha[j][j]+=1.e-10;}
    Vec direction;bool solved=GlobalHelium::Solve(ha,rhs,direction),accepted=false;
    for(int strategy=0;strategy<2&&!accepted;++strategy){
      if(strategy==1||!solved){direction.resize(active.size());for(unsigned j=0;j<active.size();++j)direction[j]=-grad[active[j]]/std::max(1.e-8,h[active[j]][active[j]]);}
      double step=1;for(unsigned j=0;j<active.size();++j)if(direction[j]<0)step=std::min(step,-p.a[active[j]]/direction[j]);
      for(int line=0;line<40&&step>1.e-16;++line){Vec next=p.a;for(unsigned j=0;j<active.size();++j)next[active[j]]=std::max(0.,p.a[active[j]]+step*direction[j]);double d=GlobalHelium::Deviance(data,t,next);
        if(d<previous){p.a=next;previous=d;accepted=true;break;}step*=.5;
      }
    }if(!accepted)break;
  }
  Vec gradient(n,0);for(unsigned b=0;b<t.size();++b){double v=0;for(int j=0;j<n;++j)v+=p.a[j]*t[b][j];v=std::max(1.e-20,v);for(int j=0;j<n;++j)gradient[j]+=t[b][j]*(1-data.bins[b].n/(data.count*v));}
  p.kkt=0;for(int j=0;j<n;++j)if(!(p.a[j]<1.e-12&&gradient[j]>0))p.kkt=std::max(p.kkt,std::fabs(gradient[j]));
  p.converged=p.kkt<Value(cfg,"amplitudeGradientTolerance");p.d=GlobalHelium::Deviance(data,t,p.a);return p;
}
class Study {
 public:
  Settings cfg;Target target;std::string phase;std::vector<Data> rows;double mid,half;
  Study(const Settings& c,const Target& t):cfg(c),target(t),phase(Text(c,"phase")){
    double lo=std::log10(Value(c,"rigidityMin")),hi=std::log10(Value(c,"rigidityMax"));mid=(lo+hi)/2;half=(hi-lo)/2;
    if(!(half>0)||(phase!="free"&&phase!="regularized"&&phase!="common"))throw std::runtime_error("invalid phase/domain");
    if((phase=="free"&&(target.present||Value(c,"penaltyStrength")!=0))||(phase!="free"&&!target.present)||(phase=="common"&&Value(c,"penaltyStrength")!=0))throw std::runtime_error("inconsistent target/penalty phase");
  }
  bool SharedFeasible(const Vec& g,double lowAxis,double highAxis)const{
    if(g.size()!=14)return false;for(unsigned j=0;j<g.size();++j)if(!TMath::Finite(g[j]))return false;
    double lo=mid-half,hi=mid+half;
    for(int sp=0;sp<2;++sp){int o=sp*7;
      if(Extreme(&g[o],3,lo,hi,true)<lowAxis-Value(cfg,"meanAxisPadding")||Extreme(&g[o],3,lo,hi,false)>highAxis+Value(cfg,"meanAxisPadding"))return false;
      if(Extreme(&g[o+4],2,lo,hi,true)<Value(cfg,"minimumSigma")||Extreme(&g[o+4],2,lo,hi,false)>Value(cfg,"maximumSigma"))return false;
    }return true;
  }
  bool ShapeFeasible(const Data& d,const Shape& shape)const{
    for(int j=0;j<3;++j)if(d.included[j]&&(!TMath::Finite(shape.mu[j])||!TMath::Finite(shape.sigma[j])||
       shape.mu[j]<d.h->GetXaxis()->GetXmin()-Value(cfg,"meanAxisPadding")||shape.mu[j]>d.h->GetXaxis()->GetXmax()+Value(cfg,"meanAxisPadding")||
       shape.sigma[j]<Value(cfg,"minimumSigma")||shape.sigma[j]>Value(cfg,"maximumSigma")))return false;
    return TMath::Finite(shape.slope)&&std::fabs(shape.slope)<=Value(cfg,"exponentialSlopeLimit");
  }
  Matrix Templates(const Data& d,const Shape& shape,const Vec& shared,Vec* acceptance=0)const{
    int n=d.seed.bg=="linear"?5:4;Matrix t(d.raw.bins.size(),Vec(n,0));Vec accepted(3,0);
    for(int sp=0;sp<3;++sp)if(d.included[sp]){
      if(phase=="common"&&sp>0){for(unsigned q=0;q<d.nodes.size();++q)accepted[sp]+=d.nodeWeights[q]*Gaussian(d.raw.fitLo,d.raw.fitHi,SharedMean(shared,sp-1,d.nodes[q]),SharedSigma(shared,sp-1,d.nodes[q]));}
      else accepted[sp]=Gaussian(d.raw.fitLo,d.raw.fitHi,shape.mu[sp],shape.sigma[sp]);
      if(!(accepted[sp]>Value(cfg,"minimumGaussianAcceptance")))return Matrix();
      for(unsigned b=0;b<t.size();++b){const GlobalHelium::Bin& bin=d.raw.bins[b];
        if(phase=="common"&&sp>0)for(unsigned q=0;q<d.nodes.size();++q)t[b][sp]+=d.nodeWeights[q]*Gaussian(bin.lo,bin.hi,SharedMean(shared,sp-1,d.nodes[q]),SharedSigma(shared,sp-1,d.nodes[q]));
        else t[b][sp]=Gaussian(bin.lo,bin.hi,shape.mu[sp],shape.sigma[sp]);
        t[b][sp]/=accepted[sp];
      }
    }
    double width=d.raw.fitHi-d.raw.fitLo;
    for(unsigned b=0;b<t.size();++b){const GlobalHelium::Bin& bin=d.raw.bins[b];
      if(d.seed.bg=="linear"){double center=(bin.lo+bin.hi)/2;t[b][3]=2*(bin.hi-bin.lo)*(d.raw.fitHi-center)/(width*width);t[b][4]=2*(bin.hi-bin.lo)*(center-d.raw.fitLo)/(width*width);}
      else{double k=shape.slope;t[b][3]=std::fabs(k)<1.e-8?(bin.hi-bin.lo)/width:std::exp(k*(bin.lo-d.raw.fitLo)/width)*std::expm1(k*(bin.hi-bin.lo)/width)/std::expm1(k);}
    }
    for(int col=0;col<n;++col){double sum=0;for(unsigned b=0;b<t.size();++b)sum+=t[b][col];if(sum>0)for(unsigned b=0;b<t.size();++b)t[b][col]/=sum;}
    if(acceptance)*acceptance=accepted;return t;
  }
  double Penalty(const Data& d,const Shape& shape)const{
    if(phase!="regularized")return 0;double loss=0;
    for(int sp=0;sp<2;++sp){double r=d.raw.r,mu=SharedMean(target.coefficients,sp,r),sigma=SharedSigma(target.coefficients,sp,r);
      double dm=(shape.mu[sp+1]-mu)/d.weight.meanScale[sp],ds=(shape.sigma[sp+1]-sigma)/d.weight.widthScale[sp];
      loss+=d.weight.info[sp]*(dm*dm+ds*ds);
    }return Value(cfg,"penaltyStrength")*loss;
  }
  double LocalLoss(const Data& d,const Shape& shape,const Vec& shared)const{
    if(!ShapeFeasible(d,shape))return invalid;Matrix t=Templates(d,shape,shared);if(t.empty())return invalid;
    Profile p=WindowProfile(d.raw,t,d.a,cfg);return p.converged?p.d+Penalty(d,shape):invalid;
  }
  double CommonLoss(const Vec& g)const{
    if(rows.empty()||!SharedFeasible(g,rows[0].h->GetXaxis()->GetXmin(),rows[0].h->GetXaxis()->GetXmax()))return invalid;
    double sum=0;for(unsigned j=0;j<rows.size();++j){double value=LocalLoss(rows[j],rows[j].shape,g);if(value>=1.e19)return invalid;sum+=value;}return sum;
  }
  void Refresh(Data& d,const Vec& shared){
    if(phase=="common")for(int sp=1;sp<3;++sp){d.shape.mu[sp]=SharedMean(shared,sp-1,d.raw.r);d.shape.sigma[sp]=SharedSigma(shared,sp-1,d.raw.r);}
    Vec acceptance;Matrix t=Templates(d,d.shape,shared,&acceptance);
    if(t.empty()){d.finite=false;d.ampConverged=false;d.dataD=invalid;return;}
    Profile p=WindowProfile(d.raw,t,d.a,cfg);d.a=p.a;d.ampConverged=p.converged;d.ampKkt=p.kkt;d.ampIterations=p.iterations;d.dataD=p.d;d.penalty=Penalty(d,d.shape);d.finite=TMath::Finite(p.d)&&p.d<1.e19;
    d.pearson=0;for(unsigned b=0;b<t.size();++b){double expected=0;for(unsigned k=0;k<d.a.size();++k)expected+=d.raw.count*d.a[k]*t[b][k];double residual=d.raw.bins[b].n-expected;d.pearson+=residual*residual/std::max(1.e-20,expected);}
  }
};
struct LocalObjective {
  Study* study;Data* data;const Vec* shared;std::vector<int> free;
  Shape Decode(const double* p)const{Shape s=data->shape;for(unsigned j=0;j<free.size();++j){int index=free[j];if(index==6)s.slope=p[j];else if(index%2)s.sigma[index/2]=std::exp(p[j]);else s.mu[index/2]=p[j];}return s;}
  double operator()(const double* p)const{return study->LocalLoss(*data,Decode(p),*shared);}
};
void FitLocal(Study& study,Data& d,const Vec& shared){
  LocalObjective objective;objective.study=&study;objective.data=&d;objective.shared=&shared;
  for(int sp=0;sp<3;++sp)if(d.included[sp]&&d.shapeFree[sp]){objective.free.push_back(2*sp);objective.free.push_back(2*sp+1);}
  if(d.seed.bg=="exponential")objective.free.push_back(6);
  d.freeCount=objective.free.size();d.boundary=false;d.covariance.clear();
  if(!objective.free.empty()){
    ROOT::Math::Minimizer* m=ROOT::Math::Factory::CreateMinimizer("Minuit2","Migrad");if(!m)throw std::runtime_error("local Minuit2 unavailable");
    ROOT::Math::Functor f(objective,objective.free.size());m->SetFunction(f);m->SetPrintLevel(-1);
    m->SetMaxFunctionCalls(int(Value(study.cfg,"localMaxCalls")));m->SetTolerance(Value(study.cfg,"localTolerance"));
    Vec initial;for(unsigned j=0;j<objective.free.size();++j){int index=objective.free[j];double value,low,high,step;
      if(index==6){value=d.shape.slope;low=-Value(study.cfg,"exponentialSlopeLimit");high=-low;step=Value(study.cfg,"backgroundSlopeStep");}
      else if(index%2){value=std::log(d.shape.sigma[index/2]);low=std::log(Value(study.cfg,"minimumSigma"));high=std::log(Value(study.cfg,"maximumSigma"));step=Value(study.cfg,"sigmaLogStep");}
      else{value=d.shape.mu[index/2];low=d.h->GetXaxis()->GetXmin()-Value(study.cfg,"meanAxisPadding");high=d.h->GetXaxis()->GetXmax()+Value(study.cfg,"meanAxisPadding");step=Value(study.cfg,"meanStep");}
      if(value<low||value>high||!TMath::Finite(value))throw std::runtime_error("from040/continuation seed outside declared numerical support");
      initial.push_back(value);m->SetLimitedVariable(j,Form("local_%d",index),value,step,low,high);
    }
    double before=objective(&initial[0]);bool ok=m->Minimize();Shape candidate=objective.Decode(m->X());double proposed=study.LocalLoss(d,candidate,shared);
    d.accepted=proposed<1.e19&&proposed<=before;if(d.accepted)d.shape=candidate;
    d.status=m->Status();d.edm=TMath::Finite(m->Edm())?m->Edm():invalid;d.covStatus=m->CovMatrixStatus();
    if(d.accepted&&d.covStatus>0){d.covariance.assign(objective.free.size(),Vec(objective.free.size()));for(unsigned j=0;j<objective.free.size();++j)for(unsigned k=0;k<objective.free.size();++k)d.covariance[j][k]=m->CovMatrix(j,k);}
    if(!ok&&d.status==0)d.status=1;
    for(unsigned j=0;j<objective.free.size();++j){int index=objective.free[j];double value=m->X()[j],lo,hi;
      if(index==6){lo=-Value(study.cfg,"exponentialSlopeLimit");hi=-lo;}
      else if(index%2){lo=std::log(Value(study.cfg,"minimumSigma"));hi=std::log(Value(study.cfg,"maximumSigma"));}
      else{lo=d.h->GetXaxis()->GetXmin()-Value(study.cfg,"meanAxisPadding");hi=d.h->GetXaxis()->GetXmax()+Value(study.cfg,"meanAxisPadding");}
      if(std::min(std::fabs(value-lo),std::fabs(hi-value))<Value(study.cfg,"axisTolerance"))d.boundary=true;
    }delete m;
  }else{d.status=0;d.covStatus=-1;d.edm=0;d.accepted=study.LocalLoss(d,d.shape,shared)<1.e19;}
  study.Refresh(d,shared);
}
// Permutations preserve data likelihood, not nuclear identity. At free phase
// only, canonicalize active FREE Gaussian labels by ascending mean.
void CanonicalizeFree(Data& d){
  std::vector<int> active;for(int sp=0;sp<3;++sp)if(d.included[sp]&&d.shapeFree[sp])active.push_back(sp);
  std::vector<int> sorted=active;for(unsigned j=0;j<sorted.size();++j)for(unsigned k=j+1;k<sorted.size();++k)if(d.shape.mu[sorted[k]]<d.shape.mu[sorted[j]])std::swap(sorted[k],sorted[j]);
  Shape old=d.shape;Vec amplitudes=d.a,permutation=d.labelPermutation;
  for(unsigned j=0;j<active.size();++j){int to=active[j],from=sorted[j];if(to!=from)d.labelSwapped=true;d.shape.mu[to]=old.mu[from];d.shape.sigma[to]=old.sigma[from];d.a[to]=amplitudes[from];d.labelPermutation[to]=permutation[from];}
  if(d.labelSwapped){d.covariance.clear();d.covStatus=-1;}
}
struct CommonObjective {
  Study* study;
  double operator()(const double* p)const{return study->CommonLoss(FromOptim(Vec(p,p+14),study->mid,study->half,Value(study->cfg,"minimumSigma")));}
};
bool AllLocalValid(const Study& s){for(unsigned j=0;j<s.rows.size();++j)if(s.rows[j].status!=0||!s.rows[j].accepted||!s.rows[j].ampConverged||!s.rows[j].finite)return false;return true;}
struct CommonResult {
  Vec g,trace;int status,covStatus;double edm,coefficientDrift;bool accepted,outerConverged;
  bool finalLocalStationary,finalGlobalStationary,stationarityEvaluable;
  double localImprovementPerCount,gradientAbsolute,gradientPerCount;
  int centralDirections,oneSidedDirections,blockedDirections,profileFailures;
  Matrix covariance;Vec covariancePoint;
  CommonResult():status(-1),covStatus(-1),edm(invalid),coefficientDrift(invalid),accepted(false),outerConverged(false),
    finalLocalStationary(false),finalGlobalStationary(false),stationarityEvaluable(false),
    localImprovementPerCount(invalid),gradientAbsolute(invalid),gradientPerCount(invalid),
    centralDirections(0),oneSidedDirections(0),blockedDirections(0),profileFailures(0){}
};
CommonResult FitCommon(Study& study){
  CommonResult result;result.g=study.target.coefficients;
  if(!study.SharedFeasible(result.g,study.rows[0].h->GetXaxis()->GetXmin(),study.rows[0].h->GetXaxis()->GetXmax()))throw std::runtime_error("shared continuation seed outside mathematical detector support");
  Vec p=ToOptim(result.g,study.mid,study.half,Value(study.cfg,"minimumSigma"));double previous=invalid;
  for(int outer=0;outer<int(Value(study.cfg,"maximumCommonOuterIterations"));++outer){
    for(unsigned j=0;j<study.rows.size();++j)FitLocal(study,study.rows[j],result.g);
    CommonObjective objective;objective.study=&study;double before=objective(&p[0]);Vec old=p;
    ROOT::Math::Minimizer* m=ROOT::Math::Factory::CreateMinimizer("Minuit2","Migrad");if(!m)throw std::runtime_error("common Minuit2 unavailable");
    ROOT::Math::Functor f(objective,14);m->SetFunction(f);m->SetPrintLevel(-1);m->SetMaxFunctionCalls(int(Value(study.cfg,"commonMaxCalls")));m->SetTolerance(Value(study.cfg,"commonTolerance"));
    for(int j=0;j<14;++j){int index=j%7;double step=index<4?Value(study.cfg,"commonMeanStep"):(index==4?Value(study.cfg,"commonWidthLogStep"):Value(study.cfg,"commonWidthShapeStep"));m->SetVariable(j,Form("shared_%d",j),p[j],step);}
    bool ok=m->Minimize();Vec next(m->X(),m->X()+14);double candidate=objective(&next[0]);
    result.accepted=candidate<1.e19&&candidate<=before;
    if(result.accepted){p=next;result.g=FromOptim(p,study.mid,study.half,Value(study.cfg,"minimumSigma"));}
    result.status=m->Status();if(!ok&&result.status==0)result.status=1;result.edm=TMath::Finite(m->Edm())?m->Edm():invalid;result.covStatus=m->CovMatrixStatus();result.covariance.clear();result.covariancePoint.clear();
    if(result.accepted&&result.covStatus>0){result.covariance.assign(14,Vec(14));result.covariancePoint=p;for(int j=0;j<14;++j)for(int k=0;k<14;++k)result.covariance[j][k]=m->CovMatrix(j,k);}
    delete m;result.coefficientDrift=0;for(int j=0;j<14;++j)result.coefficientDrift=std::max(result.coefficientDrift,std::fabs(p[j]-old[j]));
    for(unsigned j=0;j<study.rows.size();++j)study.Refresh(study.rows[j],result.g);
    double value=study.CommonLoss(result.g);result.trace.push_back(value);
    bool stable=previous<1.e19&&std::fabs(previous-value)<=Value(study.cfg,"commonRelativeTolerance")*std::max(1.,std::fabs(value));
    result.outerConverged=stable&&result.coefficientDrift<=Value(study.cfg,"commonCoefficientTolerance")&&result.status==0&&result.accepted&&AllLocalValid(study);
    std::cout<<"ITERATIVE_COMMON outer="<<outer<<" dataDeviance="<<value<<" status="<<result.status<<" drift="<<result.coefficientDrift<<" localValid="<<AllLocalValid(study)<<std::endl;
    if(result.outerConverged)break;previous=value;
  }
  // Keep the retained T/BG shape point for conditional covariance, and use a
  // separate final block only to diagnose remaining local improvement.
  Study probe=study;result.localImprovementPerCount=0;
  for(unsigned j=0;j<probe.rows.size();++j){double before=study.rows[j].dataD;FitLocal(probe,probe.rows[j],result.g);
    result.localImprovementPerCount=std::max(result.localImprovementPerCount,std::max(0.,before-probe.rows[j].dataD)/std::max(1.,study.rows[j].raw.count));
  }
  result.finalLocalStationary=AllLocalValid(probe)&&result.localImprovementPerCount<=Value(study.cfg,"commonRelativeTolerance");
  double center=study.CommonLoss(result.g),count=0,step=Value(study.cfg,"commonCoefficientTolerance");result.gradientAbsolute=0;result.stationarityEvaluable=center<1.e19&&step>0;
  for(unsigned j=0;j<study.rows.size();++j)count+=study.rows[j].raw.count;
  Vec point=ToOptim(result.g,study.mid,study.half,Value(study.cfg,"minimumSigma"));
  for(int k=0;k<14&&result.stationarityEvaluable;++k){
    Vec up=point,down=point;up[k]+=step;down[k]-=step;Vec gu=FromOptim(up,study.mid,study.half,Value(study.cfg,"minimumSigma")),gd=FromOptim(down,study.mid,study.half,Value(study.cfg,"minimumSigma"));
    bool plusConstraint=study.SharedFeasible(gu,study.rows[0].h->GetXaxis()->GetXmin(),study.rows[0].h->GetXaxis()->GetXmax());
    bool minusConstraint=study.SharedFeasible(gd,study.rows[0].h->GetXaxis()->GetXmin(),study.rows[0].h->GetXaxis()->GetXmax());
    double vp=plusConstraint?study.CommonLoss(gu):invalid,vm=minusConstraint?study.CommonLoss(gd):invalid;
    // Feasible AMP/profile failure is NOT a one-sided physical boundary.
    if((plusConstraint&&vp>=1.e19)||(minusConstraint&&vm>=1.e19)){++result.profileFailures;result.stationarityEvaluable=false;continue;}
    if(plusConstraint&&minusConstraint){++result.centralDirections;result.gradientAbsolute=std::max(result.gradientAbsolute,std::fabs(vp-vm)/(2*step));}
    else if(plusConstraint){++result.oneSidedDirections;result.gradientAbsolute=std::max(result.gradientAbsolute,std::max(0.,center-vp)/step);}
    else if(minusConstraint){++result.oneSidedDirections;result.gradientAbsolute=std::max(result.gradientAbsolute,std::max(0.,center-vm)/step);}
    else{++result.blockedDirections;result.stationarityEvaluable=false;}
  }
  result.gradientPerCount=result.gradientAbsolute/std::max(1.,count);
  result.finalGlobalStationary=result.stationarityEvaluable&&result.gradientPerCount<=Value(study.cfg,"commonRelativeTolerance");
  result.outerConverged=result.outerConverged&&result.finalLocalStationary&&result.finalGlobalStationary;
  return result;
}
void RNodes(Data& d,int order){
  const double x3[]={-.7745966692414834,0,.7745966692414834},w3[]={5./18,4./9,5./18};
  const double x5[]={-.9061798459386640,-.5384693101056831,0,.5384693101056831,.9061798459386640},w5[]={.1184634425280945,.2393143352496832,.2844444444444444,.2393143352496832,.1184634425280945};
  if(order==1){d.nodes.push_back(d.raw.r);d.nodeWeights.push_back(1);}
  else if(order==3||order==5)for(int j=0;j<order;++j){d.nodes.push_back(d.raw.r+(d.seed.hi-d.seed.lo)/2*(order==3?x3[j]:x5[j]));d.nodeWeights.push_back(order==3?w3[j]:w5[j]);}
  else throw std::runtime_error("R quadrature must be 1/3/5");
}
void LoadData(Study& study,TFile& input,const std::vector<Seed>& seeds){
  GlobalHelium::InputAudit audit;
  for(unsigned j=0;j<seeds.size();++j){Data d;d.seed=seeds[j];d.raw.index=d.seed.index;d.raw.lo=d.seed.lo;d.raw.hi=d.seed.hi;d.raw.r=(d.seed.lo+d.seed.hi)/2;
    TH1D* original=dynamic_cast<TH1D*>(input.Get(Form("raw_native_%03d",d.seed.index)));if(!original)throw std::runtime_error("missing frozen raw_native seed slice");
    d.native=static_cast<TH1D*>(original->Clone(Form("raw_native_%03d",d.seed.index)));d.native->SetDirectory(0);
    for(int b=0;b<=d.native->GetNbinsX()+1;++b)GlobalHelium::AuditCount(d.native->GetBinContent(b),std::pow(d.native->GetBinError(b),2),study.cfg,audit);
    d.h=static_cast<TH1D*>(d.native->Clone(Form("raw_%03d",d.seed.index)));d.h->SetDirectory(0);
    double factorValue=Value(study.cfg,"dedxRebinFactor");int factor=int(factorValue);if(factor<1||factor!=factorValue||d.h->GetNbinsX()%factor)throw std::runtime_error("invalid exact dE/dx rebin");d.h->Rebin(factor);
    double lo=std::max(d.h->GetXaxis()->GetXmin(),d.seed.anchorMean+Value(study.cfg,"fitNSigmaMin")*d.seed.anchorSigma),hi=std::min(d.h->GetXaxis()->GetXmax(),d.seed.anchorMean+Value(study.cfg,"fitNSigmaMax")*d.seed.anchorSigma);
    int first=1,last=d.h->GetNbinsX();while(first<=last&&d.h->GetBinCenter(first)<lo)++first;while(last>=first&&d.h->GetBinCenter(last)>=hi)--last;if(last<first)throw std::runtime_error("empty fixed fit window");
    d.raw.fitLo=d.h->GetBinLowEdge(first);d.raw.fitHi=d.h->GetBinLowEdge(last)+d.h->GetBinWidth(last);
    for(int b=first;b<=last;++b){GlobalHelium::Bin bin;bin.lo=d.h->GetBinLowEdge(b);bin.hi=bin.lo+d.h->GetBinWidth(b);bin.n=d.h->GetBinContent(b);GlobalHelium::AuditCount(bin.n,std::pow(d.h->GetBinError(b),2),study.cfg,audit);d.raw.count+=bin.n;d.raw.bins.push_back(bin);}
    if(!(d.raw.count>0))throw std::runtime_error("zero fit-window count");
    RNodes(d,int(Value(study.cfg,"rQuadratureOrder")));d.shape.slope=d.seed.bgSlope;
    for(int sp=0;sp<3;++sp){d.shape.mu[sp]=d.seed.mu[sp];d.shape.sigma[sp]=d.seed.sigma[sp];d.included[sp]=sp>0||d.seed.tPresent;d.shapeFree[sp]=study.phase!="common"||sp==0;}
    if(d.included[0]){double acceptance=Gaussian(d.raw.fitLo,d.raw.fitHi,d.shape.mu[0],d.shape.sigma[0]);if(acceptance<=Value(study.cfg,"minimumGaussianAcceptance")){d.included[0]=false;d.shapeFree[0]=false;}else if(d.seed.yield[0]<Value(study.cfg,"minimumComponentWindowYield"))d.shapeFree[0]=false;}
    if(study.phase=="regularized"){std::map<int,Weight>::const_iterator w=study.target.weights.find(d.seed.index);if(w==study.target.weights.end())throw std::runtime_error("missing fixed regularizer weights");d.weight=w->second;}
    int amplitudes=d.seed.bg=="linear"?5:4;d.a.assign(amplitudes,0);for(int sp=0;sp<3;++sp)d.a[sp]=d.included[sp]?d.seed.yield[sp]/d.raw.count:0;d.a[3]=d.seed.bgYield/d.raw.count;if(amplitudes==5){d.a[3]*=.5;d.a[4]=d.a[3];}
    if(study.phase=="common")for(int sp=1;sp<3;++sp){d.shape.mu[sp]=SharedMean(study.target.coefficients,sp-1,d.raw.r);d.shape.sigma[sp]=SharedSigma(study.target.coefficients,sp-1,d.raw.r);}
    study.rows.push_back(d);
  }
  std::cout<<"ITERATIVE_INPUT countVarianceCells="<<audit.cells<<" maximumIntegerResidual="<<audit.maximumIntegerResidual<<" maximumVarianceRelativeResidual="<<audit.maximumVarianceResidual<<std::endl;
}
void Protect(const TString& out,const char* settings,const char* seeds,const char* targets){
  const char* names[]={"fit_result.json","phase.root","slices.tsv","completed.json"};
  for(unsigned j=0;j<sizeof(names)/sizeof(names[0]);++j)if(!gSystem->AccessPathName(out+"/"+names[j]))throw std::runtime_error("existing iterative output will not be replaced");
  if(!gSystem->AccessPathName(out)){
    bool recorded=!gSystem->AccessPathName(out+"/context.json")||!gSystem->AccessPathName(out+"/../context.json");
    if(!recorded||TString(settings)!=out+"/settings.tsv"||TString(seeds)!=out+"/seeds.tsv"||TString(targets)!=out+"/targets.tsv")throw std::runtime_error("existing output is not a prepared iterative study");
  }else if(gSystem->mkdir(out,kTRUE))throw std::runtime_error("cannot create new iterative output");
}
void WriteMatrix(std::ostream& o,const Matrix& m){o<<'[';for(unsigned j=0;j<m.size();++j){if(j)o<<',';WriteVec(o,m[j]);}o<<']';}
void WriteProducts(Study& study,const TString& out,const CommonResult& common){
  TFile root(out+"/phase.root","CREATE");if(root.IsZombie())throw std::runtime_error("cannot CREATE phase ROOT");
  std::ofstream json((out+"/fit_result.json").Data()),tsv((out+"/slices.tsv").Data());if(!json||!tsv)throw std::runtime_error("cannot CREATE phase result/manifest");
  json<<std::setprecision(17);tsv<<std::setprecision(17);
  tsv<<"index\tlo\thi\tr\tanchorMean\tanchorSigma\tfitLo\tfitHi\tfirstFitBin\tlastFitBin\tstatus\tcovarianceStatus\tquality\tatBoundary\tdataDeviance\tpenaltyLoss\tndf\toldChosen\toldValid\toldQuality\trawSeparatePeak\tchi2\tfinite\tampConverged\n";
  double totalD=0,totalQ=0,totalCounts=0;bool numerical=true;int finiteRows=0,ampRows=0,status0Rows=0;
  for(unsigned j=0;j<study.rows.size();++j){totalD+=study.rows[j].dataD;totalQ+=study.rows[j].penalty;totalCounts+=study.rows[j].raw.count;finiteRows+=study.rows[j].finite; ampRows+=study.rows[j].ampConverged;status0Rows+=study.rows[j].status==0;numerical=numerical&&study.rows[j].finite&&study.rows[j].ampConverged&&study.rows[j].status==0&&study.rows[j].accepted;}
  json<<"{\"schemaVersion\":1,\"phase\":"<<Quote(study.phase)<<",\"penaltyStrength\":"<<Value(study.cfg,"penaltyStrength")
    <<",\"dataDeviance\":"<<totalD<<",\"penaltyLoss\":"<<totalQ<<",\"totalLoss\":"<<totalD+totalQ
    <<",\"executionIsNotPhysicalValidation\":true,\"numericalConverged\":"<<(numerical&&(study.phase!="common"||common.outerConverged)?"true":"false")
    <<",\"widthMeaning\":"<<Quote(study.phase=="common"?"He3/He4 intrinsic kernel with uniform R-bin GL average; T remains centre-only projection nuisance":"centre-only R-projection width (folded), NOT intrinsic width")
    <<",\"labelPolicy\":\"free-only ascending-mean canonicalization; permutation saved; no isotope truth or legacy ordering constraint\""
    <<",\"yieldDefinition\":\"nonnegative Gaussian/smooth counts CONDITIONAL on the actual fit-bin domain; full Gaussian yield derived by dividing window acceptance, not an efficiency\""
    <<",\"effectiveCountDefinition\":\"sum(componentExpectedCount^2/totalExpectedCount) heuristic information, not independently measured signal count\""
    <<",\"qualityScope\":\"numerical status/AMP/finite/boundaries only, NOT model adequacy or isotope identification\""
    <<",\"nSigmaAnchorPolicy\":\"original frozen from040 anchor fixed; counts/errors unchanged, no y-Jacobian\""
    <<",\"rQuadratureOrder\":"<<Value(study.cfg,"rQuadratureOrder")<<",\"sharedCoefficients\":";
  if(study.phase=="common")WriteVec(json,common.g);else if(study.target.present)WriteVec(json,study.target.coefficients);else json<<"null";
  json<<",\"summary\":{\"dataDeviance\":"<<totalD<<",\"penaltyLoss\":"<<totalQ<<",\"totalLoss\":"<<totalD+totalQ<<",\"slices\":"<<study.rows.size()<<",\"finiteSlices\":"<<finiteRows<<",\"ampConvergedSlices\":"<<ampRows<<",\"localStatus0Slices\":"<<status0Rows<<",\"totalFitCount\":"<<totalCounts<<",\"numericalConverged\":"<<(numerical&&(study.phase!="common"||common.outerConverged)?"true":"false")<<",\"commonStatus\":"<<common.status<<",\"commonOuterConverged\":"<<(common.outerConverged?"true":"false")<<"}";
  json<<",\"common\":{\"status\":"<<common.status<<",\"pointAccepted\":"<<(common.accepted?"true":"false")<<",\"edm\":"<<common.edm<<",\"outerConverged\":"<<(common.outerConverged?"true":"false")
    <<",\"normalizedCoefficientDrift\":"<<common.coefficientDrift
    <<",\"finalLocalStationary\":"<<(common.finalLocalStationary?"true":"false")<<",\"finalLocalImprovementPerCount\":"<<common.localImprovementPerCount
    <<",\"finalGlobalStationary\":"<<(common.finalGlobalStationary?"true":"false")<<",\"stationarityEvaluable\":"<<(common.stationarityEvaluable?"true":"false")
    <<",\"stationarityAbsoluteResidual\":"<<common.gradientAbsolute<<",\"stationarityResidualPerCount\":"<<common.gradientPerCount
    <<",\"stationarityCentralDirections\":"<<common.centralDirections<<",\"stationarityOneSidedDirections\":"<<common.oneSidedDirections
    <<",\"stationarityBlockedDirections\":"<<common.blockedDirections<<",\"stationarityProfileFailures\":"<<common.profileFailures
    <<",\"stationarityPolicy\":\"axis necessary-condition diagnostic in normalized optimizer coordinates; step=commonCoefficientTolerance; per-count tolerance=commonRelativeTolerance; NOT all feasible tangent directions\""
    <<",\"trace\":";WriteVec(json,common.trace);
  json<<",\"covarianceStatus\":"<<common.covStatus<<",\"covarianceCoordinates\":\"normalized cubic mu; log positive minimum-offset, slope and curvature for quadratic sigma\""
    <<",\"covarianceScope\":\"conditional on retained T and nonlinear background shapes; yields profiled; NOT full nuisance or smoothing uncertainty\""
    <<",\"covarianceEvaluationPoint\":";if(common.covariancePoint.empty())json<<"null";else WriteVec(json,common.covariancePoint);
  json<<",\"covariance\":";if(common.covariance.empty())json<<"null";else WriteMatrix(json,common.covariance);json<<"},\"slices\":[";
  for(unsigned j=0;j<study.rows.size();++j){
    Data& d=study.rows[j];Vec acceptance;Matrix t=study.Templates(d,d.shape,common.g,&acceptance);
    if(t.empty())throw std::runtime_error("returned point has no finite prediction templates");
    const int first=d.h->FindBin(d.raw.fitLo+Value(study.cfg,"axisTolerance")),last=first+d.raw.bins.size()-1;
    double count[3]={0,0,0},effective[3]={0,0,0},background=0;
    for(int sp=0;sp<3;++sp)count[sp]=d.raw.count*d.a[sp];for(unsigned k=3;k<d.a.size();++k)background+=d.raw.count*d.a[k];
    root.cd();d.native->Write();d.h->Write();
    TH1D projection(Form("projection_%03d",d.seed.index),"fixed old anchor",d.h->GetNbinsX(),(d.h->GetXaxis()->GetXmin()-d.seed.anchorMean)/d.seed.anchorSigma,(d.h->GetXaxis()->GetXmax()-d.seed.anchorMean)/d.seed.anchorSigma);
    projection.Sumw2();for(int b=0;b<=d.h->GetNbinsX()+1;++b){projection.SetBinContent(b,d.h->GetBinContent(b));projection.SetBinError(b,d.h->GetBinError(b));}projection.SetEntries(d.h->GetEntries());projection.Write();
    const char* labels[]={"total","smooth","t","he3","he4"};std::vector<TH1D*> predictions;
    for(int p=0;p<5;++p){TH1D* h=static_cast<TH1D*>(d.h->Clone(Form("prediction_%s_%03d",labels[p],d.seed.index)));h->Reset();h->SetDirectory(0);predictions.push_back(h);}
    TH1D* used=static_cast<TH1D*>(d.h->Clone(Form("fit_used_%03d",d.seed.index)));used->Reset();used->SetDirectory(0);
    for(unsigned b=0;b<t.size();++b){double total=0,smooth=0,c[3];for(int sp=0;sp<3;++sp){c[sp]=count[sp]*t[b][sp];total+=c[sp];}
      for(unsigned k=3;k<d.a.size();++k)smooth+=d.raw.count*d.a[k]*t[b][k];total+=smooth;
      for(int sp=0;sp<3;++sp)effective[sp]+=c[sp]*c[sp]/std::max(1.e-20,total);
      predictions[0]->SetBinContent(first+b,total);predictions[1]->SetBinContent(first+b,smooth);
      for(int sp=0;sp<3;++sp)predictions[sp+2]->SetBinContent(first+b,c[sp]);used->SetBinContent(first+b,1);
    }
    for(unsigned p=0;p<predictions.size();++p){predictions[p]->Write();delete predictions[p];}used->Write();delete used;
    bool valid=d.finite&&d.ampConverged&&d.status==0&&d.accepted;bool quality=valid&&!d.boundary&&d.covStatus==3;
    int amplitudeDof=0;for(unsigned k=0;k<d.a.size();++k){double sum=0;for(unsigned b=0;b<t.size();++b)sum+=t[b][k];if(sum>0)++amplitudeDof;}
    int ndf=int(d.raw.bins.size())-d.freeCount-amplitudeDof;
    if(j)json<<',';json<<"{\"index\":"<<d.seed.index<<",\"lo\":"<<d.seed.lo<<",\"hi\":"<<d.seed.hi<<",\"r\":"<<d.raw.r
      <<",\"anchorMean\":"<<d.seed.anchorMean<<",\"anchorSigma\":"<<d.seed.anchorSigma
      <<",\"requestedFitLo\":"<<std::max(d.h->GetXaxis()->GetXmin(),d.seed.anchorMean+Value(study.cfg,"fitNSigmaMin")*d.seed.anchorSigma)
      <<",\"requestedFitHi\":"<<std::min(d.h->GetXaxis()->GetXmax(),d.seed.anchorMean+Value(study.cfg,"fitNSigmaMax")*d.seed.anchorSigma)
      <<",\"fitLo\":"<<d.raw.fitLo<<",\"fitHi\":"<<d.raw.fitHi<<",\"firstFitBin\":"<<first<<",\"lastFitBin\":"<<last
      <<",\"fitCount\":"<<d.raw.count<<",\"oldChosen\":"<<d.seed.oldChosen<<",\"oldValid\":"<<(d.seed.oldValid?"true":"false")<<",\"oldQuality\":"<<(d.seed.oldQuality?"true":"false")
      <<",\"backgroundModel\":"<<Quote(d.seed.bg)<<",\"backgroundSlope\":"<<d.shape.slope<<",\"backgroundWindowYield\":"<<background
      <<",\"finite\":"<<(d.finite?"true":"false")<<",\"valid\":"<<(valid?"true":"false")<<",\"quality\":"<<(quality?"true":"false")
      <<",\"atBoundary\":"<<(d.boundary?"true":"false")<<",\"status\":"<<d.status<<",\"covarianceStatus\":"<<d.covStatus<<",\"edm\":"<<d.edm
      <<",\"pointAccepted\":"<<(d.accepted?"true":"false")<<",\"dataDeviance\":"<<d.dataD<<",\"penaltyLoss\":"<<d.penalty<<",\"totalLoss\":"<<d.dataD+d.penalty
      <<",\"chi2\":"<<d.pearson<<",\"ndf\":"<<ndf<<",\"ndfScope\":\"conditional nominal slice bins minus locally free shapes and nonzero amplitude templates; not global14 NDF\""
      <<",\"ampConverged\":"<<(d.ampConverged?"true":"false")<<",\"amplitudeKktResidual\":"<<d.ampKkt<<",\"amplitudeIterations\":"<<d.ampIterations
      <<",\"rawSeparatePeak\":null,\"labelSwapped\":"<<(d.labelSwapped?"true":"false")<<",\"labelPermutation\":";WriteVec(json,d.labelPermutation);
    json<<",\"localCovarianceCoordinates\":\"fitted local mean/logSigma/expSlope order, amplitudes profiled; null after relabel or common global update\""
      <<",\"localCovariance\":";if(d.covariance.empty()||study.phase=="common")json<<"null";else WriteMatrix(json,d.covariance);
    const char* component[]={"t","he3","he4"};
    for(int sp=0;sp<3;++sp){
      double scale=study.phase=="regularized"&&sp>0?d.weight.widthScale[sp-1]:d.shape.sigma[sp];
      json<<','<<Quote(component[sp])<<":{\"mean\":"<<d.shape.mu[sp]<<",\"sigma\":"<<d.shape.sigma[sp]<<",\"windowYield\":"<<count[sp]
        <<",\"effectiveCount\":"<<effective[sp]<<",\"penaltyScale\":"<<scale<<",\"windowAcceptance\":"<<acceptance[sp]<<",\"fullYield\":";
      if(acceptance[sp]>Value(study.cfg,"minimumGaussianAcceptance"))json<<count[sp]/acceptance[sp];else json<<"null";
      json<<",\"included\":"<<(d.included[sp]?"true":"false")<<",\"shapeFree\":"<<(d.shapeFree[sp]?"true":"false")
        <<",\"shapeState\":"<<Quote(!d.included[sp]?"inactive":(!d.shapeFree[sp]?(sp>0&&study.phase=="common"?"shared_fixed":"fixed_unidentified"):(study.phase=="regularized"&&sp>0?"regularized":"free")))
        <<",\"originalLabel\":"<<d.labelPermutation[sp];
      if(sp>0)json<<",\"regularizerInformation\":"<<d.weight.info[sp-1]<<",\"regularizerMeanScale\":"<<d.weight.meanScale[sp-1]<<",\"regularizerWidthScale\":"<<d.weight.widthScale[sp-1];
      json<<'}';
    }json<<'}';
    tsv<<d.seed.index<<'\t'<<d.seed.lo<<'\t'<<d.seed.hi<<'\t'<<d.raw.r<<'\t'<<d.seed.anchorMean<<'\t'<<d.seed.anchorSigma<<'\t'<<d.raw.fitLo<<'\t'<<d.raw.fitHi<<'\t'<<first<<'\t'<<last<<'\t'<<d.status<<'\t'<<d.covStatus<<'\t'<<quality<<'\t'<<d.boundary<<'\t'<<d.dataD<<'\t'<<d.penalty<<'\t'<<ndf<<'\t'<<d.seed.oldChosen<<'\t'<<d.seed.oldValid<<'\t'<<d.seed.oldQuality<<"\t-1\t"<<d.pearson<<'\t'<<d.finite<<'\t'<<d.ampConverged<<'\n';
  }json<<"]}\n";json.close();tsv.close();root.Close();
}
}  // namespace HeliumIterative
#endif
void IterativeHeliumFit(const char* cacheRoot,const char* output,const char* settings,const char* seeds,const char* targets){
#ifndef __CINT__
  using namespace HeliumIterative;
  try{
    Settings cfg=ReadSettings(settings);Target target=ReadTarget(targets);Study study(cfg,target);
    Protect(TString(output),settings,seeds,targets);if(gSystem->Load("libMinuit2")<0)throw std::runtime_error("cannot load Minuit2");
    TFile cache(cacheRoot,"READ");if(cache.IsZombie())throw std::runtime_error("cannot open frozen from040 ROOT");
    LoadData(study,cache,ReadSeeds(seeds));CommonResult common;
    if(study.phase=="common")common=FitCommon(study);
    else for(unsigned j=0;j<study.rows.size();++j){
      FitLocal(study,study.rows[j],target.coefficients);
      if(study.phase=="free"){CanonicalizeFree(study.rows[j]);study.Refresh(study.rows[j],target.coefficients);}
      std::cout<<"ITERATIVE_LOCAL index="<<study.rows[j].seed.index<<" dataDeviance="<<study.rows[j].dataD<<" penalty="<<study.rows[j].penalty<<" status="<<study.rows[j].status<<" AMP="<<study.rows[j].ampConverged<<std::endl;
    }
    WriteProducts(study,TString(output),common);for(unsigned j=0;j<study.rows.size();++j){delete study.rows[j].native;delete study.rows[j].h;}
    std::cout<<"ITERATIVE_HELIUM_COMPLETE phase="<<study.phase<<" slices="<<study.rows.size()<<std::endl;
  }catch(const std::exception& e){std::cerr<<"IterativeHeliumFit: "<<e.what()<<std::endl;gSystem->Exit(2);}
#endif
}
