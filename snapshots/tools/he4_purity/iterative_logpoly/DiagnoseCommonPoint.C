// Read-only frozen-point diagnostics. Never fit shapes, change PID, or mutate core.
void DiagnoseCommonPoint(const char* cacheRoot,const char* phaseRoot,const char* output,
                         const char* settings,const char* seeds,const char* targets,
                         const char* amplitudes,const char* variants);
#ifndef __CINT__
#include "IterativeHeliumFit.C"
namespace HeliumPointDiagnostic {
using namespace HeliumIterative;
struct Variant {std::string name;double tolerance;int iterations;};
std::vector<Variant> ReadVariants(const char* path){
  std::ifstream f(path);if(!f)throw std::runtime_error("missing explicit diagnostic variants");
  std::vector<Variant> out;Variant v;std::set<std::string> names;
  while(f>>v.name>>v.tolerance>>v.iterations){
    if(!names.insert(v.name).second||!TMath::Finite(v.tolerance)||!(v.tolerance>0)||v.iterations<1)throw std::runtime_error("invalid diagnostic variant");
    out.push_back(v);
  }
  if(!f.eof()||out.empty())throw std::runtime_error("malformed/empty diagnostic variants");return out;
}
void Scalar(std::ostream& out,double x){if(TMath::Finite(x))out<<x;else out<<"null";}
void RestoreAmplitudes(Study& s,TFile& phase,const char* path){
  std::ifstream f(path);if(!f)throw std::runtime_error("missing final Gaussian/window amplitudes");
  std::map<int,Vec> rows;int index;Vec a(4); // T, He3, He4 and total smooth yield.
  while(f>>index>>a[0]>>a[1]>>a[2]>>a[3]){if(!rows.insert(std::make_pair(index,a)).second)throw std::runtime_error("duplicate amplitude index");}
  if(!f.eof())throw std::runtime_error("malformed amplitude row");
  for(unsigned j=0;j<s.rows.size();++j){Data& d=s.rows[j];std::map<int,Vec>::const_iterator row=rows.find(d.seed.index);
    if(row==rows.end())throw std::runtime_error("missing retained final amplitude");
    for(int sp=0;sp<3;++sp)d.a[sp]=row->second[sp]/d.raw.count;
    if(d.seed.bg=="linear"){
      TH1D* h=dynamic_cast<TH1D*>(phase.Get(Form("prediction_smooth_%03d",d.seed.index)));if(!h)throw std::runtime_error("missing saved linear BG prediction");
      const int first=h->FindBin(d.raw.fitLo+Value(s.cfg,"axisTolerance")),last=first+d.raw.bins.size()-1;
      if(last<=first)throw std::runtime_error("linear background requires two saved bins");
      double width=d.raw.fitHi-d.raw.fitLo;
      double x0=h->GetBinCenter(first),x1=h->GetBinCenter(last),dx0=h->GetBinWidth(first),dx1=h->GetBinWidth(last);
      double u0=2*dx0*(d.raw.fitHi-x0)/(width*width),v0=2*dx0*(x0-d.raw.fitLo)/(width*width);
      double u1=2*dx1*(d.raw.fitHi-x1)/(width*width),v1=2*dx1*(x1-d.raw.fitLo)/(width*width);
      double det=u0*v1-u1*v0,b0=h->GetBinContent(first),b1=h->GetBinContent(last);
      if(!(det>0))throw std::runtime_error("singular saved linear BG basis");
      double n3=(b0*v1-b1*v0)/det,n4=(u0*b1-u1*b0)/det;
      if(n3< -Value(s.cfg,"countIntegerTolerance")||n4< -Value(s.cfg,"countIntegerTolerance"))throw std::runtime_error("negative reconstructed linear BG yield");
      d.a[3]=std::max(0.,n3)/d.raw.count;d.a[4]=std::max(0.,n4)/d.raw.count;
    }else d.a[3]=row->second[3]/d.raw.count;
  }
}
Vec Acceptances(const Study& s,const Data& d,const Vec& g){
  Vec a(3,0);for(int sp=0;sp<3;++sp)if(d.included[sp]){
    if(sp==0)a[sp]=Gaussian(d.raw.fitLo,d.raw.fitHi,d.shape.mu[sp],d.shape.sigma[sp]);
    else for(unsigned q=0;q<d.nodes.size();++q)a[sp]+=d.nodeWeights[q]*Gaussian(d.raw.fitLo,d.raw.fitHi,SharedMean(g,sp-1,d.nodes[q]),SharedSigma(g,sp-1,d.nodes[q]));
  }return a;
}
void Evaluate(const Study& s,const Vec& g,const std::string& name,std::ostream& out){
  const bool shared=s.SharedFeasible(g,s.rows[0].h->GetXaxis()->GetXmin(),s.rows[0].h->GetXaxis()->GetXmax());
  std::map<std::string,int> reasons;double total=0,maxKkt=0;int failed=0;
  out<<"{\"case\":"<<Quote(name)<<",\"sharedFeasible\":"<<(shared?"true":"false")<<",\"coefficients\":";WriteVec(out,g);out<<",\"slices\":[";
  for(unsigned j=0;j<s.rows.size();++j){const Data& d=s.rows[j];Vec a=Acceptances(s,d,g);std::string reason="ok";Profile p;
    bool local=s.ShapeFeasible(d,d.shape);Matrix t;
    if(!shared)reason="shared_support";
    else if(!local)reason="local_shape_support";
    else {t=s.Templates(d,d.shape,g);if(t.empty())reason="gaussian_acceptance_or_templates";
      else {p=WindowProfile(d.raw,t,d.a,s.cfg);if(!p.converged)reason="amplitude_kkt";}}
    ++reasons[reason];if(reason!="ok")++failed;else total+=p.d;
    if(!t.empty()&&TMath::Finite(p.kkt))maxKkt=std::max(maxKkt,p.kkt);
    if(j)out<<',';out<<"{\"index\":"<<d.seed.index<<",\"r\":"<<d.raw.r<<",\"reason\":"<<Quote(reason)
      <<",\"shapeFeasible\":"<<(local?"true":"false")<<",\"fitCount\":"<<d.raw.count<<",\"gaussianAcceptances\":";WriteVec(out,a);
    out<<",\"included\":[";for(int sp=0;sp<3;++sp){if(sp)out<<',';out<<(d.included[sp]?"true":"false");}out<<']';
    out<<",\"amplitudeIterations\":"<<p.iterations<<",\"amplitudeConverged\":"<<(p.converged?"true":"false")<<",\"amplitudeKktResidual\":";
    if(!t.empty())Scalar(out,p.kkt);else out<<"null";
    out<<",\"dataDeviance\":";if(!t.empty())Scalar(out,p.d);else out<<"null";out<<'}';
  }
  out<<"],\"valid\":"<<(!failed?"true":"false")<<",\"objective\":";if(failed)out<<"null";else Scalar(out,total);
  out<<",\"validSliceDevianceSum\":"<<total<<",\"failedSlices\":"<<failed<<",\"maximumAmplitudeKktResidual\":"<<maxKkt<<",\"reasonCounts\":{";
  int k=0;for(std::map<std::string,int>::const_iterator i=reasons.begin();i!=reasons.end();++i){if(k++)out<<',';out<<Quote(i->first)<<':'<<i->second;}out<<"}}";
}
}
#endif
void DiagnoseCommonPoint(const char* cacheRoot,const char* phaseRoot,const char* output,
                         const char* settings,const char* seeds,const char* targets,
                         const char* amplitudes,const char* variants){
#ifndef __CINT__
  using namespace HeliumPointDiagnostic;
  try{
    const TString out(output);if(!gSystem->AccessPathName(out+"/diagnostic.json"))throw std::runtime_error("diagnostic output already exists");
    if(gSystem->AccessPathName(out+"/context.json"))throw std::runtime_error("explicit prepared diagnostic context required");
    Settings base=ReadSettings(settings);Target target=ReadTarget(targets);std::vector<Seed> seed=ReadSeeds(seeds);std::vector<Variant> variant=ReadVariants(variants);
    if(Text(base,"phase")!="common"||!target.present||Value(base,"penaltyStrength")!=0)throw std::runtime_error("only frozen unpenalized common-point diagnosis");
    TFile cache(cacheRoot,"READ"),phase(phaseRoot,"READ");if(cache.IsZombie()||phase.IsZombie())throw std::runtime_error("unreadable frozen cache/phase");
    std::ofstream json((out+"/diagnostic.json").Data());if(!json)throw std::runtime_error("cannot create diagnostic JSON");
    json<<std::setprecision(17)<<"{\"schemaVersion\":1,\"shapeRefit\":false,\"productionEnabled\":false,\"probeCoordinates\":\"exact normalized optimizer coordinates from frozen ToOptim/FromOptim\",\"stepFromSetting\":\"commonCoefficientTolerance\",\"variants\":[";
    for(unsigned v=0;v<variant.size();++v){Settings cfg=base;std::ostringstream tolerance,iterations;tolerance<<std::setprecision(17)<<variant[v].tolerance;iterations<<variant[v].iterations;
      cfg["amplitudeGradientTolerance"]=tolerance.str();cfg["amplitudeMaxIterations"]=iterations.str();Study study(cfg,target);LoadData(study,cache,seed);RestoreAmplitudes(study,phase,amplitudes);
      Vec point=ToOptim(target.coefficients,study.mid,study.half,Value(cfg,"minimumSigma"));const double step=Value(cfg,"commonCoefficientTolerance");
      if(v)json<<',';json<<"{\"name\":"<<Quote(variant[v].name)<<",\"amplitudeTolerance\":"<<variant[v].tolerance<<",\"amplitudeMaxIterations\":"<<variant[v].iterations<<",\"cases\":[";
      Evaluate(study,target.coefficients,"base",json);
      for(int k=0;k<14;++k)for(int sign=-1;sign<=1;sign+=2){Vec shifted=point;shifted[k]+=sign*step;json<<',';Evaluate(study,FromOptim(shifted,study.mid,study.half,Value(cfg,"minimumSigma")),Form("coordinate_%d_%s",k,sign<0?"minus":"plus"),json);}
      json<<"]}";for(unsigned j=0;j<study.rows.size();++j){delete study.rows[j].native;delete study.rows[j].h;}
      std::cout<<"COMMON_POINT_DIAGNOSTIC_VARIANT "<<variant[v].name<<" cases=29"<<std::endl;
    }json<<"]}\n";json.close();std::cout<<"COMMON_POINT_DIAGNOSTIC_COMPLETE"<<std::endl;
  }catch(const std::exception& e){std::cerr<<"DiagnoseCommonPoint: "<<e.what()<<std::endl;gSystem->Exit(2);}
#endif
}
