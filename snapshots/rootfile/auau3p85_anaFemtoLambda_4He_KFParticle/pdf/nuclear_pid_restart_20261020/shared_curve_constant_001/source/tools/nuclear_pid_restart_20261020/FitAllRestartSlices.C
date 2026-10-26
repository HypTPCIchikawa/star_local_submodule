// Fixed common polynomial responses; ROOT5 batch entry, study only.
void FitAllRestartSlices(const char*,const char*,const char*,const char*,const char*);
#ifndef __CINT__
#include "FitStorageSlices.C"
#include "AllFixedResponseFit.h"
namespace RestartAllIO {
typedef std::vector<RestartAllFixed::Shape> Shapes;
inline const char* Species(unsigned i){const char* names[]={"d","t","he3","he4"};RestartPid::Check(i<4,"Species index");return names[i];}
inline unsigned SpeciesIndex(const std::string& name){for(unsigned j=0;j<4;++j)if(name==Species(j))return j;throw std::runtime_error("Unknown species: "+name);}
inline std::map<std::string,Shapes> ReadNodes(const char* path){
  std::ifstream input(path);RestartPid::Check(input.good(),"Cannot read fixed nodes");std::string line;
  RestartPid::Check(bool(std::getline(input,line))&&line=="id\tspecies\tweight\tmean\tsigma\twidthRatio\tactive","Unexpected fixed-node header");
  std::map<std::string,Shapes> result;
  while(std::getline(input,line)){if(line.empty())continue;std::vector<std::string> row=RestartPid::Split(line,'\t');RestartPid::Check(row.size()==7,"Fixed node needs seven columns");
    RestartPid::SafeId(row[0]);unsigned j=SpeciesIndex(row[1]);double w=RestartPid::Number(row[2],"weight"),m=RestartPid::Number(row[3],"mean"),s=RestartPid::Number(row[4],"sigma"),ratio=RestartPid::Number(row[5],"ratio");
    RestartPid::Check(w>0.&&s>0.&&ratio>0.&&(row[6]=="0"||row[6]=="1"),"Invalid fixed node");if(result[row[0]].empty())result[row[0]].resize(4);
    RestartAllFixed::Shape& shape=result[row[0]][j];bool active=row[6]=="1";RestartPid::Check(!shape.supportProvided||shape.component.active==active,"Inconsistent activity");
    shape.supportProvided=true;shape.component.name=Species(j);shape.component.active=active;shape.component.freeMean=false;shape.component.freeSigma=false;
    shape.component.nodes.push_back(NuclearFit::Node(w,m,s));shape.ratios.push_back(ratio);
  }
  RestartPid::Check(input.eof()&&!result.empty(),"Nodes empty/read failure");
  for(std::map<std::string,Shapes>::iterator row=result.begin();row!=result.end();++row)for(unsigned j=0;j<4;++j){RestartAllFixed::Shape& s=row->second[j];
    s.component.name=Species(j);s.component.freeMean=false;s.component.freeSigma=false;if(!s.supportProvided){s.component.active=false;continue;}
    double sum=0.;s.component.mean=0.;s.component.sigma=0.;for(unsigned q=0;q<s.ratios.size();++q)sum+=s.component.nodes[q].weight;RestartPid::Check(sum>0.,"Invalid node weights");
    for(unsigned q=0;q<s.ratios.size();++q){NuclearFit::Node& n=s.component.nodes[q];n.weight/=sum;s.component.mean+=n.weight*n.mean;s.component.sigma+=n.weight*n.sigma;}}
  return result;
}
inline NuclearFit::Request Request(const TH1D& native,const RestartPid::Slice& slice,const Shapes& shapes,const NuclearFit::Settings& cfg){
  NuclearFit::Request r;r.bins=RestartPid::NativeBins(native,slice);r.nativeAxisMin=native.GetXaxis()->GetXmin();r.nativeAxisMax=native.GetXaxis()->GetXmax();
  r.backgroundModel=NuclearFit::Text(cfg,"backgroundModel");r.backgroundSlope=NuclearFit::Number(cfg,"backgroundSlope");r.freeBackgroundSlope=false;
  RestartPid::Check(NuclearFit::Integer(cfg,"freeBackgroundSlope")==0,"Background slope must be fixed");
  double fraction=NuclearFit::Number(cfg,"startTargetFraction");RestartPid::Check(fraction>0.&&fraction<1.,"Invalid amplitude start");
  unsigned active=0;for(unsigned j=0;j<4;++j)active+=shapes[j].component.active;unsigned target=SpeciesIndex(slice.target);
  RestartPid::Check(shapes[target].supportProvided&&shapes[target].component.active&&active,"Unsupported/inactive target");
  double count=NuclearFit::TotalCount(r.bins);r.backgroundWindowYield=count*(1.-fraction);
  for(unsigned j=0;j<4;++j){NuclearFit::Component c=shapes[j].component;c.startWindowYield=c.active?count*fraction/active:0.;r.components.push_back(c);}return r;
}
inline RestartPid::Shapes LegacyShapes(const Shapes& shapes){RestartPid::Shapes result;for(unsigned j=0;j<4;++j){RestartPid::ShapeInput s;s.nodes=shapes[j].component.nodes;s.mean=shapes[j].component.mean;s.sigma=shapes[j].component.sigma;result[Species(j)]=s;}return result;}
inline void IntervalsJson(std::ostream& out,const std::vector<RestartAllFixed::Intervals>& nodes){out<<'[';for(unsigned q=0;q<nodes.size();++q){if(q)out<<',';out<<'[';for(unsigned k=0;k<nodes[q].size();++k){if(k)out<<',';const RestartAllFixed::Interval& x=nodes[q][k];out<<"{\"lo\":";RestartPid::JsonNumber(out,x.lo);out<<",\"hi\":";RestartPid::JsonNumber(out,x.hi);out<<",\"includeLo\":"<<(x.includeLo?"true":"false")<<",\"includeHi\":"<<(x.includeHi?"true":"false")<<'}';}out<<']';}out<<']';}
inline void GateJson(std::ostream& out,const RestartAllFixed::Gate& g){
  out<<"{\"completeSupport\":"<<(g.completeSupport?"true":"false")<<",\"insideFit\":"<<(g.insideFit?"true":"false")<<",\"covarianceValid\":"<<(g.covarianceValid?"true":"false")<<",\"signal\":";RestartPid::JsonNumber(out,g.signal);
  out<<",\"background\":";RestartPid::JsonNumber(out,g.background);out<<",\"total\":";RestartPid::JsonNumber(out,g.total);
  out<<",\"candidatePurity\":";if(g.purity>=0.)RestartPid::JsonNumber(out,g.purity);else out<<"null";
  out<<",\"conditionalError\":";if(g.error>=0.)RestartPid::JsonNumber(out,g.error);else out<<"null";
  out<<",\"fractions\":";RestartPid::JsonVector(out,g.fractions);out<<",\"gateYields\":";RestartPid::JsonVector(out,g.gateYields);out<<",\"gradient\":";RestartPid::JsonVector(out,g.gradient);
  out<<",\"signalBackgroundCovariance\":";RestartPid::JsonMatrix(out,g.countCovariance);out<<",\"nominalIntervalsByNode\":";IntervalsJson(out,g.nominal);out<<",\"effectiveIntervalsByNode\":";IntervalsJson(out,g.effective);out<<'}';
}
inline void ReplaceOne(std::string& text,const std::string& old,const std::string& replacement){size_t at=text.find(old);RestartPid::Check(at!=std::string::npos&&text.find(old,at+old.size())==std::string::npos,"Nonunique JSON metadata replacement");text.replace(at,old.size(),replacement);}
inline std::string ResultJson(const RestartPid::Slice& slice,const Shapes& shapes,const NuclearFit::Request& request,const RestartAllFixed::Result& result,const TH1D& native,size_t index,double rLo,double rHi,const std::string& mapping){
  std::ostringstream base;RestartPid::WriteResult(base,slice,LegacyShapes(shapes),request,result.fit,native,index,rLo,rHi);std::string row=base.str();
  ReplaceOne(row,"\"gateSource\":\"initial_git_curve_diagnostic_only\"","\"gateSource\":\"shared_polynomial_per_GL3_node_strict_nsigma_and_best_species\"");
  std::ostringstream purity;purity<<"\"curveAdoptionAllowed\":false,\"purityAvailable\":"<<(result.purityAvailable?"true":"false")<<",\"purity\":";
  if(result.purityAvailable)RestartPid::JsonNumber(purity,result.winner.purity);else purity<<"null";
  ReplaceOne(row,"\"curveAdoptionAllowed\":false,\"purityAvailable\":false,\"purity\":null",purity.str());RestartPid::Check(!row.empty()&&row[row.size()-1]=='}',"Invalid base JSON");
  std::ostringstream out;out<<row.substr(0,row.size()-1)<<",\"storageMapping\":"<<RestartPid::Quote(mapping)<<",\"effectiveFitLo\":"<<std::setprecision(17)<<request.bins.front().lo<<",\"effectiveFitHi\":"<<request.bins.back().hi;
  out<<",\"scaledFisherInfinityCondition\":";RestartPid::JsonNumber(out,result.informationCondition);out<<",\"targetRelativeConditionalYieldError\":";RestartPid::JsonNumber(out,result.targetRelativeError);
  out<<",\"conditionalDropDeltaD\":";RestartPid::JsonVector(out,result.fit.conditionalDropDeltaD);out<<",\"purityConditionalError\":";if(result.purityAvailable)RestartPid::JsonNumber(out,result.winner.error);else out<<"null";
  out<<",\"purityCandidate\":";if(result.winner.purity>=0.)RestartPid::JsonNumber(out,result.winner.purity);else out<<"null";
  out<<",\"purityQualityReasons\":[";for(unsigned k=0;k<result.reasons.size();++k){if(k)out<<',';out<<RestartPid::Quote(result.reasons[k]);}
  out<<"],\"winnerGate\":";GateJson(out,result.winner);out<<",\"nSigmaGateDiagnostic\":";GateJson(out,result.nsigma);out<<",\"fixedResponseNodes\":[";
  for(unsigned j=0;j<4;++j){if(j)out<<',';out<<"{\"species\":"<<RestartPid::Quote(Species(j))<<",\"supportProvided\":"<<(shapes[j].supportProvided?"true":"false")<<",\"nodes\":[";
    for(unsigned q=0;q<shapes[j].ratios.size();++q){if(q)out<<',';const NuclearFit::Node& n=shapes[j].component.nodes[q];out<<"{\"weight\":"<<n.weight<<",\"mean\":"<<n.mean<<",\"sigma\":"<<n.sigma<<",\"widthRatio\":"<<shapes[j].ratios[q]<<'}';}out<<"]}";}
  out<<"],\"uncertaintyScope\":\"conditional_fixed_curves_and_GL3_R_mixture; local_boundary_approximation; excludes_PID_curve_and_model_uncertainty\"}";return out.str();
}
inline void Run(const char* cache,const char* settingsPath,const char* slicesPath,const char* nodesPath,const char* outputDir){
  RestartPid::Check(cache&&*cache&&settingsPath&&*settingsPath&&slicesPath&&*slicesPath&&nodesPath&&*nodesPath&&outputDir&&*outputDir,"Five explicit paths required");
  const NuclearFit::Settings cfg=NuclearFit::ReadSettings(settingsPath);RestartAllFixed::ValidateSettings(cfg);std::string mapping=RestartStorage::Mapping(cfg);
  const bool transferDiagnostic=NuclearFit::Integer(cfg,"tofTransferDiagnostic")!=0;
  const std::vector<RestartPid::Slice> slices=RestartPid::ReadSlices(slicesPath);const std::map<std::string,Shapes> shapes=ReadNodes(nodesPath);RestartPid::Check(shapes.size()==slices.size(),"Node IDs differ");
  std::string rootPath=std::string(outputDir)+"/native_spectra.root",jsonPath=std::string(outputDir)+"/fits.json",mapPath=std::string(outputDir)+"/native_effective_bins.tsv";
  RestartPid::Check(gSystem->AccessPathName(rootPath.c_str())&&gSystem->AccessPathName(jsonPath.c_str())&&gSystem->AccessPathName(mapPath.c_str()),"Fresh outputs required");
  TFile input(cache,"READ");RestartPid::Check(!input.IsZombie()&&!input.IsWritable(),"Cannot open immutable cache");TFile output(rootPath.c_str(),"CREATE");RestartPid::Check(!output.IsZombie(),"Cannot create fit ROOT");
  RestartStorage::FreshText jsonFile(jsonPath),mapFile(mapPath);mapFile.Write(RestartStorage::MapHeader()+"\n");
  std::vector<std::string> rows,mapRows;std::vector<bool> fitted;unsigned failures=0,accepted=0;
  for(size_t index=0;index<slices.size();++index){const RestartPid::Slice& slice=slices[index];std::map<std::string,Shapes>::const_iterator found=shapes.find(slice.id);RestartPid::Check(found!=shapes.end(),"Missing nodes "+slice.id);
    TH2* source=RestartPid::InputHistogram(input,cfg,slice.mode);RestartPid::Check(slice.rBin<=source->GetNbinsX(),"R bin out of range");double rLo=source->GetXaxis()->GetBinLowEdge(slice.rBin),rHi=source->GetXaxis()->GetBinUpEdge(slice.rBin),split=NuclearFit::Number(cfg,"tofTransitionR");
    RestartPid::Check(transferDiagnostic?slice.mode!="ALL":((rHi<=split+1.e-10&&slice.mode=="ALL")||(rLo>=split-1.e-10&&slice.mode!="ALL")),"R sample transition mismatch");
    TH1D* native=source->ProjectionY(RestartPid::HistogramName("raw_native",index).c_str(),slice.rBin,slice.rBin,"e");RestartPid::Check(native!=0,"Projection failed");native->SetDirectory(0);output.cd();native->Write();bool complete=false;std::string row;
    try{NuclearFit::Request request=Request(*native,slice,found->second,cfg);std::vector<double> edges=RestartStorage::EffectiveEdges(request,mapping);
      for(unsigned b=0;b<request.bins.size();++b){std::string line=RestartStorage::MapRow(slice,*native,request.bins[b],edges[b],edges[b+1],mapping);mapFile.Write(line+"\n");mapRows.push_back(line);}
      RestartStorage::ApplyEdges(request,edges);RestartAllFixed::Result result=RestartAllFixed::Fit(request,found->second,SpeciesIndex(slice.target),cfg);
      if(transferDiagnostic){result.purityAvailable=false;result.reasons.push_back("calibration_transfer_diagnostic_only");}
      RestartPid::SavePrediction(output,*native,slice,result.fit,index,cfg);row=ResultJson(slice,found->second,request,result,*native,index,rLo,rHi,mapping);complete=true;accepted+=result.purityAvailable;
      std::cout<<"ALL_RESTART_FIT id="<<slice.id<<" status="<<result.fit.status<<" purityAvailable="<<result.purityAvailable<<" pearson="<<result.fit.pearson<<" ndf="<<result.fit.conditionalNdf<<std::endl;
    }catch(const std::exception& e){++failures;std::ostringstream error;error<<"{\"id\":"<<RestartPid::Quote(slice.id)<<",\"mode\":"<<RestartPid::Quote(slice.mode)<<",\"target\":"<<RestartPid::Quote(slice.target)<<",\"rBin\":"<<slice.rBin<<",\"rLo\":"<<std::setprecision(17)<<rLo<<",\"rHi\":"<<rHi<<",\"fitLo\":"<<slice.fitLo<<",\"fitHi\":"<<slice.fitHi<<",\"fitCompleted\":false,\"purityAvailable\":false,\"purity\":null,\"nativeHistogram\":"<<RestartPid::Quote(RestartPid::HistogramName("raw_native",index))<<",\"error\":"<<RestartPid::Quote(e.what())<<'}';row=error.str();std::cerr<<"ALL_RESTART_FIT_FAILED id="<<slice.id<<" error="<<e.what()<<std::endl;}
    rows.push_back(row);fitted.push_back(complete);delete native;
  }
  mapFile.Close();output.Write();RestartPid::Check(!output.TestBit(TFile::kWriteError),"ROOT write error");output.Close();RestartPid::Readback(rootPath.c_str(),input,slices,cfg,fitted);
  std::ifstream mapRead(mapPath.c_str());std::string line;RestartPid::Check(bool(std::getline(mapRead,line))&&line==RestartStorage::MapHeader(),"Mapping readback header");
  for(unsigned i=0;i<mapRows.size();++i)RestartPid::Check(bool(std::getline(mapRead,line))&&line==mapRows[i],"Mapping readback mismatch");RestartPid::Check(!std::getline(mapRead,line)&&mapRead.eof(),"Unexpected mapping rows");input.Close();
  std::ostringstream json;json<<"{\"status\":\"completed\",\"study\":\"fixed_shared_nuclear_curves_20261020\",\"productionEnabled\":false,\"curveAdoptionAllowed\":false,\"model\":\"fixed_splitnormal_nonnegative_amplitudes\",\"gatePolicy\":\"strict_nsigma_and_best_species_order_d_t_he3_he4\",\"RIntegration\":\"supplied_GL3_conditional_approximation\",\"storageMapping\":"<<RestartPid::Quote(mapping)<<",\"onFileEncodingAsserted\":false,\"rawCountsModified\":false,\"rawReadback\":true,\"predictionReadback\":true,\"mappingReadback\":true,\"totalSlices\":"<<slices.size()<<",\"fitFailures\":"<<failures<<",\"qualityAcceptedSlices\":"<<accepted<<",\"cachePath\":"<<RestartPid::Quote(cache)<<",\"results\":[";
  for(unsigned i=0;i<rows.size();++i){if(i)json<<',';json<<rows[i];}json<<"]}\n";jsonFile.Write(json.str());jsonFile.Close();std::ifstream reread(jsonPath.c_str());std::ostringstream bytes;bytes<<reread.rdbuf();RestartPid::Check(!reread.bad()&&bytes.str()==json.str(),"JSON exact readback mismatch");
  std::cout<<"ALL_RESTART_FITS_COMPLETE slices="<<slices.size()<<" failures="<<failures<<" accepted="<<accepted<<" rawReadback=1 predictionReadback=1 mappingReadback=1"<<std::endl;
}
}
#endif
void FitAllRestartSlices(const char* cache,const char* settings,const char* slices,const char* nodes,const char* outputDir){
#ifndef __CINT__
  try{RestartAllIO::Run(cache,settings,slices,nodes,outputDir);}catch(const std::exception& e){std::cerr<<"ALL_RESTART_FATAL: "<<e.what()<<std::endl;gSystem->Exit(1);}
#endif
}
