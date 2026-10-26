// Actual production-profile checks. Expected accepted ranges arrive separately.
int test_femto_nuclear_pid_actual_domain(const char*,const char*,const char*,const char*);
#ifndef __CINT__
#include "test_femto_nuclear_pid_domain.C"
#include <map>
#include <set>
#include <vector>
#include <cstdlib>
namespace ActualDomain {
struct Range{double lo,hi;Range(double a,double b):lo(a),hi(b){}};
typedef std::map<std::string,std::string> Row;
typedef std::vector<Row> Rows;
std::vector<std::string> Split(const std::string& line){std::vector<std::string> v;std::string cell;std::istringstream s(line);while(std::getline(s,cell,'\t'))v.push_back(cell);return v;}
Rows Read(const char* path){std::ifstream input(path);DomainRequire(bool(input),"expected TSV readable");std::string line;DomainRequire(bool(std::getline(input,line)),"expected header exists");std::vector<std::string> names=Split(line);Rows rows;while(std::getline(input,line)){std::vector<std::string> fields=Split(line);DomainRequire(fields.size()==names.size(),"expected field count");Row row;for(unsigned i=0;i<names.size();++i)DomainRequire(row.insert(std::make_pair(names[i],fields[i])).second,"unique expected column");rows.push_back(row);}DomainRequire(input.eof()&&!rows.empty(),"expected rows present");return rows;}
double Num(const Row& row,const char* key){Row::const_iterator it=row.find(key);DomainRequire(it!=row.end(),"expected numeric column exists");char* end=0;double v=std::strtod(it->second.c_str(),&end);DomainRequire(end&&!*end&&std::isfinite(v),"finite expected numeric field");return v;}
bool Allowed(const std::vector<Range>& ranges,double p){for(unsigned i=0;i<ranges.size();++i)if(p>=ranges[i].lo&&p<ranges[i].hi)return true;return false;}
void Check(FemtoLambdaLegacy& helper,const std::vector<Range> ranges[4],const double lo[4],const double hi[4]){
  // These regression expectations are the authorized contract, not runtime cuts.
  const double threshold=2.,cap=5.,inf=std::numeric_limits<double>::infinity(),nan=std::numeric_limits<double>::quiet_NaN();
  std::set<double> probes;for(int i=1;i<=1020;++i)probes.insert(i*.005);
  probes.insert(0);probes.insert(-1);probes.insert(threshold);probes.insert(cap);
  probes.insert(std::nextafter(threshold,0.));probes.insert(std::nextafter(threshold,inf));
  probes.insert(std::nextafter(cap,0.));probes.insert(std::nextafter(cap,inf));
  for(int s=0;s<4;++s)for(unsigned i=0;i<ranges[s].size();++i){const double edges[]={ranges[s][i].lo,ranges[s][i].hi};for(int k=0;k<2;++k){probes.insert(edges[k]);probes.insert(std::nextafter(edges[k],0.));probes.insert(std::nextafter(edges[k],inf));}}
  for(int species=0;species<4;++species){double pulls[]={999.,999.,999.,999.};pulls[species]=.1;const double inside=.5*(lo[species]+hi[species]);const double values[]={inside,lo[species],hi[species],std::nextafter(lo[species],hi[species]),std::nextafter(hi[species],lo[species]),nan,inf,-inf};
    for(std::set<double>::const_iterator p=probes.begin();p!=probes.end();++p)for(int tof=0;tof<2;++tof)for(unsigned m=0;m<8;++m){bool pass=*p>0&&*p<cap&&Allowed(ranges[species],*p)&&(*p<threshold||(tof&&std::isfinite(values[m])&&values[m]>lo[species]&&values[m]<hi[species]));DomainRequire(helper.SelectNuclearType(*p,pulls,tof!=0,values[m])==(pass?species:-1),"actual accepted-range/strict-TOF/cap agreement");}
    for(int alternate=0;alternate<4;++alternate)if(alternate!=species){pulls[alternate]=.2;for(std::set<double>::const_iterator p=probes.begin();p!=probes.end();++p)if(*p>0&&*p<cap&&!Allowed(ranges[species],*p)&&Allowed(ranges[alternate],*p))DomainRequire(helper.SelectNuclearType(*p,pulls,true,inside)==-1,"rejected winner never promotes alternate species");pulls[alternate]=999.;}
    const double invalid[]={nan,inf,-inf};for(unsigned i=0;i<3;++i)DomainRequire(helper.SelectNuclearType(invalid[i],pulls,true,inside)==-1,"nonfinite actual rigidity fails closed");
  }
  const double ties[]={.1,.1,.1,.1};for(std::set<double>::const_iterator p=probes.begin();p!=probes.end();++p)if(*p>0&&*p<cap)DomainRequire(helper.SelectNuclearType(*p,ties,true,.5*(lo[0]+hi[0]))==(Allowed(ranges[0],*p)?0:-1),"TPC ties choose d before post-winner rejection");
}
}
#endif
int test_femto_nuclear_pid_actual_domain(const char* mainconf,const char* expectedRanges,const char* expectedWindows,const char* outputDirectory){
#ifndef __CINT__
  using namespace ActualDomain;
  try{DomainChecks=0;DomainRequire(mainconf&&*mainconf&&expectedRanges&&*expectedRanges&&expectedWindows&&*expectedWindows&&outputDirectory&&*outputDirectory,"all actual regression inputs explicit");std::string source(mainconf),output(outputDirectory);if(source[0]!='/')source=std::string(gSystem->WorkingDirectory())+"/"+source;if(output[0]!='/')output=std::string(gSystem->WorkingDirectory())+"/"+output;DomainRequire(gSystem->AccessPathName(output.c_str()),"test output is new");DomainRequire(gSystem->mkdir(output.c_str(),kTRUE)==0,"create output");
    Rows accepted=Read(expectedRanges),windows=Read(expectedWindows);std::vector<Range> ranges[4];double lo[4],hi[4];bool found[4]={false,false,false,false};for(unsigned i=0;i<accepted.size();++i){int s=-1;for(int k=0;k<4;++k)if(accepted[i]["species"]==DomainSpecies[k])s=k;DomainRequire(s>=0,"known expected species");double a=Num(accepted[i],"min"),b=Num(accepted[i],"max");DomainRequire(a>0&&b>a,"valid half-open approved range");ranges[s].push_back(Range(a,b));}for(unsigned i=0;i<windows.size();++i)for(int s=0;s<4;++s)if(windows[i]["species"]==DomainSpecies[s]){lo[s]=Num(windows[i],"studiedMin");hi[s]=Num(windows[i],"studiedMax");DomainRequire(!found[s]&&hi[s]>lo[s],"unique expected TOF window");found[s]=true;}
    const size_t marker=source.rfind("/config/");DomainRequire(marker!=std::string::npos,"mainconf beneath config");YAML::Node main=YAML::LoadFile(source);std::string ref=main["nuclearid"].as<std::string>();YAML::Node nuclear=YAML::LoadFile(ref[0]=='/'?ref:source.substr(0,marker)+"/config/"+ref);DomainRequire(nuclear["maxPOverQ"].as<double>()==5.&&!nuclear["highRigidityTof"]["extendHe4BeyondMaxPOverQ"].as<bool>()&&nuclear["highRigidityTof"]["thresholdPOverQ"].as<double>()==2.,"cap5/nonextended/TOF-at2 authorized contract");for(int s=0;s<4;++s){DomainRequire(!ranges[s].empty()&&found[s],"all expected species present");YAML::Node w=nuclear["highRigidityTof"]["species"][DomainSpecies[s]];DomainRequire(w["mass2Min"].as<double>()==lo[s]&&w["mass2Max"].as<double>()==hi[s],"actual TOF cut equals independently saved study window");}
    FemtoLambdaLegacy helper;DomainRequire(helper.Init(source.c_str(),std::cerr),"actual mainconf Init");Check(helper,ranges,lo,hi);
    {TFile file((output+"/actual_profile.root").c_str(),"CREATE");DomainRequire(!file.IsZombie()&&helper.Write(&file),"write actual metadata ROOT");}
    {TFile file((output+"/actual_profile.root").c_str(),"READ");TNamed* meta=dynamic_cast<TNamed*>(file.Get("FemtoNuclearPidDomain"));DomainRequire(meta&&std::string(meta->GetTitle()).find("never reassign")!=std::string::npos,"actual domain metadata readback");}
    std::cout<<"PASS actual nuclear PID domain: "<<DomainChecks<<" checks; "<<source<<std::endl;return 0;
  }catch(const std::exception& e){std::cerr<<"FAIL actual nuclear PID domain after "<<DomainChecks<<" checks: "<<e.what()<<std::endl;return 1;}
#else
  return 1;
#endif
}
