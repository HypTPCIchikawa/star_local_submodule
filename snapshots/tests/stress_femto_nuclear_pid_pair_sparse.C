int stress_femto_nuclear_pid_pair_sparse(const char*,const char*,long long);
#ifndef __CINT__
#include "StMaker/common/FemtoNuclearPairPidSparse.h"
#include "yaml-cpp/yaml.h"
#include "TAxis.h"
#include "TFile.h"
#include "TNamed.h"
#include "TSystem.h"
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace SparseCapacityTest {
unsigned long long checks=0;
void Need(bool ok,const char* text) { ++checks;if(!ok)throw std::runtime_error(text); }
unsigned long long Audit(THnSparseD& h,long long count) {
  const long long extra=count?6:0,expected=count+extra;
  Need(h.GetEntries()==expected,"entries mismatch");
  Need(h.GetNbins()==expected,"unique allocated bin count mismatch");
  std::vector<bool> seen(static_cast<size_t>(expected),false);
  const int nk=h.GetAxis(0)->GetNbins(),nr=h.GetAxis(1)->GetNbins();
  unsigned long long flows=0;
  for(Long64_t b=0;b<h.GetNbins();++b) {
    Int_t c[8];
    Need(h.GetBinContent(b,c)==1.,"unit content mismatch");
    Need(h.GetBinError(b)==1.,"unit error/variance mismatch");
    Need(c[3]==1 && c[4]==1 && c[5]==1 && c[6]==2 && c[7]==2,"fixed category mismatch");
    long long id=-1;int flowAxis=-1,flowSide=-1;
    for(int a=0;a<3;++a) {
      if(c[a]==0 || c[a]==h.GetAxis(a)->GetNbins()+1) {
        Need(flowAxis<0,"unexpected multidimensional flow tuple");
        flowAxis=a;flowSide=c[a]==0?0:1;
      } else Need(c[a]>=1 && c[a]<=h.GetAxis(a)->GetNbins(),"invalid numeric coordinate");
    }
    if(flowAxis>=0) {
      for(int a=0;a<3;++a)if(a!=flowAxis)Need(c[a]==1,"unexpected flow companion coordinate");
      id=count+2*flowAxis+flowSide;++flows;
    } else id=(c[0]-1)+static_cast<long long>(nk)*(c[1]-1)+static_cast<long long>(nk)*nr*(c[2]-1);
    Need(id>=0 && id<expected,"grid cell outside intended coverage");
    Need(!seen[static_cast<size_t>(id)],"duplicate grid cell");seen[static_cast<size_t>(id)]=true;
  }
  for(size_t i=0;i<seen.size();++i) Need(seen[i],"missing intended grid cell");
  Need(flows==static_cast<unsigned long long>(extra),"flow count mismatch");
  return flows;
}
}

int stress_femto_nuclear_pid_pair_sparse(const char* mainconf,const char* output,long long count) {
  using namespace SparseCapacityTest;
  try {
    checks=0;Need(mainconf && *mainconf && output && *output && count>=0,"explicit valid arguments required");
    Need(gSystem->AccessPathName(output),"refuse preexisting sample directory");
    Need(gSystem->mkdir(output,true)==0,"cannot create sample directory");
    const std::string mainPath=mainconf;
    const size_t at=mainPath.rfind("/config/");Need(at!=std::string::npos,"mainconf must resolve below config");
    const YAML::Node main=YAML::LoadFile(mainPath);
    Need(bool(main["nuclearPidPairQaHist"]),"explicit QA YAML required");
    std::string yaml=main["nuclearPidPairQaHist"].as<std::string>();
    if(yaml.empty() || yaml[0]!='/')yaml=mainPath.substr(0,at)+"/config/"+yaml;
    FemtoNuclearPairPidSparse helper;Need(helper.Init(yaml.c_str(),std::cerr),"helper Init failed");
    THnSparseD& h=*helper.Histogram();
    const long long nk=h.GetAxis(0)->GetNbins(),nr=h.GetAxis(1)->GetNbins(),ns=h.GetAxis(2)->GetNbins();
    Need(count<=nk*nr*ns,"requested unique occupancy exceeds configured grid");
    const std::chrono::steady_clock::time_point started=std::chrono::steady_clock::now();
    for(long long i=0;i<count;++i)
      helper.Fill(h.GetAxis(0)->GetBinCenter(i%nk+1),h.GetAxis(1)->GetBinCenter((i/nk)%nr+1),
                  h.GetAxis(2)->GetBinCenter(i/(nk*nr)+1),0,0,false,1,true);
    if(count)for(int axis=0;axis<3;++axis)for(int side=0;side<2;++side) {
      double v[]={h.GetAxis(0)->GetBinCenter(1),h.GetAxis(1)->GetBinCenter(1),h.GetAxis(2)->GetBinCenter(1)};
      v[axis]=side?h.GetAxis(axis)->GetXmax()+h.GetAxis(axis)->GetBinWidth(1):
                        h.GetAxis(axis)->GetXmin()-h.GetAxis(axis)->GetBinWidth(1);
      helper.Fill(v[0],v[1],v[2],0,0,false,1,true);
    }
    const double fillSeconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();
    const unsigned long long flows=Audit(h,count);
    const std::string rootPath=std::string(output)+"/capacity.root";
    {TFile file(rootPath.c_str(),"CREATE");Need(!file.IsZombie(),"ROOT CREATE failed");
      Need(helper.Write(&file),"helper ROOT Write failed");file.Close();}
    {TFile file(rootPath.c_str(),"READ");Need(!file.IsZombie(),"closed ROOT READ failed");
      THnSparseD* read=dynamic_cast<THnSparseD*>(file.Get("hNuclearPidPairJoint"));Need(read!=0,"missing stored joint");
      Need(Audit(*read,count)==flows,"readback flow mismatch");
      TNamed* metadata=dynamic_cast<TNamed*>(file.Get("NuclearPidPairJointMetadata"));Need(metadata!=0,"missing stored metadata");
      std::ostringstream accepted;accepted<<"accepted: "<<count+(count?6:0)<<"\n";
      const std::string text=metadata->GetTitle();Need(text.find(accepted.str())!=std::string::npos,"metadata accepted count mismatch");
      Need(text.find("rejectedNonfinite: 0")!=std::string::npos && text.find("rejectedCategory: 0")!=std::string::npos,"metadata rejects mismatch");}
    const double totalSeconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();
    std::ofstream summary((std::string(output)+"/summary.json").c_str());
    summary.precision(17);
    summary << "{\"passed\":true,\"uniqueRegularBins\":"<<count<<",\"flowBins\":"<<flows
      <<",\"entries\":"<<count+(count?6:0)<<",\"allocatedBins\":"<<h.GetNbins()
      <<",\"checks\":"<<checks<<",\"fillSeconds\":"<<fillSeconds<<",\"fillAuditWriteReadbackSeconds\":"<<totalSeconds
      <<",\"physicalDistribution\":false,\"category\":\"prePID/he4/SE/cent0/signalLow\"}\n";
    summary.close();Need(bool(summary),"cannot write stress summary");
    std::cout<<"PAIR_PID_SPARSE_STRESS_PASS unique="<<count<<" checks="<<checks<<std::endl;return 0;
  } catch(const std::exception& e) {std::cerr<<"PAIR_PID_SPARSE_STRESS_FAIL "<<e.what()<<std::endl;return 1;}
}
#endif
