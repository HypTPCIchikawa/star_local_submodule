#include <TFile.h>
#include <THnSparse.h>
#include <TH1D.h>
#include <TSystem.h>
#include <fstream>
#include <string>
void project_joint_parts(const char* listPath,const char* outputPath,int axis,int cent9,int representation,int mixingMode,int pidStage,int massLo,int massHi){
 if(axis<0||axis>2||cent9< -1||cent9>8||representation<0||representation>2||mixingMode<0||mixingMode>1||pidStage<0||pidStage>1||massLo<0||massHi>3||massLo>massHi||!gSystem->AccessPathName(outputPath)){gSystem->Exit(81);return;}
 TH1::AddDirectory(kFALSE);std::ifstream in(listPath);std::string path;TH1D* total=0;int files=0;
 while(std::getline(in,path)){
  if(path.empty())continue;TFile f(path.c_str(),"READ");if(f.IsZombie()||f.TestBit(TFile::kRecovered)){gSystem->Exit(82);return;}
  THnSparse* h=dynamic_cast<THnSparse*>(f.Get("hNuclearPidPairJoint"));if(!h||h->GetNdimensions()!=8){gSystem->Exit(83);return;}
  for(int a=0;a<3;++a)h->GetAxis(a)->SetRange(0,h->GetAxis(a)->GetNbins()+1);
  if(cent9>=0)h->GetAxis(3)->SetRange(cent9+1,cent9+1);
  h->GetAxis(4)->SetRange(representation+1,representation+1);h->GetAxis(5)->SetRange(mixingMode+1,mixingMode+1);
  h->GetAxis(6)->SetRange(massLo+1,massHi+1);h->GetAxis(7)->SetRange(pidStage+1,pidStage+1);
  TH1D* projected=h->Projection(axis,"EA");if(!projected){gSystem->Exit(84);return;}
  projected->SetDirectory(0);if(!total){total=projected;total->SetName("hProjectedNuclearPidPairJoint");}else{total->Add(projected);delete projected;}
  delete h;f.Close();++files;
 }
 if(!files||!total){gSystem->Exit(85);return;}TFile output(outputPath,"CREATE");if(output.IsZombie()||total->Write()<=0){gSystem->Exit(86);return;}output.Close();delete total;gSystem->Exit(0);
}
