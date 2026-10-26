#include <TFile.h>
#include <TKey.h>
#include <TSystem.h>
#include <TObjArray.h>
#include <iostream>
#include <fstream>
#include <map>
#include <string>
void key_sizes(const char* list){
 std::ifstream in(list);std::string path;
 while(std::getline(in,path)){
  TFile f(path.c_str(),"READ");if(f.IsZombie()){std::cout<<"FAILED "<<path<<std::endl;continue;}
  std::cout<<"FILE "<<path<<std::endl;
  std::map<std::string,int> classes;int large=0;TIter next(f.GetListOfKeys());TKey* k=0;
  while((k=(TKey*)next())){
   classes[k->GetClassName()]++;
   if(k->GetObjlen()>100000000){std::cout<<"LARGE "<<k->GetName()<<";"<<k->GetCycle()<<" "<<k->GetClassName()<<" uncompressed="<<k->GetObjlen()<<" compressed="<<k->GetNbytes()<<std::endl;++large;}
  }
  std::cout<<"KEYCOUNT "<<f.GetListOfKeys()->GetSize()<<" LARGECOUNT "<<large<<std::endl;
  for(std::map<std::string,int>::const_iterator it=classes.begin();it!=classes.end();++it)std::cout<<"CLASS "<<it->first<<" "<<it->second<<std::endl;
 }
 gSystem->Exit(0);
}
