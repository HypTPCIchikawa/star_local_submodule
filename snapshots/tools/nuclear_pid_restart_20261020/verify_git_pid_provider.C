// ROOT5 ACLiC: compare the actual mainconf provider with an independent oracle.
// No event loop, histogram fill, library rebuild or production fit is performed.
int verify_git_pid_provider(const char*,const char*,double,double);
#ifndef __CINT__
#include "StMaker/common/FemtoLambdaLegacy.h"
#include "StMaker/common/NuclearIdDeDxVsMom.h"
#include "TMath.h"
#include "TSystem.h"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
void RestartRequire(bool condition,const std::string& message) {
  if(!condition)throw std::runtime_error(message);
}
std::vector<std::string> RestartFields(const std::string& line) {
  std::vector<std::string> result;std::string item;std::istringstream input(line);
  while(std::getline(input,item,'\t'))result.push_back(item);
  return result;
}
double RestartNumber(const std::string& text) {
  std::istringstream input(text);double value=0.;char extra;
  RestartRequire(bool(input>>value) && !(input>>extra) && TMath::Finite(value),
                 "Oracle contains an invalid numerical field: "+text);
  return value;
}
bool RestartNear(double actual,double expected,double relative,double absolute) {
  return TMath::Finite(actual) && TMath::Finite(expected) &&
    std::fabs(actual-expected)<=absolute+relative*std::fabs(expected);
}
std::string RestartQuote(const std::string& value) {
  std::string result="\"";
  for(size_t i=0;i<value.size();++i) {
    if(value[i]=='\\' || value[i]=='\"')result+='\\';
    RestartRequire(value[i]!='\n' && value[i]!='\r',"Unexpected path newline");
    result+=value[i];
  }
  return result+'\"';
}
}

int verify_git_pid_provider(const char* oracle,const char* report,
                            double relative,double absolute) {
  FemtoLambdaLegacy* helper=0;
  try {
    RestartRequire(oracle && oracle[0] && report && report[0],
                   "Explicit oracle and report paths are required");
    RestartRequire(TMath::Finite(relative) && relative>=0. &&
                   TMath::Finite(absolute) && absolute>=0.,
                   "Explicit finite nonnegative tolerances are required");
    RestartRequire(gSystem->AccessPathName(report),"Runtime report must be new");
    std::ifstream input(oracle);RestartRequire(input.good(),"Cannot read oracle TSV");
    std::string line;
    RestartRequire(bool(std::getline(input,line)) &&
                   line=="mainconf\tspecies\trigidity\tdedx\texpectedPull",
                   "Oracle TSV header mismatch");
    unsigned rows=0,checks=0;std::string active;std::set<std::string> mains;
    while(std::getline(input,line)) {
      RestartRequire(!line.empty(),"Empty oracle row is not allowed");
      const std::vector<std::string> fields=RestartFields(line);
      RestartRequire(fields.size()==5,"Oracle TSV must have exactly five columns");
      if(fields[0]!=active) {
        RestartRequire(mains.find(fields[0])==mains.end(),"Mainconf rows must be contiguous");
        delete helper;helper=new FemtoLambdaLegacy();
        RestartRequire(helper->Init(fields[0].c_str(),std::cerr),
                       "Actual provider Init failed: "+fields[0]);
        active=fields[0];mains.insert(active);
      }
      const double rawSpecies=RestartNumber(fields[1]);
      const int species=static_cast<int>(rawSpecies);
      RestartRequire(rawSpecies==species && species>=0 && species<4,"Invalid species index");
      const double rigidity=RestartNumber(fields[2]),dedx=RestartNumber(fields[3]);
      const double expected=RestartNumber(fields[4]);
      RestartRequire(rigidity>0. && dedx>0.,"Oracle must use positive rigidity and dE/dx");
      const NuclearIdDeDxVsMom::SpeciesIndex sp=
        static_cast<NuclearIdDeDxVsMom::SpeciesIndex>(species);
      const double legacy=NuclearIdDeDxVsMom::GetNSigma(sp,rigidity,dedx);
      const double actual=helper->NuclearDedxPull(species,rigidity,dedx);
      RestartRequire(RestartNear(legacy,expected,relative,absolute),
                     "Git header and independent Python oracle disagree: "+line);
      ++checks;
      RestartRequire(RestartNear(actual,expected,relative,absolute),
                     "Actual mainconf provider differs from Git PID: "+line);
      ++checks;
      RestartRequire(RestartNear(actual,legacy,relative,absolute),
                     "Actual mainconf provider differs from the C++ legacy function: "+line);
      ++checks;++rows;
    }
    RestartRequire(input.eof() && rows>0 && !mains.empty(),"Incomplete or empty oracle input");
    delete helper;helper=0;
    std::ofstream output(report);RestartRequire(output.good(),"Cannot write runtime report");
    output<<std::setprecision(17)
          <<"{\n  \"status\": \"passed\",\n  \"runtimeProviderVerified\": true,"
          <<"\n  \"noEventLoop\": true,\n  \"productionFitsExecuted\": false,"
          <<"\n  \"oracle\": "<<RestartQuote(oracle)
          <<",\n  \"mainconfCount\": "<<mains.size()
          <<",\n  \"oracleRows\": "<<rows<<",\n  \"comparisonChecks\": "<<checks
          <<",\n  \"comparisonTolerance\": {\"relative\": "<<relative
          <<", \"absolute\": "<<absolute<<"}\n}\n";
    output.close();RestartRequire(output.good(),"Runtime report write failed");
    std::cout<<"GIT_PID_PROVIDER_PASS mainconfs="<<mains.size()
             <<" rows="<<rows<<" checks="<<checks<<std::endl;
    return 0;
  } catch(const std::exception& error) {
    delete helper;
    std::cerr<<"GIT_PID_PROVIDER_FAIL "<<error.what()<<std::endl;
    return 1;
  }
}
#endif
