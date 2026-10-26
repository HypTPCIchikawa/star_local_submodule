// Read-only single-run pilot audit; no physics cuts are changed.
#include "auditFemtoLambdaOutputs.C"

namespace FemtoPilot {
void Count(TDirectory* d, int& hist, int& mt) {
  TIter next(d->GetListOfKeys()); TKey* key=0;
  while ((key=(TKey*)next())) {
    TObject* obj=key->ReadObj();
    FemtoLambdaAudit::Require(obj!=0,"unreadable ROOT key");
    if(obj->InheritsFrom(TDirectory::Class())) Count((TDirectory*)obj,hist,mt);
    else {
      if(obj->InheritsFrom(TH1::Class())) {
        ++hist;
        if(std::string(key->GetName()).find("hKstarMt_")==0) {
          ++mt;
          FemtoLambdaAudit::Require(obj->InheritsFrom(TH2::Class()),"mT map is not TH2");
        }
      }
      delete obj;
    }
  }
}
}

void check_femto_lambda_pilot(const char* inputList, const char* reportPath,
                            int expectedHistograms, int expectedMtMaps,
                            int minRun, int maxRun) {
  using namespace FemtoLambdaAudit;
  try {
    std::ifstream list(inputList);
    Require(list.good(),"cannot open pilot ROOT list");
    Require(gSystem->AccessPathName(reportPath),"refusing existing pilot report");
    std::ofstream report(reportPath);
    Require(report.good(),"cannot create pilot report");
    report<<"file,input_events,reconstructed_events,selected_lambda,histograms,mt_maps,signal_SE,signal_ME,mass_peak_bin\n";
    std::string path; int files=0;
    while(std::getline(list,path)) {
      Require(!path.empty(),"empty input path");
      TFile f(path.c_str(),"READ");
      Require(!f.IsZombie()&&!f.TestBit(TFile::kRecovered),"invalid ROOT file");
      TNamed* status=dynamic_cast<TNamed*>(f.Get("FemtoLambdaRunStatus"));
      Require(status&&std::string(status->GetTitle())=="completed","missing/not completed status");
      std::vector<FemtoLambdaAudit::Event> events=Events(f,true);
      Long64_t reco=0,selected=0;
      for(size_t i=0;i<events.size();++i) {
        Require(events[i].run>=minRun&&events[i].run<=maxRun,"wrong run period");
        if(events[i].status==4)++reco;
        selected+=events[i].selected;
      }
      Require(Parameter(f,"inputEvents")==Long64_t(events.size()),"input/event ledger mismatch");
      Require(Parameter(f,"reconstructedEvents")==reco,"reconstructed counter mismatch");
      Require(Parameter(f,"selectedLambdaCandidates")==selected,"selected counter mismatch");
      std::vector<FemtoLambdaAudit::Candidate> candidates=Candidates(f,false,&events);
      Require(Long64_t(candidates.size())==selected,"candidate ledger mismatch");
      int hist=0,mt=0; FemtoPilot::Count(&f,hist,mt);
      Require(hist==expectedHistograms&&mt==expectedMtMaps,"histogram inventory mismatch");
      const char* acc[]={"hLambda_PtVsYLab_signal","hLambdaProton_PtVsYLab_signal",
                        "hLambdaPion_PtVsYLab_signal","hNucleus_PtVsYLab_he4"};
      for(int a=0;a<4;++a) Require(Hist(f,acc[a])->GetDimension()==2,"acceptance is not TH2");
      const char* modes[]={"SE","ME"};
      for(int mode=0;mode<2;++mode) {
        const std::string base=std::string("hKstar")+modes[mode]+"_lambda_he4_signal";
        TH1 *all=Hist(f,base),*low=Hist(f,base+"Low"),*high=Hist(f,base+"High");
        for(int b=0;b<CellCount(*all);++b)
          Require(all->GetBinContent(b)==low->GetBinContent(b)+high->GetBinContent(b),
                  "low/high signal closure failure");
        Hist(f,std::string("hKstar")+modes[mode]+"VsCent_lambda_he4_signalLow");
        Hist(f,std::string("hKstar")+modes[mode]+"VsCent_lambda_he4_signalHigh");
      }
      TH1* mass=Hist(f,"hLambda_InvMass");
      Require(mass->GetEntries()==selected,"mass histogram/candidate count mismatch");
      report<<path<<","<<events.size()<<","<<reco<<","<<selected<<","<<hist<<","<<mt<<","
            <<Hist(f,"hKstarSE_lambda_he4_signal")->GetEntries()<<","
            <<Hist(f,"hKstarME_lambda_he4_signal")->GetEntries()<<","
            <<mass->GetXaxis()->GetBinCenter(mass->GetMaximumBin())<<"\n";
      ++files;
    }
    Require(files>0,"no pilot files checked");
    report.close();
    Require(report.good(),"pilot report write failed");
    std::cout<<"PASS pilot files="<<files<<": status, event/candidate ledgers, run period, "
             <<"histogram inventory, acceptance, mT and signal-half closure"<<std::endl;
    gSystem->Exit(0);
  } catch(const std::exception& e) {
    std::cerr<<"FAIL pilot: "<<e.what()<<std::endl; gSystem->Exit(1);
  }
}
