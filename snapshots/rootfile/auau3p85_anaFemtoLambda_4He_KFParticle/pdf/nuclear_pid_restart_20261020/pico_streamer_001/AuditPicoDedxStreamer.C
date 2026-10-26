// Read-only, metadata-only audit of one cached 2021 P24iy PicoDst.
#include <TFile.h>
#include <TList.h>
#include <TObjArray.h>
#include <TStreamerInfo.h>
#include <TStreamerElement.h>
#include <TSystem.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

void AuditPicoDedxStreamer(const char* inputPath, const char* reportPath) {
  if (!inputPath || !inputPath[0] || !reportPath || !reportPath[0] ||
      !gSystem->AccessPathName(reportPath)) {
    std::cerr << "Explicit input and NEW report path required" << std::endl;
    gSystem->Exit(20); return;
  }
  TFile input(inputPath, "READ");
  if (input.IsZombie() || input.TestBit(TFile::kRecovered)) {
    std::cerr << "Missing, zombie or recovered input; no metadata adopted" << std::endl;
    gSystem->Exit(21); return;
  }
  TList* infos = input.GetStreamerInfoList();
  if (!infos) { gSystem->Exit(22); return; }
  std::ofstream output(reportPath);
  if (!output.good()) { gSystem->Exit(23); return; }
  output << std::setprecision(17);
  output << "{\"metadataOnly\":true,\"eventsRead\":0,\"input\":\""
         << inputPath << "\",\"fileWrittenRootVersion\":" << input.GetVersion()
         << ",\"streamers\":[";
  int count = 0;
  TIter next(infos);
  TObject* object = 0;
  while ((object = next())) {
    TStreamerInfo* info = dynamic_cast<TStreamerInfo*>(object);
    if (!info || std::string(info->GetName()) != "StPicoTrack") continue;
    TObjArray* elements = info->GetElements();
    if (!elements) continue;
    TStreamerElement* element = dynamic_cast<TStreamerElement*>(elements->FindObject("mDedx"));
    if (!element) continue;
    if (count) output << ',';
    output << "{\"classVersion\":" << info->GetClassVersion()
           << ",\"checkSum\":" << info->GetCheckSum()
           << ",\"member\":\"mDedx\",\"typeName\":\"" << element->GetTypeName()
           << "\",\"typeCode\":" << element->GetType()
           << ",\"factor\":" << element->GetFactor()
           << ",\"xmin\":" << element->GetXmin()
           << ",\"xmax\":" << element->GetXmax()
           << ",\"arrayDim\":" << element->GetArrayDim() << '}';
    ++count;
  }
  output << "],\"matchingStreamers\":" << count << "}\n";
  output.close();
  input.Close();
  if (count < 1 || !output.good()) { gSystem->Exit(24); return; }
  std::cout << "PICO_DEDX_STREAMER_AUDIT_PASS matchingStreamers=" << count
            << " eventsRead=0" << std::endl;
}
