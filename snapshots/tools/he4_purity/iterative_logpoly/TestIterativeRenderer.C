// Synthetic numeric/display fixtures ONLY. No real event data or PID settings.
#include "RenderIterativeHelium.C"

#ifndef __CINT__
namespace IterativeRendererTest {
using namespace IterativeHeliumViews;
double GaussianBin(double a,double mu,double sigma,double low,double high) {
  return a*.5*(TMath::Erf((high-mu)/(std::sqrt(2.)*sigma))-TMath::Erf((low-mu)/(std::sqrt(2.)*sigma)));
}
Fields SettingsFixture() {
  Fields c;
  const char* keys[]={"schemaVersion","panelsPerPage","expectedSlices","dedxRebinFactor","displayNSigmaMin","displayNSigmaMax","diagnosticRigidityMin","diagnosticRigidityMax","axisRelativeTolerance","componentRelativeTolerance","predictionNonnegativeTolerance","linearMaximumFactor","logMaximumFactor","logMinimum","linearCanvasWidth","linearCanvasHeight","diagnosticCanvasWidth","diagnosticCanvasHeight","paperWidthCm","paperHeightCm","marginLeft","marginRight","marginBottom","marginTop","axisFont","axisTitleSize","axisLabelSize","xTitleOffset","yTitleOffset","markerSize","legendXMin","legendXMax","legendYMin","legendYMax","legendColumns","legendTextSize","headerY","headerTextSize","statusTextSize","statusFirstY","statusSecondY","captionTextSize","captionFirstY","captionSecondY","captionThirdY"};
  const char* values[]={"1","4","2","5","-8","4","0.4","0.5","1e-9","1e-9","1e-9","1.35","8","0.5","1600","1200","1800","1500","30","22","0.18","0.04","0.26","0.33","42","0.030","0.026","1.35","1.8","0.25","0.18","0.96","0.74","0.87","2","0.020","0.966","0.030","0.019","0.916","0.888","0.015","0.090","0.055","0.022"};
  Require(sizeof(keys)==sizeof(values),"Synthetic setting dimensions mismatch");
  for(unsigned i=0;i<sizeof(keys)/sizeof(keys[0]);++i)c[keys[i]]=values[i];
  c["allowSparseSlices"]="0";
  c["nativeHistogramPattern"]="raw_native_%03d";c["rawHistogramPattern"]="raw_%03d";c["nsigmaHistogramPattern"]="projection_%03d";c["fitMaskPattern"]="fit_used_%03d";
  for(int i=0;i<5;++i)c[std::string("prediction_")+components[i]+"Pattern"]=std::string("prediction_")+components[i]+"_%03d";
  return c;
}
Fields ManifestFixture(int index) {
  Fields f;f["index"]=Form("%d",index);f["lo"]=Form("%.17g",.4+.01*index);f["hi"]=Form("%.17g",.41+.01*index);f["r"]=Form("%.17g",.405+.01*index);
  f["anchorMean"]="106.05761652176085";f["anchorSigma"]="6.265485936208043";f["fitLoDedx"]="54.8";f["fitHiDedx"]="119.8";f["firstFitBin"]="111";f["lastFitBin"]="240";
  f["stageMode"]=index?"shared_fixed":"local_free";f["model"]="linear";f["valid"]="1";f["quality"]=index?"-1":"0";f["status"]="0";f["covarianceStatus"]="3";f["atBoundary"]="0";f["chi2"]=index?"-1":"201";f["deviance"]="190";f["penalty"]="2.5";f["ndf"]="120";f["rawLocalMaximum"]=index?"-1":"0";f["oldAccepted"]="0";
  f["tShapeState"]="free";f["he3ShapeState"]=index?"shared_fixed":"free";f["he4ShapeState"]=index?"shared_fixed":"unidentified";return f;
}
void WriteFixtures(const std::string& output,const Fields& c,const std::vector<Fields>& rows) {
  TFile root((output+"/phase.root").c_str(),"NEW");Require(!root.IsZombie(),"Refusing to overwrite synthetic ROOT");
  for(unsigned index=0;index<rows.size();++index) {
    TH1D native(Form("raw_native_%03d",index),"synthetic observed counts",1200,0,120);native.Sumw2();
    for(int b=1;b<=native.GetNbinsX();++b) {
      const double low=native.GetBinLowEdge(b),high=low+native.GetBinWidth(b);
      const double count=GaussianBin(60000,64,4,low,high)+GaussianBin(45000,72,4.5,low,high)+GaussianBin(9000,83,8,low,high)+15;
      native.SetBinContent(b,std::floor(count+.5));native.SetBinError(b,std::sqrt(native.GetBinContent(b)));
    }
    native.SetBinContent(0,11);native.SetBinError(0,std::sqrt(11.));native.SetBinContent(1201,37);native.SetBinError(1201,std::sqrt(37.));native.SetEntries(native.Integral(0,1201));native.Write();
    TH1D* raw=static_cast<TH1D*>(native.Clone(Form("raw_%03d",index)));raw->SetDirectory(0);raw->Rebin(5);raw->Write();
    TH1D nsigma(Form("projection_%03d",index),"synthetic fixed-anchor transform",240,(0-Number(rows[index],"anchorMean"))/Number(rows[index],"anchorSigma"),(120-Number(rows[index],"anchorMean"))/Number(rows[index],"anchorSigma"));nsigma.Sumw2();
    for(int b=0;b<=241;++b){nsigma.SetBinContent(b,raw->GetBinContent(b));nsigma.SetBinError(b,raw->GetBinError(b));}nsigma.SetEntries(raw->GetEntries());nsigma.Write();
    TH1D mask(Form("fit_used_%03d",index),"used synthetic fit bins",240,0,120);TH1D* prediction[5];
    for(int i=0;i<5;++i){prediction[i]=new TH1D(Form("prediction_%s_%03d",components[i],index),"synthetic bin-integrated expectations",240,0,120);prediction[i]->SetDirectory(0);}
    for(int b=111;b<=240;++b) {
      mask.SetBinContent(b,1);const double low=raw->GetBinLowEdge(b),high=low+raw->GetBinWidth(b);
      prediction[1]->SetBinContent(b,75);prediction[2]->SetBinContent(b,GaussianBin(60000,64,4,low,high));prediction[3]->SetBinContent(b,GaussianBin(45000,72,4.5,low,high));prediction[4]->SetBinContent(b,GaussianBin(9000,83,8,low,high));double total=0;
      for(int i=1;i<5;++i)total+=prediction[i]->GetBinContent(b);prediction[0]->SetBinContent(b,total);
    }mask.Write();for(int i=0;i<5;++i){prediction[i]->Write();delete prediction[i];}delete raw;
  }root.Close();
  std::ofstream settings((output+"/settings.tsv").c_str());for(Fields::const_iterator i=c.begin();i!=c.end();++i)settings<<i->first<<'\t'<<i->second<<'\n';settings.close();Require(settings.good(),"Synthetic settings write failed");
  std::ofstream manifest((output+"/manifest.tsv").c_str());for(Fields::const_iterator i=rows[0].begin();i!=rows[0].end();++i){if(i!=rows[0].begin())manifest<<'\t';manifest<<i->first;}manifest<<'\n';
  for(unsigned row=0;row<rows.size();++row){for(Fields::const_iterator i=rows[row].begin();i!=rows[row].end();++i){if(i!=rows[row].begin())manifest<<'\t';manifest<<i->second;}manifest<<'\n';}manifest.close();Require(manifest.good(),"Synthetic manifest write failed");
}
}
#endif

void TestIterativeRenderer(const char* outputDirectory) {
#ifndef __CINT__
  using namespace IterativeRendererTest;
  try {
    Require(!gSystem->AccessPathName(outputDirectory),"Explicit synthetic output directory required");const std::string output(outputDirectory);const Fields c=SettingsFixture();std::vector<Fields> rows;rows.push_back(ManifestFixture(0));rows.push_back(ManifestFixture(1));ValidateSettings(c);WriteFixtures(output,c,rows);
    TFile* source=TFile::Open((output+"/phase.root").c_str(),"READ");Require(source&&!source->IsZombie(),"Cannot read synthetic ROOT");
    Slice s=LoadSlice(source,rows[0],c);Require(s.firstDisplay>=1&&s.lastDisplay<=240,"Display includes flow bins");Require(s.raw->GetBinContent(0)==11&&s.raw->GetBinContent(241)==37,"Synthetic flow preservation failed");int passed=1;
    const double original=s.nsigma->GetBinContent(0);s.nsigma->SetBinContent(0,original+1);bool rejected=false;try{LoadSlice(source,rows[0],c);}catch(const std::exception&){rejected=true;}Require(rejected,"Failed to reject changed transformed underflow");s.nsigma->SetBinContent(0,original);++passed;
    const double total=s.prediction[0]->GetBinContent(140);s.prediction[0]->SetBinContent(140,total+1);rejected=false;try{LoadSlice(source,rows[0],c);}catch(const std::exception&){rejected=true;}Require(rejected,"Failed to reject component closure change");s.prediction[0]->SetBinContent(140,total);++passed;
    s.mask->SetBinContent(150,0);rejected=false;try{LoadSlice(source,rows[0],c);}catch(const std::exception&){rejected=true;}Require(rejected,"Failed to reject changed used-fit mask");s.mask->SetBinContent(150,1);++passed;
    Fields unknown=rows[1];LoadSlice(source,unknown,c);++passed;source->Close();delete source;
    std::ofstream report((output+"/numeric_selftest.json").c_str());report<<"{\"syntheticOnly\":true,\"productionChanged\":false,\"passed\":"<<passed<<",\"checks\":[\"normal-bin clipping, observed count/error/flow and native variance closure\",\"reject affine-transform flow mutation\",\"reject model component nonclosure\",\"reject saved used-fit mask inconsistency\",\"unknown quality/peak/chi2 remains explicit\"]}\n";report.close();Require(report.good(),"Numeric selftest report write failed");
    RenderIterativeHelium((output+"/phase.root").c_str(),(output+"/manifest.tsv").c_str(),(output+"/settings.tsv").c_str(),outputDirectory,"SYNTHETIC TEST (not data)");
    std::cout<<"ITERATIVE_RENDERER_SELFTEST_PASS checks="<<passed<<std::endl;
  }catch(const std::exception& e){std::cerr<<"ITERATIVE_RENDERER_SELFTEST_ERROR "<<e.what()<<std::endl;gSystem->Exit(1);}
#endif
}
