#include "FemtoNuclearDedxCalibration.h"
#include "NuclearIdDeDxVsMom.h"
#include "yaml-cpp/yaml.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>
#include <stdexcept>

namespace {
// Numerical equality tolerance only, not a physics interval or cut value.
const double kBoundaryTolerance = 1.e-12;
template<class T> T Require(const YAML::Node& n, const char* key) {
  if (!n[key]) throw std::runtime_error(std::string("Missing calibration key: ")+key);
  return n[key].as<T>();
}
double Finite(const YAML::Node& n, const char* key) {
  const double v=Require<double>(n,key);
  if (!std::isfinite(v)) throw std::runtime_error(std::string("Nonfinite calibration key: ")+key);
  return v;
}
void KnownKeys(const YAML::Node& n, const char* const* keys, unsigned count, const char* label) {
  if (!n.IsMap()) throw std::runtime_error(std::string(label)+" must be a map");
  for (YAML::const_iterator i=n.begin();i!=n.end();++i) {
    const std::string key=i->first.as<std::string>();
    bool known=false;
    for(unsigned j=0;j<count;++j) known=known || key==keys[j];
    if(!known) throw std::runtime_error(std::string("Unknown ")+label+" key: "+key);
  }
}
std::string Text(const YAML::Node& n, const char* key) {
  const std::string value=Require<std::string>(n,key);
  if(value.empty()) throw std::runtime_error(std::string("Empty calibration key: ")+key);
  return value;
}
void Sha256(const YAML::Node& n,const char* key) {
  const std::string value=Text(n,key);
  if(value.size()!=64) throw std::runtime_error(std::string("Invalid SHA256: ")+key);
  for(size_t i=0;i<value.size();++i)
    if(!((value[i]>='0'&&value[i]<='9')||(value[i]>='a'&&value[i]<='f')||
         (value[i]>='A'&&value[i]<='F')))
      throw std::runtime_error(std::string("Invalid SHA256: ")+key);
}
bool GoodSpecies(int s) { return s>=0 && s<4; }
}

FemtoNuclearDedxCalibration::FemtoNuclearDedxCalibration() { ResetLegacy(); }
void FemtoNuclearDedxCalibration::ResetLegacy() {
  mEnabled=false;mMaxInterpolationGap=0;
  mCurve[0]=Curve();mCurve[1]=Curve();
  mMetadata="mode: legacy\nsource: StMaker/common/NuclearIdDeDxVsMom.h\n"
            "dedxCalibration: absent or explicitly disabled\n";
}

bool FemtoNuclearDedxCalibration::Load(const std::string& path,std::ostream& diagnostics) {
  ResetLegacy();
  try {
    std::ifstream input(path.c_str());
    if(!input) throw std::runtime_error("Cannot read calibration: "+path);
    const std::string contents((std::istreambuf_iterator<char>(input)),std::istreambuf_iterator<char>());
    if(input.bad() || contents.empty()) throw std::runtime_error("Empty/unreadable calibration: "+path);
    const YAML::Node doc=YAML::Load(contents);
    const char* top[]={"schemaVersion","profile","rigidityUnit","dedxUnit","interpolation",
      "edgePolicy","outsideRangePolicy","maxInterpolationGap","source","species"};
    KnownKeys(doc,top,sizeof(top)/sizeof(top[0]),"dedx calibration");
    if(Require<int>(doc,"schemaVersion")!=1) throw std::runtime_error("Unsupported dedx schemaVersion");
    const std::string profile=Text(doc,"profile");
    if(Text(doc,"rigidityUnit")!="GeV/c" || Text(doc,"dedxUnit")!="keV/cm")
      throw std::runtime_error("Unsupported dedx calibration units");
    if(Text(doc,"interpolation")!="linear" || Text(doc,"edgePolicy")!="clampWithinDomain" ||
       Text(doc,"outsideRangePolicy")!="legacy")
      throw std::runtime_error("Unsupported dedx interpolation/edge/outside policy");
    mMaxInterpolationGap=Finite(doc,"maxInterpolationGap");
    if(!(mMaxInterpolationGap>0)) throw std::runtime_error("maxInterpolationGap must be positive");
    const YAML::Node source=doc["source"];
    if(!source || !source.IsMap()) throw std::runtime_error("Calibration source must be a map");
    if(Text(source,"collision")!="auau3p85" || Require<int>(source,"year")!=2021 ||
       Text(source,"production")!="P24iy")
      throw std::runtime_error("Calibration source is not auau3p85/2021/P24iy");
    Text(source,"sourceRoot");
    if(!source["sourceRootFingerprint"] || !source["sourceRootFingerprint"].IsMap() ||
       source["sourceRootFingerprint"].size()==0)
      throw std::runtime_error("Missing sourceRootFingerprint map");
    const char* paths[]={"fitSummary","fitSource","settings"};
    const char* hashes[]={"fitSummarySha256","fitSourceSha256","settingsSha256"};
    for(unsigned i=0;i<3;++i) { Text(source,paths[i]);Sha256(source,hashes[i]); }
    const char* speciesKeys[]={"he3","he4"};
    const YAML::Node species=doc["species"];
    KnownKeys(species,speciesKeys,2,"calibration species");
    for(int s=0;s<2;++s) {
      const YAML::Node node=species[speciesKeys[s]];
      const char* curveKeys[]={"domain","supportBinHalfWidth","knots"};
      KnownKeys(node,curveKeys,3,"species curve");
      const YAML::Node domain=node["domain"];
      if(!domain || !domain.IsSequence() || domain.size()!=2)
        throw std::runtime_error("Species domain must have exactly two edges");
      Curve& c=mCurve[s];
      c.minimum=domain[0].as<double>();c.maximum=domain[1].as<double>();
      c.supportHalfWidth=Finite(node,"supportBinHalfWidth");
      if(!std::isfinite(c.minimum)||!std::isfinite(c.maximum)||
         !(c.minimum>0 && c.maximum>c.minimum && c.supportHalfWidth>0) ||
         !(2*c.supportHalfWidth<=mMaxInterpolationGap+kBoundaryTolerance))
        throw std::runtime_error("Invalid species domain/supportBinHalfWidth");
      const YAML::Node knots=node["knots"];
      if(!knots || !knots.IsSequence() || knots.size()<2)
        throw std::runtime_error("Species calibration needs at least two knots");
      for(size_t i=0;i<knots.size();++i) {
        const YAML::Node n=knots[i];
        const char* knotKeys[]={"rigidity","mean","sigma","provenance"};
        KnownKeys(n,knotKeys,4,"calibration knot");
        Knot k;k.rigidity=Finite(n,"rigidity");k.mean=Finite(n,"mean");k.sigma=Finite(n,"sigma");
        if(!(k.rigidity>=c.minimum-kBoundaryTolerance && k.rigidity<=c.maximum+kBoundaryTolerance &&
             k.mean>0 && k.sigma>0) || (!c.knots.empty() && !(k.rigidity>c.knots.back().rigidity)))
          throw std::runtime_error("Unordered/out-of-domain/nonpositive calibration knot");
        if(n["provenance"] && !n["provenance"].IsMap())
          throw std::runtime_error("Knot provenance must be a map");
        c.knots.push_back(k);
      }
      if(c.knots.front().rigidity-c.minimum>c.supportHalfWidth+kBoundaryTolerance ||
         c.maximum-c.knots.back().rigidity>c.supportHalfWidth+kBoundaryTolerance)
        throw std::runtime_error("Domain extends farther than supported endpoint half-bin");
    }
    mEnabled=true;
    std::ostringstream metadata;
    metadata << "mode: explicit_yaml\nresolvedCalibrationFile: " << path
             << "\nprofile: " << profile << "\nlegacySpecies: deuteron,triton\n"
             << "gapPolicy: supported endpoint half-bin clamp; unsupported gap uses legacy\n"
             << "----- exact consumed calibration YAML -----\n" << contents;
    mMetadata=metadata.str();
    diagnostics << "[FemtoNuclearDedxCalibration] calibration(argument)=" << path
                << " profile=" << profile << "\n  interpolation=linear; outsideRangePolicy=legacy; "
                << "edgePolicy=clampWithinDomain; maxInterpolationGap=" << mMaxInterpolationGap
                << "\n  he3 knots=" << mCurve[0].knots.size() << "; he4 knots=" << mCurve[1].knots.size()
                << "; d/t retain legacy tables" << std::endl;
    return true;
  } catch(const std::exception& error) {
    ResetLegacy();
    diagnostics << "[FemtoNuclearDedxCalibration] Load failed: " << error.what() << std::endl;
    return false;
  }
}

bool FemtoNuclearDedxCalibration::Evaluate(int species,double r,double& mean,double& sigma) const {
  if(!mEnabled || species<2 || species>3 || !std::isfinite(r) || r<=0) return false;
  const Curve& c=mCurve[species-2];
  if(r<c.minimum-kBoundaryTolerance || r>c.maximum+kBoundaryTolerance) return false;
  // Lower-bound lookup is logarithmic in the number of exported fit points.
  size_t first=0,last=c.knots.size();
  while(first<last) {
    const size_t middle=first+(last-first)/2;
    if(c.knots[middle].rigidity<r) first=middle+1;else last=middle;
  }
  const size_t i=first;
  if(i<c.knots.size() && std::fabs(c.knots[i].rigidity-r)<=kBoundaryTolerance) {
    mean=c.knots[i].mean;sigma=c.knots[i].sigma;return true;
  }
  if(i>0 && std::fabs(c.knots[i-1].rigidity-r)<=kBoundaryTolerance) {
    mean=c.knots[i-1].mean;sigma=c.knots[i-1].sigma;return true;
  }
  if(i==0) { mean=c.knots.front().mean;sigma=c.knots.front().sigma;return true; }
  if(i==c.knots.size()) { mean=c.knots.back().mean;sigma=c.knots.back().sigma;return true; }
  const Knot& right=c.knots[i];const Knot& left=c.knots[i-1];
  const double gap=right.rigidity-left.rigidity;
  if(gap>mMaxInterpolationGap+kBoundaryTolerance) {
    if(r<=left.rigidity+c.supportHalfWidth+kBoundaryTolerance) {
      mean=left.mean;sigma=left.sigma;return true;
    }
    if(r>=right.rigidity-c.supportHalfWidth-kBoundaryTolerance) {
      mean=right.mean;sigma=right.sigma;return true;
    }
    return false;
  }
  const double fraction=(r-left.rigidity)/gap;
  mean=left.mean+fraction*(right.mean-left.mean);
  sigma=left.sigma+fraction*(right.sigma-left.sigma);
  return true;
}
bool FemtoNuclearDedxCalibration::UsesCalibration(int s,double r) const {
  double mean=0,sigma=0;return Evaluate(s,r,mean,sigma);
}
double FemtoNuclearDedxCalibration::Mean(int s,double r) const {
  if(!GoodSpecies(s)||!std::isfinite(r)||r<=0) return 0;
  double mean=0,sigma=0;if(Evaluate(s,r,mean,sigma)) return mean;
  return NuclearIdDeDxVsMom::GetMean(static_cast<NuclearIdDeDxVsMom::SpeciesIndex>(s),r);
}
double FemtoNuclearDedxCalibration::Sigma(int s,double r) const {
  if(!GoodSpecies(s)||!std::isfinite(r)||r<=0) return 0;
  double mean=0,sigma=0;if(Evaluate(s,r,mean,sigma)) return sigma;
  return NuclearIdDeDxVsMom::GetSigma(static_cast<NuclearIdDeDxVsMom::SpeciesIndex>(s),r);
}
double FemtoNuclearDedxCalibration::Pull(int s,double r,double dedx) const {
  if(!GoodSpecies(s)||!std::isfinite(r)||r<=0||!std::isfinite(dedx)||dedx<=0) return 999.;
  const double mean=Mean(s,r),sigma=Sigma(s,r);
  if(!std::isfinite(mean)||!std::isfinite(sigma)||!(mean>0 && sigma>0)) return 999.;
  const double pull=(dedx-mean)/sigma;
  return std::isfinite(pull)?pull:999.;
}
