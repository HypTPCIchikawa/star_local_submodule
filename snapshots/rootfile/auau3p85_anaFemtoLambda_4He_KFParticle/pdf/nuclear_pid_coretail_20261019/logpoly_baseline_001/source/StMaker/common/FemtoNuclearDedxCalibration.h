#ifndef FEMTO_NUCLEAR_DEDX_CALIBRATION_H
#define FEMTO_NUCLEAR_DEDX_CALIBRATION_H

#include <iosfwd>
#include <string>
#include <vector>

// Opt-in calibration owned by the Lambda-nucleus compatibility provider.
// The shared NuclearIdDeDxVsMom tables remain the legacy source of truth.
class FemtoNuclearDedxCalibration {
 public:
  FemtoNuclearDedxCalibration();
  void ResetLegacy();
  bool Load(const std::string& path, std::ostream& diagnostics);
  double Mean(int species, double rigidity) const;
  double Sigma(int species, double rigidity) const;
  double Pull(int species, double rigidity, double dedx) const;
  bool UsesCalibration(int species, double rigidity) const;
  const std::string& Metadata() const { return mMetadata; }

 private:
  struct Knot { double rigidity, mean, sigma; };
  struct Curve {
    double minimum, maximum, supportHalfWidth;
    std::vector<Knot> knots;
    Curve() : minimum(0), maximum(0), supportHalfWidth(0) {}
  };
  bool Evaluate(int species, double rigidity, double& mean, double& sigma) const;
  bool mEnabled;
  double mMaxInterpolationGap;
  Curve mCurve[2];  // He3, He4; d and t deliberately retain legacy calibration.
  std::string mMetadata;
};

#endif
