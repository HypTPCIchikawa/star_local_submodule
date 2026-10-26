#ifndef STAR_RESTART_FLOAT16_BIN_MAPPING_H
#define STAR_RESTART_FLOAT16_BIN_MAPPING_H
// Diagnostic hypothesis ONLY: default ROOT Float16_t, no range annotation,
// twelve stored mantissa bits. This is not evidence of the input file's writer
// metadata and is not enabled by any production PID configuration.
//
// ROOT's default writer rounds the retained mantissa but CLAMPS an overflow to
// 4095 without carrying into the exponent. Consequently, a boundary at a power
// of two must NOT be mapped to the midpoint of adjacent reconstructed values.
// The exact positive Float32 encoder below, plus an ordered-bit binary search,
// handles that nonuniform last cell directly. ContinuumEdge also includes the
// initial IEEE Float32 round-to-nearest operation. Tie membership at a single
// continuum endpoint has zero probability in the integrated continuous model.
//
// Scope: finite nonnegative Float32 values and finite nonnegative histogram
// edges. Signed values, range-coded Float16_t and nondefault bit counts are not
// modeled. The caller must retain nominal native counts/axes and record the
// nominal AND effective edges when testing this response hypothesis.

#include <cmath>
#include <cstring>
#include <fenv.h>
#include <limits>
#include <stdexcept>
#include <vector>

namespace RestartFloat16 {

inline void Require(bool condition, const char* message) {
  if (!condition) throw std::runtime_error(message);
}

inline void CheckRepresentation() {
  Require(sizeof(float) == 4 && sizeof(unsigned int) == 4 &&
          std::numeric_limits<float>::is_iec559 &&
          std::numeric_limits<float>::digits == 24,
          "Float16 diagnostic requires IEEE binary32 and 32-bit unsigned int");
}

inline unsigned int FloatBits(float value) {
  unsigned int result = 0;
  std::memcpy(&result, &value, sizeof(result));
  return result;
}

inline float FromBits(unsigned int value) {
  float result = 0;
  std::memcpy(&result, &value, sizeof(result));
  return result;
}

// This implements the default ROOT WriteFloat16/ReadFloat16 round trip, not
// standard IEEE half precision. nbits=12 is a representation contract, not a cut.
inline float QuantizeDefault12(float value) {
  CheckRepresentation();
  const unsigned int bits = FloatBits(value);
  Require((bits & 0x80000000U) == 0 && (bits & 0x7f800000U) != 0x7f800000U,
          "Float16 diagnostic accepts only finite nonnegative Float32 input");
  unsigned int mantissa = (bits >> (23 - 12 - 1)) & ((1U << (12 + 1)) - 1U);
  ++mantissa;
  mantissa >>= 1;
  if (mantissa & (1U << 12)) mantissa = (1U << 12) - 1U;
  return FromBits((bits & 0x7f800000U) | (mantissa << (23 - 12)));
}

// The first nonnegative input Float32 for which the reconstructed measurement
// is >= the nominal edge. This is monotone even at the mantissa-overflow clamp.
inline unsigned int FirstInputBits(double nominalEdge) {
  CheckRepresentation();
  Require(std::isfinite(nominalEdge) && nominalEdge >= 0.,
          "Float16 nominal edge must be finite and nonnegative");
  const unsigned int lastFinite = 0x7f7fffffU;
  Require(nominalEdge <= double(QuantizeDefault12(FromBits(lastFinite))),
          "Float16 nominal edge exceeds the largest reconstructed finite value");
  unsigned int lo = 0, hi = lastFinite;
  while (lo < hi) {
    const unsigned int mid = lo + (hi - lo) / 2U;
    if (double(QuantizeDefault12(FromBits(mid))) >= nominalEdge) hi = mid;
    else lo = mid + 1U;
  }
  return lo;
}

inline double ContinuumEdge(double nominalEdge) {
  Require(fegetround() == FE_TONEAREST,
          "Float16 continuum mapping requires round-to-nearest Float32 conversion");
  const unsigned int first = FirstInputBits(nominalEdge);
  if (!first) return 0.;
  // The two Float32 values are exactly representable in double. Their midpoint
  // gives the continuous-input boundary before the initial Float32 conversion.
  return .5 * (double(FromBits(first - 1U)) + double(FromBits(first)));
}

inline std::vector<double> IntegrationEdges(const std::vector<double>& nominal) {
  Require(nominal.size() >= 2, "Float16 mapping needs at least two native edges");
  std::vector<double> result;
  result.reserve(nominal.size());
  for (unsigned i = 0; i < nominal.size(); ++i) {
    if (i) Require(nominal[i] > nominal[i - 1],
                   "Nominal Float16 histogram edges must increase strictly");
    const double effective = ContinuumEdge(nominal[i]);
    if (i) Require(effective >= result.back(),
                   "Effective Float16 preimage edges are not monotone");
    result.push_back(effective);
  }
  return result;
}

}  // namespace RestartFloat16
#endif
