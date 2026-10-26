// In-memory ROOT5 oracle tests. No output ROOT/data/PID files are written.
int TestFloat16BinMapping();
#ifndef __CINT__
#include "Float16BinMapping.h"
#include <TBufferFile.h>
#include <TAxis.h>
#include <TMath.h>
#include <TSystem.h>
#include <algorithm>
#include <iostream>
#include <string>

namespace {
unsigned Float16Checks = 0;

void Float16Check(bool condition, const std::string& label) {
  ++Float16Checks;
  if (!condition) throw std::runtime_error(label);
}

void Float16Oracle(const std::vector<float>& values) {
  TBufferFile output(TBuffer::kWrite);
  for (unsigned i = 0; i < values.size(); ++i) {
    Float_t value = values[i];
    output.WriteFloat16(&value, 0);
  }
  Float16Check(output.Length() == int(3 * values.size()),
               "ROOT default Float16 writer uses three bytes per input");
  TBufferFile input(TBuffer::kRead, output.Length(), output.Buffer(), kFALSE);
  for (unsigned i = 0; i < values.size(); ++i) {
    Float_t decoded = -1;
    input.ReadFloat16(&decoded, 0);
    const float predicted = RestartFloat16::QuantizeDefault12(values[i]);
    Float16Check(RestartFloat16::FloatBits(decoded) == RestartFloat16::FloatBits(predicted),
                 "Bit-exact encoder disagrees with ROOT WriteFloat16/ReadFloat16");
  }
}

void Float16EncoderTests() {
  std::vector<float> values;
  values.push_back(0.f);
  values.push_back(RestartFloat16::FromBits(1U));
  values.push_back(RestartFloat16::FromBits(0x007fffffU));
  values.push_back(RestartFloat16::FromBits(0x00800000U));
  values.push_back(RestartFloat16::FromBits(0x7f7fffffU));
  for (int exponent = -20; exponent <= 20; ++exponent) {
    const float base = float(std::ldexp(1., exponent));
    const unsigned int boundary = RestartFloat16::FloatBits(base);
    for (int offset = -8; offset <= 8; ++offset)
      values.push_back(RestartFloat16::FromBits(boundary + offset));
    const double step = std::ldexp(1., exponent - 12);
    const int levels[] = {0, 1, 2, 1024, 2048, 4094, 4095};
    for (unsigned j = 0; j < sizeof(levels) / sizeof(levels[0]); ++j) {
      const double transition = double(base) + (levels[j] + .5) * step;
      const unsigned int bits = RestartFloat16::FloatBits(float(transition));
      for (int offset = -2; offset <= 2; ++offset)
        values.push_back(RestartFloat16::FromBits(bits + offset));
    }
    const float lastInput = RestartFloat16::FromBits(
      RestartFloat16::FloatBits(float(2. * base)) - 1U);
    const float clamped = float(2. * base - step);
    Float16Check(RestartFloat16::QuantizeDefault12(lastInput) == clamped,
                 "Mantissa overflow must clamp within the old exponent");
    Float16Check(RestartFloat16::QuantizeDefault12(float(2. * base)) == float(2. * base),
                 "The next exponent must begin exactly at the power of two");
  }
  unsigned int state = 0x93b14a27U;
  for (unsigned i = 0; i < 20000; ++i) {
    state = 1664525U * state + 1013904223U;
    values.push_back(RestartFloat16::FromBits(state % 0x7f800000U));
  }
  Float16Oracle(values);
}

void Float16MappingTests() {
  TAxis axis(1200, 0., 120.);
  std::vector<double> nominal;
  for (int bin = 1; bin <= axis.GetNbins(); ++bin) nominal.push_back(axis.GetBinLowEdge(bin));
  nominal.push_back(axis.GetBinUpEdge(axis.GetNbins()));
  const std::vector<double> effective = RestartFloat16::IntegrationEdges(nominal);
  Float16Check(effective.size() == nominal.size(), "Mapping preserves native edge cardinality");
  for (unsigned i = 0; i < nominal.size(); ++i) {
    const unsigned int first = RestartFloat16::FirstInputBits(nominal[i]);
    Float16Check(double(RestartFloat16::QuantizeDefault12(RestartFloat16::FromBits(first))) >= nominal[i],
                 "First qualifying Float32 must enter or pass the measurement edge");
    if (first) {
      Float16Check(double(RestartFloat16::QuantizeDefault12(RestartFloat16::FromBits(first - 1U))) < nominal[i],
                   "Predecessor must remain below the measurement edge");
      Float16Check(double(RestartFloat16::FromBits(first - 1U)) < effective[i] &&
                    effective[i] < double(RestartFloat16::FromBits(first)),
                   "Effective edge must bisect neighboring input Float32 values");
    }
    if (i) Float16Check(effective[i] >= effective[i - 1], "Mapped native bins are contiguous and monotone");
  }
  double widthSum = 0.;
  for (unsigned i = 1; i < effective.size(); ++i) widthSum += effective[i] - effective[i - 1];
  Float16Check(std::fabs(widthSum - (effective.back() - effective.front())) < 1.e-12,
               "All effective widths telescope without overlap or holes");
  for (int exponent = -10; exponent <= 6; ++exponent) {
    const double power = std::ldexp(1., exponent);
    const unsigned int bits = RestartFloat16::FloatBits(float(power));
    Float16Check(RestartFloat16::FirstInputBits(power) == bits,
                 "Exponent boundary must not collect a rounded previous-exponent cell");
    const double mapped = RestartFloat16::ContinuumEdge(power);
    const double midpoint = .5 * (double(RestartFloat16::FromBits(bits - 1U)) + power);
    Float16Check(mapped == midpoint, "Power-of-two boundary includes only the initial Float32 rounding");
  }
  unsigned int state = 0x471bc23aU;
  for (unsigned i = 0; i < 30000; ++i) {
    state = 1664525U * state + 1013904223U;
    const double x = 120. * (double(state) + .37) / 4294967296.;
    const double stored = double(RestartFloat16::QuantizeDefault12(float(x)));
    const int observed = axis.FindFixBin(stored);
    const int mapped = int(std::upper_bound(effective.begin(), effective.end(), x) - effective.begin());
    Float16Check(observed == mapped, "Native ROOT bin membership disagrees with continuum preimage");
  }
  // A half-keV block away from exponent changes covers an exact number of
  // quantization cells. This is a mathematical representation check, not rebinning data.
  for (int i = 1; i < 16; ++i) {
    const double lo = 32. + .5 * i, hi = lo + .5;
    Float16Check(std::fabs(RestartFloat16::ContinuumEdge(hi) -
                           RestartFloat16::ContinuumEdge(lo) - .5) < 1.e-12,
                 "Half-keV effective widths away from exponent changes are equal");
  }
  bool rejected = false;
  try { RestartFloat16::QuantizeDefault12(-1.f); } catch (const std::exception&) { rejected = true; }
  Float16Check(rejected, "Signed inputs are outside this hypothesis");
  rejected = false;
  try { RestartFloat16::ContinuumEdge(-.1); } catch (const std::exception&) { rejected = true; }
  Float16Check(rejected, "Negative native edges are rejected");
  rejected = false;
  try {
    std::vector<double> reversed(2, 1.);
    RestartFloat16::IntegrationEdges(reversed);
  } catch (const std::exception&) { rejected = true; }
  Float16Check(rejected, "Repeated nominal edges are rejected");
}
}  // namespace
#endif

int TestFloat16BinMapping() {
#ifndef __CINT__
  try {
    Float16Checks = 0;
    RestartFloat16::CheckRepresentation();
    Float16EncoderTests();
    Float16MappingTests();
    std::cout << "FLOAT16_BIN_MAPPING_TEST_PASS checks=" << Float16Checks
              << " defaultMantissaBits=12 hypothesisOnly=1 filesWritten=0" << std::endl;
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "FLOAT16_BIN_MAPPING_TEST_FAIL " << error.what() << std::endl;
    return 1;
  }
#else
  return 1;
#endif
}
