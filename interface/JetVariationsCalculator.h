#ifndef JETVARIATIONSCALCULATOR_H
#define JETVARIATIONSCALCULATOR_H
#include "JetMETVariationsCalculatorBase.h"
#include "Utility.h"

class JetVariationsCalculator : public JetMETVariationsCalculatorBase {
public:
  using result_t = rdfhelpers::ModifiedPtMCollection;

  JetVariationsCalculator() = default;

  std::vector<std::string> available(const std::string& attr = {}) const;
  // interface for NanoAOD
  result_t produce(
      const p4compv_t& jet_pt, const p4compv_t& jet_eta, const p4compv_t& jet_phi, const p4compv_t& jet_mass,
      const p4compv_t& jet_rawcorr, const p4compv_t& jet_area, const ROOT::VecOps::RVec<int>& jet_jetId, const float rho,
      // MC-only
      const std::uint32_t seed,
      const p4compv_t& genjet_pt, const p4compv_t& genjet_eta, const p4compv_t& genjet_phi, const p4compv_t& genjet_mass
      ) const;
};

#endif
