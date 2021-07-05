#ifndef TYPE1METVARIATIONSCALCULATOR_H
#define TYPE1METVARIATIONSCALCULATOR_H
#include "JetMETVariationsCalculatorBase.h"
#include "Utility.h"
class Type1METVariationsCalculator : public JetMETVariationsCalculatorBase {
public:
  using result_t = rdfhelpers::ModifiedMET;

  Type1METVariationsCalculator() = default;

  // additional settings: L1-only JEC
  void setL1JEC(const std::vector<JetCorrectorParameters>& jecParams);
  void setUnclusteredEnergyTreshold(float threshold) { m_unclEnThreshold = threshold; }
  void setIsT1SmearedMET(bool isT1SmearedMET) { m_isT1SmearedMET = isT1SmearedMET; }

  std::vector<std::string> available(const std::string& attr = {}) const;
  // interface for NanoAOD
  result_t produce(
      const p4compv_t& jet_pt, const p4compv_t& jet_eta, const p4compv_t& jet_phi, const p4compv_t& jet_mass,
      const p4compv_t& jet_rawcorr, const p4compv_t& jet_area,
      const p4compv_t& jet_muonSubtrFactor, const p4compv_t& jet_neEmEF, const p4compv_t& jet_chEmEF, const ROOT::VecOps::RVec<int>& jet_jetId,
      const float rho,
      // MC-only
      const std::uint32_t seed,
      const p4compv_t& genjet_pt, const p4compv_t& genjet_eta, const p4compv_t& genjet_phi, const p4compv_t& genjet_mass,
      // MET-specific
      const float rawmet_phi, const float rawmet_pt,
      const float met_unclustenupdx, const float met_unclustenupdy,
      const p4compv_t& lowptjet_rawpt, const p4compv_t& lowptjet_eta, const p4compv_t& lowptjet_phi, const p4compv_t& lowptjet_area,
      const p4compv_t& lowptjet_muonSubtrFactor, const p4compv_t& lowptjet_neEmEF, const p4compv_t& lowptjet_chEmEF
      ) const;
  //protected:
  float m_unclEnThreshold = 15.;
  bool m_isT1SmearedMET = false;
  //std::unique_ptr<FactorizedJetCorrectorCalculator,jetcorrdeleter> m_jetCorrectorL1;
  std::unique_ptr<FactorizedJetCorrectorCalculator> m_jetCorrectorL1;
  void addVariations(Type1METVariationsCalculator::result_t& out,
      const p4compv_t& jet_pt, const p4compv_t& jet_eta, const p4compv_t& jet_phi, const p4compv_t& jet_mass,
      const p4compv_t& jet_rawcorr, const p4compv_t& jet_area, const p4compv_t& jet_muonSubtrFactor,
      const p4compv_t& jet_neEmEF, const p4compv_t& jet_chEmEF,
      const ROOT::VecOps::RVec<int>& jet_jetId, const ROOT::VecOps::RVec<bool>& jet_mask,
      const float rho,
      const p4compv_t& genjet_pt, const p4compv_t& genjet_eta, const p4compv_t& genjet_phi, const p4compv_t& genjet_mass,
      TRandom3& rg) const;
};
#endif
