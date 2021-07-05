#ifndef FIXEE2017TYPE1METVARIATIONSCALCULATOR_H
#define FIXEE2017TYPE1METVARIATIONSCALCULATOR_H
#include "JetMETVariationsCalculatorBase.h"
#include "Type1METVariationsCalculator.h"
#include "Utility.h"
class FixEE2017Type1METVariationsCalculator : public Type1METVariationsCalculator {
public:
  FixEE2017Type1METVariationsCalculator() = default;

  // additional settings: full and L1-only JEC used in production
  void setJECProd(const std::vector<JetCorrectorParameters>& jecParams);
  void setL1JECProd(const std::vector<JetCorrectorParameters>& jecParams);

  // interface for NanoAOD
  result_t produce(
      const p4compv_t& jet_pt, const p4compv_t& jet_eta, const p4compv_t& jet_phi, const p4compv_t& jet_mass,
      const p4compv_t& jet_rawcorr, const p4compv_t& jet_area,
      const p4compv_t& jet_muonSubtrFactor, const p4compv_t& jet_neEmEF, const p4compv_t& jet_chEmEF,
      const float rho,
      // MC-only
      const std::uint32_t seed,
      const p4compv_t& genjet_pt, const p4compv_t& genjet_eta, const p4compv_t& genjet_phi, const p4compv_t& genjet_mass,
      // MET-specific
      const float rawmet_phi, const float rawmet_pt, // "RawMET"
      const float met_unclustenupdx, const float met_unclustenupdy,
      const p4compv_t& lowptjet_rawpt, const p4compv_t& lowptjet_eta, const p4compv_t& lowptjet_phi, const p4compv_t& lowptjet_area,
      const p4compv_t& lowptjet_muonSubtrFactor, const p4compv_t& lowptjet_neEmEF, const p4compv_t& lowptjet_chEmEF,
      const float defmet_phi, const float defmet_pt, // "MET"
      const float t1met_phi, const float t1met_pt    // "METFixEE2017"
      ) const;
  //protected:
  //std::unique_ptr<FactorizedJetCorrectorCalculator,jetcorrdeleter> m_jetCorrectorProd;
  //std::unique_ptr<FactorizedJetCorrectorCalculator,jetcorrdeleter> m_jetCorrectorL1Prod;
  std::unique_ptr<FactorizedJetCorrectorCalculator> m_jetCorrectorProd;
  std::unique_ptr<FactorizedJetCorrectorCalculator> m_jetCorrectorL1Prod;
  std::array<double,4> calculateFixEE2017Offset(ROOT::VecOps::RVec<bool>& jet_mask,
      const p4compv_t& jet_pt, const p4compv_t& jet_eta, const p4compv_t& jet_phi, const p4compv_t& jet_mass,
      const p4compv_t& jet_rawcorr, const p4compv_t& jet_area, const p4compv_t& jet_muonSubtrFactor,
      const float rho) const;
};

#endif
