#ifndef FATJETVARIATIONSCALCULATOR_H
#define FATJETVARIATIONSCALCULATOR_H
#include "JetMETVariationsCalculatorBase.h"
#include "Utility.h"
class FatJetVariationsCalculator : public JetMETVariationsCalculatorBase {
public:
  using result_t = rdfhelpers::ModifiedPtMMsdCollection;

  FatJetVariationsCalculator() = default;

  void setJMRValues(double nominal, double up=1., double down=1.) { m_jmrVals = {{ nominal, up, down }}; }
  void setGMRValues(double nominal, double up=1., double down=1.) { m_gmrVals = {{ nominal, up, down }}; }
  void setJMSValues(double nominal, double up=1., double down=1.) { m_jmsVals = {{ nominal, up/nominal, down/nominal }}; }
  void setGMSValues(double nominal, double up=1., double down=1.) { m_gmsVals = {{ nominal, up/nominal, down/nominal }}; }
  void setPuppiCorrections(const std::string& genFormula, const std::array<double, 6>& reco_cen_params, const std::array<double, 6>& reco_for_params, const std::array<double, 6>& resol_cen_params, const std::array<double, 6>& resol_for_params);

  std::vector<std::string> available(const std::string& attr = {}) const;
  // interface for NanoAOD
  result_t produce(
      const p4compv_t& jet_pt, const p4compv_t& jet_eta, const p4compv_t& jet_phi, const p4compv_t& jet_mass,
      const p4compv_t& jet_rawcorr, const p4compv_t& jet_area,
      const p4compv_t& jet_msoftdrop, const ROOT::VecOps::RVec<int>& jet_subJetIdx1, const ROOT::VecOps::RVec<int>& jet_subJetIdx2,
      const p4compv_t& subjet_pt, const p4compv_t& subjet_eta, const p4compv_t& subjet_phi, const p4compv_t& subjet_mass,
      const ROOT::VecOps::RVec<int>& jet_jetId, const float rho,
      // MC-only
      const std::uint32_t seed,
      const p4compv_t& genjet_pt, const p4compv_t& genjet_eta, const p4compv_t& genjet_phi, const p4compv_t& genjet_mass,
      const p4compv_t& gensubjet_pt, const p4compv_t& gensubjet_eta, const p4compv_t& gensubjet_phi, const p4compv_t& gensubjet_mass
      ) const;
private:
  std::unique_ptr<reco::FormulaEvaluator> m_puppiCorrGen, m_puppiPoly5;
  std::array<double,6> m_puppiCorrRecoParam_cen, m_puppiCorrRecoParam_for, m_puppiResolParam_cen, m_puppiResolParam_for;
  std::array<double,3> m_jmrVals = {{1., 1., 1.}}; // nominal, up, down
  std::array<double,3> m_gmrVals = {{1., 1., 1.}}; // nominal, up, down
  std::array<double,3> m_jmsVals = {{1., 1., 1.}}; // nominal, up/nominal, down/nominal
  std::array<double,3> m_gmsVals = {{1., 1., 1.}}; // nominal, up/nominal, down/nominal
};
#endif
