#ifndef CMSJMECalculators_JMESystematicsCalculators_H
#define CMSJMECalculators_JMESystematicsCalculators_H

#include <map>
#include <ROOT/RVec.hxx>
#include "JetResolution.h"
#include "JetCorrectorParameters.h"
#include "SimpleJetCorrectionUncertainty.h"
#include "FormulaEvaluator.h"

class FactorizedJetCorrectorCalculator;
class TRandom3;

namespace rdfhelpers {

class ModifiedPtCollection {
public:
  using compv_t = ROOT::VecOps::RVec<float>;

  ModifiedPtCollection() = default;
  ModifiedPtCollection(std::size_t n, const compv_t& pt) : m_pt(n, pt) {}

  std::size_t size() const { return m_pt.size(); }
  const compv_t& pt(std::size_t i) const { return m_pt[i]; }

  void set(std::size_t i, const compv_t& pt) { m_pt[i] = pt; }
  void set(std::size_t i, compv_t&& pt) { m_pt[i] = std::move(pt); }
private:
  std::vector<compv_t> m_pt;
};

class ModifiedPtMCollection { // map of variation collections
public:
  using compv_t = ROOT::VecOps::RVec<float>;

  ModifiedPtMCollection() = default;
  ModifiedPtMCollection(std::size_t n, const compv_t& pt, const compv_t& mass)
    : m_pt(n, pt), m_mass(n, mass) {}

  std::size_t size() const { return m_pt.size(); }

  const compv_t& pt(std::size_t i) const { return m_pt[i]; }
  const compv_t& mass(std::size_t i) const { return m_mass[i]; }

  void set(std::size_t i, const compv_t& pt, const compv_t& mass) {
    m_pt[i] = pt;
    m_mass[i] = mass;
  }
  void set(std::size_t i, compv_t&& pt, compv_t&& mass) {
    m_pt[i] = std::move(pt);
    m_mass[i] = std::move(mass);
  }
private:
  std::vector<compv_t> m_pt;
  std::vector<compv_t> m_mass;
};

class ModifiedPtMMsdCollection { // for fat jets
public:
  using compv_t = ROOT::VecOps::RVec<float>;

  ModifiedPtMMsdCollection() = default;
  ModifiedPtMMsdCollection(std::size_t nPt, const compv_t& pt, const std::size_t nM, const compv_t& mass, const compv_t& msd)
    : m_pt(nPt, pt), m_mass(nM, mass), m_msd(nM, msd) {}

  std::size_t size() const { return m_pt.size(); }
  std::size_t sizeM() const { return m_mass.size(); }

  const compv_t& pt(std::size_t i) const { return m_pt[i]; }
  const compv_t& mass(std::size_t i) const { return m_mass[i]; }
  const compv_t& msoftdrop(std::size_t i) const { return m_msd[i]; }

  void set(std::size_t i, const compv_t& pt, const compv_t& mass, const compv_t& msd) {
    m_pt[i] = pt;
    m_mass[i] = mass;
    m_msd[i] = msd;
  }
  void setM(std::size_t i, const compv_t& mass, const compv_t& msd) {
    m_mass[i] = mass;
    m_msd[i] = msd;
  }
  void set(std::size_t i, compv_t&& pt, compv_t&& mass, compv_t&& msd) {
    m_pt[i] = std::move(pt);
    m_mass[i] = std::move(mass);
    m_msd[i] = std::move(msd);
  }
  void setM(std::size_t i, compv_t&& mass, compv_t&& msd) {
    m_mass[i] = std::move(mass);
    m_msd[i] = std::move(msd);
  }
private:
  std::vector<compv_t> m_pt;
  std::vector<compv_t> m_mass;
  std::vector<compv_t> m_msd;
};

class ModifiedMET {
public:
  using compv_t = ROOT::VecOps::RVec<double>;

  ModifiedMET() = default;
  // initialize with the nominal value for all variations
  ModifiedMET(std::size_t n, double px_nom, double py_nom)
    : m_px(n, px_nom), m_py(n, py_nom) {}

  std::size_t size() const { return m_px.size(); }
  const compv_t& px() const { return m_px; }
  const compv_t& py() const { return m_py; }
  double px (std::size_t i) const { return m_px[i]; }
  double py (std::size_t i) const { return m_py[i]; }
  double pt (std::size_t i) const { return std::sqrt(m_px[i]*m_px[i]+m_py[i]*m_py[i]); }
  double phi(std::size_t i) const { return std::atan2(m_py[i], m_px[i]); }

  void setXY(std::size_t i, double dpx, double dpy) {
    m_px[i] = dpx;
    m_py[i] = dpy;
  }
  void addR_proj(std::size_t i, double cosphi, double sinphi, double dp) {
    m_px[i] += dp*cosphi;
    m_py[i] += dp*sinphi;
  }
private:
  compv_t m_px;
  compv_t m_py;
};
}

class JetMETVariationsCalculatorBase {
public:
  using p4compv_t = ROOT::VecOps::RVec<float>;

  JetMETVariationsCalculatorBase() = default;
  ~JetMETVariationsCalculatorBase();

  // set up smearing (and JER systematics)
  void setSmearing(const std::string& ptResolution, const std::string& ptResolutionSF, bool splitJER, bool doGenMatch, float genMatch_maxDR=-1., float genMatch_maxDPT=-1.)
  {
    m_doSmearing = true;
    m_jetPtRes   = JME::JetResolution(ptResolution);
    m_jetEResSF = JME::JetResolutionScaleFactor(ptResolutionSF);
    m_splitJER = splitJER;
    m_smearDoGenMatch = doGenMatch;
    m_genMatch_dR2max = genMatch_maxDR*genMatch_maxDR;
    m_genMatch_dPtmax = genMatch_maxDPT;
  }

  void setJEC(const std::vector<JetCorrectorParameters>& jecParams);
  void setAddHEM2018Issue(bool enable) { m_addHEM2018Issue = enable; }

  void addJESUncertainty(const std::string& name, const JetCorrectorParameters& params)
  {
    m_jesUncSources.emplace(std::piecewise_construct,
        std::forward_as_tuple(name),
        std::forward_as_tuple(params));
  }
protected:
  std::size_t findGenMatch(const double pt, const float eta, const float phi, const ROOT::VecOps::RVec<float>& gen_pt, const ROOT::VecOps::RVec<float>& gen_eta, const ROOT::VecOps::RVec<float>& gen_phi, const double resolution ) const;

  // config options
  bool m_doSmearing{false}, m_smearDoGenMatch;      // default: yes, yes
  bool m_addHEM2018Issue{false}, m_splitJER{false}; // default: no, no
  float m_genMatch_dR2max, m_genMatch_dPtmax;       // default: R/2 (0.2) and 3
  // parameters and helpers
  JME::JetResolution m_jetPtRes;
  JME::JetResolutionScaleFactor m_jetEResSF;
  struct jetcorrdeleter { void operator()(FactorizedJetCorrectorCalculator*) const; };
  std::unique_ptr<FactorizedJetCorrectorCalculator,jetcorrdeleter> m_jetCorrector;
  std::unordered_map<std::string,SimpleJetCorrectionUncertainty> m_jesUncSources;
};

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
protected:
  float m_unclEnThreshold = 15.;
  bool m_isT1SmearedMET = false;
  std::unique_ptr<FactorizedJetCorrectorCalculator,jetcorrdeleter> m_jetCorrectorL1;
  void addVariations(Type1METVariationsCalculator::result_t& out,
      const p4compv_t& jet_pt, const p4compv_t& jet_eta, const p4compv_t& jet_phi, const p4compv_t& jet_mass,
      const p4compv_t& jet_rawcorr, const p4compv_t& jet_area, const p4compv_t& jet_muonSubtrFactor,
      const p4compv_t& jet_neEmEF, const p4compv_t& jet_chEmEF,
      const ROOT::VecOps::RVec<int>& jet_jetId, const ROOT::VecOps::RVec<bool>& jet_mask,
      const float rho,
      const p4compv_t& genjet_pt, const p4compv_t& genjet_eta, const p4compv_t& genjet_phi, const p4compv_t& genjet_mass,
      TRandom3& rg) const;
};

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
protected:
  std::unique_ptr<FactorizedJetCorrectorCalculator,jetcorrdeleter> m_jetCorrectorProd;
  std::unique_ptr<FactorizedJetCorrectorCalculator,jetcorrdeleter> m_jetCorrectorL1Prod;
  std::array<double,4> calculateFixEE2017Offset(ROOT::VecOps::RVec<bool>& jet_mask,
      const p4compv_t& jet_pt, const p4compv_t& jet_eta, const p4compv_t& jet_phi, const p4compv_t& jet_mass,
      const p4compv_t& jet_rawcorr, const p4compv_t& jet_area, const p4compv_t& jet_muonSubtrFactor,
      const float rho) const;
};

class FatJetVariationsCalculator : public JetMETVariationsCalculatorBase {
public:
  using result_t = rdfhelpers::ModifiedPtMMsdCollection;

  FatJetVariationsCalculator() = default;

  void setJMRValues(double nominal, double up=1., double down=1.) { m_jmrVals = { nominal, up, down }; }
  void setGMRValues(double nominal, double up=1., double down=1.) { m_gmrVals = { nominal, up, down }; }
  void setJMSValues(double nominal, double up=1., double down=1.) { m_jmsVals = { nominal, up/nominal, down/nominal }; }
  void setGMSValues(double nominal, double up=1., double down=1.) { m_gmsVals = { nominal, up/nominal, down/nominal }; }
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
  std::array<double,3> m_jmrVals = {1., 1., 1.}; // nominal, up, down
  std::array<double,3> m_gmrVals = {1., 1., 1.}; // nominal, up, down
  std::array<double,3> m_jmsVals = {1., 1., 1.}; // nominal, up/nominal, down/nominal
  std::array<double,3> m_gmsVals = {1., 1., 1.}; // nominal, up/nominal, down/nominal
};

#endif // CMSJMECalculators_JMESystematicsCalculators_H
