#ifndef JetMETVariationsCalculatorBase_H
#define JetMETVariationsCalculatorBase_H

#include <map>
#include <ROOT/RVec.hxx>
#if defined(PROJECT_NAME) && defined(CMSSW_GIT_HASH)
#include "JetMETCorrections/Modules/interface/JetResolution.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/SimpleJetCorrectionUncertainty.h"
#include "CommonTools/Utils/interface/FormulaEvaluator.h"
#else
#include "JetResolution.h"
#include "JetCorrectorParameters.h"
#include "SimpleJetCorrectionUncertainty.h"
#include "FormulaEvaluator.h"
#include "FactorizedJetCorrectorCalculator.h"
#endif

#include "Utility.h"
#ifdef BAMBOO_JME_DEBUG
#define LogDebug std::cout
#else
#define LogDebug if (false) std::cout
#endif



class TRandom3;

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
  //protected:
  std::size_t findGenMatch(const double pt, const float eta, const float phi, const ROOT::VecOps::RVec<float>& gen_pt, const ROOT::VecOps::RVec<float>& gen_eta, const ROOT::VecOps::RVec<float>& gen_phi, const double resolution ) const;

  // config options
  bool m_doSmearing{false}, m_smearDoGenMatch;      // default: yes, yes
  bool m_addHEM2018Issue{false}, m_splitJER{false}; // default: no, no
  float m_genMatch_dR2max, m_genMatch_dPtmax;       // default: R/2 (0.2) and 3
  // parameters and helpers
  JME::JetResolution m_jetPtRes;
  JME::JetResolutionScaleFactor m_jetEResSF;
  struct jetcorrdeleter { void operator()(FactorizedJetCorrectorCalculator*) const; };
  //std::unique_ptr<FactorizedJetCorrectorCalculator,jetcorrdeleter> m_jetCorrector;
  std::unique_ptr<FactorizedJetCorrectorCalculator> m_jetCorrector;
  std::unordered_map<std::string,SimpleJetCorrectionUncertainty> m_jesUncSources;
};







#endif // CMSJMECalculators_JMESystematicsCalculators_H
