#if defined(PROJECT_NAME) && defined(CMSSW_GIT_HASH)
#include "UserCode/CMSJMECalculators/interface/JetMETVariationsCalculatorBase.h"
#else
#include "JetMETVariationsCalculatorBase.h"
#endif
#include <Math/GenVector/LorentzVector.h>
#include <Math/GenVector/PtEtaPhiM4D.h>
#include "Math/VectorUtil.h"
#include "TRandom3.h"
#include <cassert>

JetMETVariationsCalculatorBase::~JetMETVariationsCalculatorBase()
{}

void JetMETVariationsCalculatorBase::setJEC(const std::vector<JetCorrectorParameters>& jecParams)
{
  if ( ! jecParams.empty() ) {
    //m_jetCorrector = std::unique_ptr<FactorizedJetCorrectorCalculator,jetcorrdeleter>{new FactorizedJetCorrectorCalculator(jecParams)};
    m_jetCorrector = std::unique_ptr<FactorizedJetCorrectorCalculator>{new FactorizedJetCorrectorCalculator(jecParams)};
  }
}

void JetMETVariationsCalculatorBase::jetcorrdeleter::operator() (FactorizedJetCorrectorCalculator* ptr) const
{ delete ptr; }

// TODO with orig MET and jets (sumpx,sumpy): calc modif MET(sig), produce bigger results type

std::size_t JetMETVariationsCalculatorBase::findGenMatch(const double pt, const float eta, const float phi, const ROOT::VecOps::RVec<float>& gen_pt, const ROOT::VecOps::RVec<float>& gen_eta, const ROOT::VecOps::RVec<float>& gen_phi, const double resolution ) const
{
  auto dr2Min = std::numeric_limits<float>::max();
  std::size_t igBest{gen_pt.size()};
  LogDebug << "(DRs: ";
  for ( std::size_t ig{0}; ig != gen_pt.size(); ++ig ) {
    const auto dphi = phi_mpi_pi(gen_phi[ig]-phi);
    const auto deta = (gen_eta[ig]-eta);
    const auto dr2 = dphi*dphi + deta*deta;
    LogDebug << "dr2=" << dr2;
    if ( ( dr2 < dr2Min ) && ( dr2 < m_genMatch_dR2max ) ) {
      LogDebug << "->dpt=" << std::abs(gen_pt[ig]-pt) << ",res=" << resolution;
      if ( std::abs(gen_pt[ig]-pt) < m_genMatch_dPtmax*resolution ) {
        LogDebug << "->best:" << ig;
        dr2Min = dr2;
        igBest = ig;
      }
    }
    LogDebug << ", ";
  }
  LogDebug << ")";
  return igBest;
}
