#if defined(PROJECT_NAME) && defined(CMSSW_GIT_HASH)
#include "UserCode/CMSJMECalculators/interface/FixEE2017Type1METVariationsCalculator.h"
#else
#include "FixEE2017Type1METVariationsCalculator.h"
#endif
#include <Math/GenVector/LorentzVector.h>
#include <Math/GenVector/PtEtaPhiM4D.h>
#include "Math/VectorUtil.h"
#include "TRandom3.h"
#include <cassert>

void FixEE2017Type1METVariationsCalculator::setJECProd(const std::vector<JetCorrectorParameters>& jecParams)
{
  if ( ! jecParams.empty() ) {
    //m_jetCorrectorProd = std::unique_ptr<FactorizedJetCorrectorCalculator,jetcorrdeleter>{new FactorizedJetCorrectorCalculator(jecParams)};
    m_jetCorrectorProd = std::unique_ptr<FactorizedJetCorrectorCalculator>{new FactorizedJetCorrectorCalculator(jecParams)};
  }
}
void FixEE2017Type1METVariationsCalculator::setL1JECProd(const std::vector<JetCorrectorParameters>& jecParams)
{
  if ( ! jecParams.empty() ) {
    //m_jetCorrectorL1Prod = std::unique_ptr<FactorizedJetCorrectorCalculator,jetcorrdeleter>{new FactorizedJetCorrectorCalculator(jecParams)};
    m_jetCorrectorL1Prod = std::unique_ptr<FactorizedJetCorrectorCalculator>{new FactorizedJetCorrectorCalculator(jecParams)};
  }
}

FixEE2017Type1METVariationsCalculator::result_t FixEE2017Type1METVariationsCalculator::produce(
    const p4compv_t& jet_pt, const p4compv_t& jet_eta, const p4compv_t& jet_phi, const p4compv_t& jet_mass,
    const p4compv_t& jet_rawcorr, const p4compv_t& jet_area,
    const p4compv_t& jet_muonSubtrFactor, const p4compv_t& jet_neEmEF, const p4compv_t& jet_chEmEF,
    const float rho,
    const std::uint32_t seed,
    const p4compv_t& genjet_pt, const p4compv_t& genjet_eta, const p4compv_t& genjet_phi, const p4compv_t& genjet_mass,
    const float rawmet_phi, const float rawmet_pt,
    const float met_unclustenupdx, const float met_unclustenupdy,
    const p4compv_t& lowptjet_rawpt, const p4compv_t& lowptjet_eta, const p4compv_t& lowptjet_phi, const p4compv_t& lowptjet_area,
    const p4compv_t& lowptjet_muonSubtrFactor, const p4compv_t& lowptjet_neEmEF, const p4compv_t& lowptjet_chEmEF,
    const float defmet_phi, const float defmet_pt, const float t1met_phi, const float t1met_pt
    ) const
{
  LogDebug << "JME:: hello from Type1METVariations produce with 2017 EE Fix. Got " << jet_pt.size() << " jets and " << lowptjet_rawpt.size() << " low-PT jets" << std::endl;
  assert(!m_addHEM2018Issue);
  const auto nVariations = 3+( m_doSmearing ? 2*( m_splitJER ? 6 : 1 ) : 0 )+2*m_jesUncSources.size(); // 1(nom)+2(unclust)+3(JER)+2*len(JES)
  p4compv_t lowptjet_zero(lowptjet_rawpt.size(), 0.);
  // the actual MET fix
  auto jet_mask = ROOT::VecOps::RVec<bool>(jet_pt.size(), true);
  auto lowptjet_mask = ROOT::VecOps::RVec<bool>(lowptjet_rawpt.size(), true);
  LogDebug << "JME:: First the (vetoed) jets in the noisy region" << std::endl;
  const auto offset_jets = calculateFixEE2017Offset(jet_mask,
      jet_pt, jet_eta, jet_phi, jet_mass,
      jet_rawcorr, jet_area, jet_muonSubtrFactor,
      rho);
  const auto offset_lowptjets = calculateFixEE2017Offset(lowptjet_mask,
      lowptjet_rawpt, lowptjet_eta, lowptjet_phi, lowptjet_zero,
      lowptjet_zero, lowptjet_area, lowptjet_muonSubtrFactor,
      rho);
  const auto delta_x_T1Jet  = offset_jets[0]+offset_lowptjets[0];
  const auto delta_y_T1Jet  = offset_jets[1]+offset_lowptjets[1];
  const auto delta_x_RawJet = offset_jets[2]+offset_lowptjets[2];
  const auto delta_y_RawJet = offset_jets[3]+offset_lowptjets[3];
  // default MET : add delta_T1Jet
  // unclustered EE : default MET (with delta T1Jet) - T1MET
  // correction for all MET variations: add delta_RawJet - unclustered EE
  const auto dx = delta_x_RawJet - ( defmet_pt*std::cos(defmet_phi) + delta_x_T1Jet - t1met_pt*std::cos(t1met_phi) );
  const auto dy = delta_y_RawJet - ( defmet_pt*std::sin(defmet_phi) + delta_y_T1Jet - t1met_pt*std::sin(t1met_phi) );
#ifdef BAMBOO_JME_DEBUG
  LogDebug << "JME:: T1      MET px=" << t1met_pt*std::cos(t1met_phi) << " py=" << t1met_pt*std::sin(t1met_phi) << std::endl;
  LogDebug << "JME:: default MET px=" << defmet_pt*std::cos(defmet_phi) << " py=" << defmet_pt*std::sin(defmet_phi) << std::endl;
  LogDebug << "JME:: raw     MET px=" << rawmet_pt*std::cos(rawmet_phi) << " py=" << rawmet_pt*std::sin(rawmet_phi) << std::endl;
  LogDebug << "JME:: deltas T1Jet x=" << delta_x_T1Jet << " y=" << delta_y_T1Jet << " RawJet x=" << delta_x_RawJet << " y= " << delta_y_RawJet << std::endl;
  LogDebug << "JME:: MET offset from jets in the noisy region: dx=" << dx << " and dy=" << dy << std::endl; // NO these are just minus the unclustered
#endif
  result_t out{nVariations, rawmet_pt*std::cos(rawmet_phi)+dx, rawmet_pt*std::sin(rawmet_phi)+dy};
  // usual variations, with jets that are in "unclustered EE" now vetoed
  auto& rg = getTRandom3(seed);
  // normal jets
  addVariations(out, jet_pt, jet_eta, jet_phi, jet_mass,
      jet_rawcorr, jet_area, jet_muonSubtrFactor, jet_neEmEF, jet_chEmEF, ROOT::VecOps::RVec<int>{},
      jet_mask, rho, genjet_pt, genjet_eta, genjet_phi, genjet_mass, rg);
  // low-PT jets
  addVariations(out, lowptjet_rawpt, lowptjet_eta, lowptjet_phi, lowptjet_zero,
      lowptjet_zero, lowptjet_area, lowptjet_muonSubtrFactor,
      ( lowptjet_neEmEF.empty() ? lowptjet_zero : lowptjet_neEmEF  ),
      ( lowptjet_chEmEF.empty() ? lowptjet_zero : lowptjet_chEmEF  ),
      ROOT::VecOps::RVec<int>{}, lowptjet_mask, rho,
      genjet_pt, genjet_eta, genjet_phi, genjet_mass, rg);
  // unclustered energy, based on nominal (0)
  out.setXY(nVariations-2, out.px(0)+met_unclustenupdx, out.py(0)+met_unclustenupdy);
  out.setXY(nVariations-1, out.px(0)-met_unclustenupdx, out.py(0)-met_unclustenupdy);

#ifdef BAMBOO_JME_DEBUG
  LogDebug << "JME:: returning " << out.size() << " modified METs" << std::endl;
  const auto varNames = available();
  assert(varNames.size() == nVariations);
  for ( std::size_t i{0}; i != nVariations; ++i ) {
    LogDebug << "JME:: MET_" << varNames[i] << ": PT=" << out.pt(i) << ", PHI=" << out.phi(i) << std::endl;
  }
#endif
  return out;
}

std::array<double,4> FixEE2017Type1METVariationsCalculator::calculateFixEE2017Offset(ROOT::VecOps::RVec<bool>& jet_mask,
    const p4compv_t& jet_pt, const p4compv_t& jet_eta, const p4compv_t& jet_phi, const p4compv_t& jet_mass,
    const p4compv_t& jet_rawcorr, const p4compv_t& jet_area, const p4compv_t& jet_muonSubtrFactor,
    const float rho
    ) const
{
  double delta_x_T1Jet{0.}, delta_y_T1Jet{0.};
  double delta_x_rawJet{0.}, delta_y_rawJet{0.};
  FactorizedJetCorrectorCalculator::VariableValues vals, valsL1;
  const auto nJets = jet_pt.size();
  for ( std::size_t i{0}; i != nJets; ++i ) {
    if ( ( 2.65 < std::abs(jet_eta[i]) ) && ( std::abs(jet_eta[i]) < 3.14 ) ) {
      const double jet_pt_raw = jet_pt[i]*(1-jet_rawcorr[i]);
      if ( jet_pt_raw < 50. ) {
        jet_mask[i] = false; // these are the jets to veto for the nominal variations
        // L1 and full (L1L2L3) JEC for muon-subtracted jet
        vals.setJetEta(jet_eta[i]);
        vals.setJetPt(jet_pt_raw);
        vals.setJetA(jet_area[i]);
        vals.setRho(rho);
        LogDebug << "Jet #" << i << " ETA=" << jet_eta[i] << ", PT_raw=" << jet_pt_raw << ", area=" << jet_area[i] << std::endl;
        auto jecL1L2L3 = ( m_jetCorrectorProd ? m_jetCorrectorProd : m_jetCorrector )->getCorrection(vals);
        if ( jecL1L2L3 <= 0. ) { jecL1L2L3 = 1.; }
        valsL1.setJetEta(jet_eta[i]);
        valsL1.setJetPt(jet_pt_raw);
        valsL1.setJetA(jet_area[i]);
        valsL1.setRho(rho);
        auto jecL1 = ( m_jetCorrectorL1Prod ? m_jetCorrectorL1Prod : m_jetCorrectorL1 )->getCorrection(valsL1);
        if ( jecL1     <= 0. ) { jecL1     = 1.; }
        const auto jet_pt_raw_nomu = jet_pt_raw*(1-jet_muonSubtrFactor[i]);
        const auto muon_pt = jet_pt_raw*jet_muonSubtrFactor[i];
        const auto jet_pt_nomuL1L2L3 = jet_pt_raw_nomu*jecL1L2L3;
        const auto jet_pt_nomuL1     = jet_pt_raw_nomu*jecL1;
        const auto jet_pt_L1L2L3 = jet_pt_nomuL1L2L3 + muon_pt;
        const auto jet_pt_L1     = jet_pt_nomuL1     + muon_pt;
        LogDebug << "jecL1L2L3=" << jecL1L2L3 << ", jecL1=" << jecL1 << "; PT_L1L2L3=" << jet_pt_L1L2L3 << ", PT_L1=" << jet_pt_L1 << ", PT_mu=" << muon_pt << std::endl;
        if ( jet_pt_nomuL1L2L3 > m_unclEnThreshold ) {
          const auto cosphi = std::cos(jet_phi[i]);
          const auto sinphi = std::sin(jet_phi[i]);
          delta_x_T1Jet += (jet_pt_L1L2L3-jet_pt_L1+jet_pt_raw)*cosphi;
          delta_y_T1Jet += (jet_pt_L1L2L3-jet_pt_L1+jet_pt_raw)*sinphi;
          delta_x_rawJet += jet_pt_raw*cosphi;
          delta_y_rawJet += jet_pt_raw*sinphi;
        }
      }
    }
  }
  return {{ delta_x_T1Jet, delta_y_T1Jet, delta_x_rawJet, delta_y_rawJet }};
}
