#if defined(PROJECT_NAME) && defined(CMSSW_GIT_HASH)
#include "UserCode/CMSJMECalculators/interface/Type1METVariationsCalculator.h"
#else
#include "Type1METVariationsCalculator.h"
#endif
#include <Math/GenVector/LorentzVector.h>
#include <Math/GenVector/PtEtaPhiM4D.h>
#include "Math/VectorUtil.h"
#include "TRandom3.h"
#include <cassert>

void Type1METVariationsCalculator::setL1JEC(const std::vector<JetCorrectorParameters>& jecParams)
{
  if ( ! jecParams.empty() ) {
    //m_jetCorrectorL1 = std::unique_ptr<FactorizedJetCorrectorCalculator,jetcorrdeleter>{new FactorizedJetCorrectorCalculator(jecParams)};
    m_jetCorrectorL1 = std::unique_ptr<FactorizedJetCorrectorCalculator>{new FactorizedJetCorrectorCalculator(jecParams)};
  }
}

Type1METVariationsCalculator::result_t Type1METVariationsCalculator::produce(
    const p4compv_t& jet_pt, const p4compv_t& jet_eta, const p4compv_t& jet_phi, const p4compv_t& jet_mass,
    const p4compv_t& jet_rawcorr, const p4compv_t& jet_area,
    const p4compv_t& jet_muonSubtrFactor, const p4compv_t& jet_neEmEF, const p4compv_t& jet_chEmEF, const ROOT::VecOps::RVec<int>& jet_jetId,
    const float rho,
    const std::uint32_t seed,
    const p4compv_t& genjet_pt, const p4compv_t& genjet_eta, const p4compv_t& genjet_phi, const p4compv_t& genjet_mass,
    const float rawmet_phi, const float rawmet_pt,
    const float met_unclustenupdx, const float met_unclustenupdy,
    const p4compv_t& lowptjet_rawpt, const p4compv_t& lowptjet_eta, const p4compv_t& lowptjet_phi, const p4compv_t& lowptjet_area,
    const p4compv_t& lowptjet_muonSubtrFactor, const p4compv_t& lowptjet_neEmEF, const p4compv_t& lowptjet_chEmEF
    ) const
{
  const auto nVariations = 3+( m_doSmearing ? 2*( m_splitJER ? 6 : 1 ) : 0 )+2*m_jesUncSources.size()+( m_addHEM2018Issue ? 2 : 0 ); // 1(nom)+2(unclust)+2(JER[*6])+2*len(JES)[+2(HEM)]
  result_t out{nVariations, rawmet_pt*std::cos(rawmet_phi), rawmet_pt*std::sin(rawmet_phi)};
  auto& rg = getTRandom3(seed);
  LogDebug << "JME:: hello from Type1METVariations produce. Got " << jet_pt.size() << " jets and " << lowptjet_rawpt.size() << " low-PT jets" << std::endl;
  // normal jets
  addVariations(out, jet_pt, jet_eta, jet_phi, jet_mass,
      jet_rawcorr, jet_area, jet_muonSubtrFactor, jet_neEmEF, jet_chEmEF, jet_jetId,
      ROOT::VecOps::RVec<bool>(), rho, genjet_pt, genjet_eta, genjet_phi, genjet_mass, rg);
  // low-PT jets
  p4compv_t lowptjet_zero(lowptjet_rawpt.size(), 0.);
  addVariations(out, lowptjet_rawpt, lowptjet_eta, lowptjet_phi, lowptjet_zero,
      lowptjet_zero, lowptjet_area, lowptjet_muonSubtrFactor,
      ( lowptjet_neEmEF.empty() ? lowptjet_zero : lowptjet_neEmEF  ),
      ( lowptjet_chEmEF.empty() ? lowptjet_zero : lowptjet_chEmEF  ),
      ROOT::VecOps::RVec<int>(lowptjet_rawpt.size(), 0), ROOT::VecOps::RVec<bool>(), rho,
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

// for a single jet collection
void Type1METVariationsCalculator::addVariations(Type1METVariationsCalculator::result_t& out,
    const p4compv_t& jet_pt, const p4compv_t& jet_eta, const p4compv_t& jet_phi, const p4compv_t& jet_mass,
    const p4compv_t& jet_rawcorr, const p4compv_t& jet_area, const p4compv_t& jet_muonSubtrFactor,
    const p4compv_t& jet_neEmEF, const p4compv_t& jet_chEmEF, const ROOT::VecOps::RVec<int>& jet_jetId, const ROOT::VecOps::RVec<bool>& jet_mask, const float rho,
    const p4compv_t& genjet_pt, const p4compv_t& genjet_eta, const p4compv_t& genjet_phi, const p4compv_t& genjet_mass,
    TRandom3& rg) const
{
  FactorizedJetCorrectorCalculator::VariableValues vals, valsL1;
  const auto nJets = jet_pt.size();
  for ( std::size_t i{0}; i != nJets; ++i ) {
    // L1 and full (L1L2L3) JEC for muon-subtracted jet
    double jet_pt_nom = jet_pt[i];
    double jet_mass_nom = jet_mass[i];
    const auto jet_pt_raw = jet_pt_nom*(1-jet_rawcorr[i]);
    vals.setJetEta(jet_eta[i]);
    vals.setJetPt(jet_pt_raw);
    vals.setJetA(jet_area[i]);
    vals.setRho(rho);
    LogDebug << "Jet #" << i << " ETA=" << jet_eta[i] << ", PT_raw=" << jet_pt_raw << ", area=" << jet_area[i] << std::endl;
    auto jecL1L2L3 = m_jetCorrector->getCorrection(vals);
    if ( jecL1L2L3 <= 0. ) {
      jecL1L2L3 = 1.;
    } else {
      jet_pt_nom = jet_pt_raw*jecL1L2L3;
      jet_mass_nom = jet_mass_nom*(1-jet_rawcorr[i])*jecL1L2L3;
    }
    valsL1.setJetEta(jet_eta[i]);
    valsL1.setJetPt(jet_pt_raw);
    valsL1.setJetA(jet_area[i]);
    valsL1.setRho(rho);
    auto jecL1 = m_jetCorrectorL1->getCorrection(valsL1);
    if ( jecL1     <= 0. ) { jecL1     = 1.; }
    const double jet_pt_raw_nomu = jet_pt_raw*(1-jet_muonSubtrFactor[i]);
    const double muon_pt = jet_pt_raw*jet_muonSubtrFactor[i];
    const auto jet_pt_nomuL1L2L3 = jet_pt_raw_nomu*jecL1L2L3;
    const auto jet_pt_nomuL1     = jet_pt_raw_nomu*jecL1;
    const auto jet_pt_L1L2L3 = jet_pt_nomuL1L2L3 + muon_pt;
    const auto jet_pt_L1     = jet_pt_nomuL1     + muon_pt;
    LogDebug << "jecL1L2L3=" << jecL1L2L3 << ", jecL1=" << jecL1 << "; PT_L1L2L3=" << jet_pt_L1L2L3 << ", PT_L1=" << jet_pt_L1 << ", PT_mu=" << muon_pt << std::endl;
    // JER / smearing
    double jerF_nom{1.}, jerF_up{1.}, jerF_down{1.};
    if ( m_doSmearing ) {
      const auto eOrig = ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>(jet_pt_nom, jet_phi[i], jet_eta[i], jet_mass_nom).E();
      if ( jet_pt_nom > 0. ) {
        JME::JetParameters jPar{
            {JME::Binning::JetPt , jet_pt_nom},
            {JME::Binning::JetEta, jet_eta[i]},
            {JME::Binning::Rho   , rho} };
        const auto ptRes  = m_jetPtRes.getResolution(jPar);
        LogDebug << "JME:: JetParameters: pt=" << jet_pt_nom << ", eta=" << jet_eta[i] << ", rho=" << rho << "; ptRes=" << ptRes << std::endl;
        LogDebug << "JME:: ";
        float genPt = -1;
        if ( m_smearDoGenMatch ) {
          const auto iGen = findGenMatch(jet_pt_nom, jet_eta[i], jet_phi[i], genjet_pt, genjet_eta, genjet_phi, ptRes*jet_pt_nom);
          if ( iGen != genjet_pt.size() ) {
            genPt = genjet_pt[iGen];
            LogDebug << "genPt=" << genPt << " ";
          }
        }
        const auto rand = ( genPt < 0. ) ? rg.Gaus(0, ptRes) : -1.;
        LogDebug << "jet_pt_resolution: " << ptRes << ", rand: " << rand << std::endl;
        jerF_nom  = jetESmearFactor(jet_pt_nom, eOrig, genPt, ptRes, m_jetEResSF.getScaleFactor(jPar, Variation::NOMINAL), rand);
        jerF_down = jetESmearFactor(jet_pt_nom, eOrig, genPt, ptRes, m_jetEResSF.getScaleFactor(jPar, Variation::DOWN   ), rand);
        jerF_up   = jetESmearFactor(jet_pt_nom, eOrig, genPt, ptRes, m_jetEResSF.getScaleFactor(jPar, Variation::UP     ), rand);
        // LogDebug << "  scalefactors are NOMINAL=" << m_jetEResSF.getScaleFactor(jPar, Variation::NOMINAL) << ", DOWN=" << m_jetEResSF.getScaleFactor(jPar, Variation::DOWN) << ", UP=" << m_jetEResSF.getScaleFactor(jPar, Variation::UP) << std::endl;
        // LogDebug << "  smearfactors are NOMINAL=" << jerF_nom << ", DOWN=" << jerF_down << ", UP=" << jerF_up << std::endl;
        jet_pt_nom *= jerF_nom; // for consistency with jet case (used for split JER and HEM issue below)
      }
    }
    if ( ( jet_mask.empty() || jet_mask[i] ) && ( jet_pt_nomuL1L2L3 > m_unclEnThreshold ) && ( (jet_neEmEF[i]+jet_chEmEF[i]) < 0.9 ) ) {
      std::size_t iVar = 0;
      const auto jet_cosPhi = std::cos(jet_phi[i]);
      const auto jet_sinPhi = std::sin(jet_phi[i]);
      if ( ! ( m_doSmearing && m_isT1SmearedMET ) ) {
        out.addR_proj(iVar++, jet_cosPhi, jet_sinPhi, jet_pt_L1 - jet_pt_L1L2L3);           // nominal
      }
      auto jet_pt_L1p = jet_pt_L1; // with optional offset for JES uncertainty calculation if nominal is smeared
      if ( m_doSmearing ) {
        const auto dr_jernom  = jet_pt_L1 - jet_pt_L1L2L3*jerF_nom;
        if ( m_isT1SmearedMET ) {
          const auto dr_jerup   = jet_pt_L1 - jet_pt_L1L2L3*jerF_up;
          const auto dr_jerdown = jet_pt_L1 - jet_pt_L1L2L3*jerF_down;
          out.addR_proj(iVar++, jet_cosPhi, jet_sinPhi, dr_jernom);                         // smeared nominal
          if ( m_splitJER ) {
            const auto jerBin = jerSplitID(jet_pt_nom, jet_eta[i]);
            for ( int k{0}; k != 6; ++k ) {
              if ( jerBin == k ) { // vary
                out.addR_proj(iVar++, jet_cosPhi, jet_sinPhi, dr_jerup);                    // JER[k]-up
                out.addR_proj(iVar++, jet_cosPhi, jet_sinPhi, dr_jerdown);                  // JER[k]-down
              } else { // keep nominal
                out.addR_proj(iVar++, jet_cosPhi, jet_sinPhi, dr_jernom);                   // JER[k]-up
                out.addR_proj(iVar++, jet_cosPhi, jet_sinPhi, dr_jernom);                   // JER[k]-down
              }
            }
          } else {
            out.addR_proj(iVar++, jet_cosPhi, jet_sinPhi, dr_jerup);                        // JER-up
            out.addR_proj(iVar++, jet_cosPhi, jet_sinPhi, dr_jerdown);                      // JER-down
          }
          jet_pt_L1p += jet_pt_L1L2L3*(1.-jerF_nom); // offset for JES uncertainties, since the nominal is smeared
        } else {
          for ( std::size_t k{0}; k != ( m_splitJER ? 6 : 1 ); ++k ) {
            out.addR_proj(iVar++, jet_cosPhi, jet_sinPhi, dr_jernom);                     // JER[k]-up
            out.addR_proj(iVar++, jet_cosPhi, jet_sinPhi, dr_jernom);                     // JER[k]-down
          }
        }
      }
      if ( m_addHEM2018Issue ) {
        const auto delta = deltaHEM2018Issue(jet_pt_nom, jet_jetId[i], jet_phi[i], jet_eta[i]);
        out.addR_proj(iVar++, jet_cosPhi, jet_sinPhi, jet_pt_L1p - jet_pt_L1L2L3);           // up = nominal
        out.addR_proj(iVar++, jet_cosPhi, jet_sinPhi, jet_pt_L1p - jet_pt_L1L2L3*delta);     // down
      }
      for ( auto& jesUnc : m_jesUncSources ) {
        const auto delta = getUncertainty(jesUnc.second, { {JME::Binning::JetPt, jet_pt_L1L2L3}, {JME::Binning::JetEta, jet_eta[i]} }, true);
        out.addR_proj(iVar++, jet_cosPhi, jet_sinPhi, jet_pt_L1p - jet_pt_L1L2L3*(1+delta)); // JES_i-up
        out.addR_proj(iVar++, jet_cosPhi, jet_sinPhi, jet_pt_L1p - jet_pt_L1L2L3*(1-delta)); // JES_i-down
      }
#ifdef BAMBOO_JME_DEBUG
      assert(iVar+2 == out.size()); // last two are unclustered energy up and down
#endif
    }
  }
}

std::vector<std::string> Type1METVariationsCalculator::available(const std::string&) const
{
  std::vector<std::string> products = { "nominal" };
  if ( m_doSmearing ) {
    if ( m_splitJER ) {
      for ( int i = 0; i != 6; ++i ) {
        products.emplace_back("jer"+std::to_string(i)+"up");
        products.emplace_back("jer"+std::to_string(i)+"down");
      }
    } else {
      products.emplace_back("jerup");
      products.emplace_back("jerdown");
    }
  }
  if ( m_addHEM2018Issue ) {
    products.emplace_back("jesHEMIssueup");
    products.emplace_back("jesHEMIssuedown");
  }
  for ( const auto& src : m_jesUncSources ) {
    products.emplace_back("jes"+src.first+"up");
    products.emplace_back("jes"+src.first+"down");
  }
  products.emplace_back("unclustEnup");
  products.emplace_back("unclustEndown");
  return products;
}

