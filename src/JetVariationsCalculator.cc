#if defined(PROJECT_NAME) && defined(CMSSW_GIT_HASH)
#include "UserCode/CMSJMECalculators/interface/JetVariationsCalculator.h"
#else
#include "JetVariationsCalculator.h"
#endif
#include <Math/GenVector/LorentzVector.h>
#include <Math/GenVector/PtEtaPhiM4D.h>
#include "Math/VectorUtil.h"
#include "TRandom3.h"
#include <cassert>

JetVariationsCalculator::result_t JetVariationsCalculator::produce(
    const p4compv_t& jet_pt, const p4compv_t& jet_eta, const p4compv_t& jet_phi, const p4compv_t& jet_mass,
    const p4compv_t& jet_rawcorr, const p4compv_t& jet_area, const ROOT::VecOps::RVec<int>& jet_jetId, const float rho,
    const std::uint32_t seed,
    const p4compv_t& genjet_pt, const p4compv_t& genjet_eta, const p4compv_t& genjet_phi, const p4compv_t& genjet_mass ) const
{
  const auto nVariations = 1+( m_doSmearing ? 2*( m_splitJER ? 6 : 1 ) : 0 )+2*m_jesUncSources.size()+( m_addHEM2018Issue ? 2 : 0 ); // 1(nom)+2(JER)+2*len(JES)[+2(HEM)]
  LogDebug << "JME:: hello from JetVariations produce. Got " << jet_pt.size() << " jets" << std::endl;
  const auto nJets = jet_pt.size();
  result_t out{nVariations, jet_pt, jet_mass};
  ROOT::VecOps::RVec<double> pt_nom{jet_pt}, mass_nom{jet_mass};
  if ( m_jetCorrector ) {
    LogDebug << "JME:: reapplying JEC" << std::endl;
    FactorizedJetCorrectorCalculator::VariableValues vals;
    for ( std::size_t i{0}; i != nJets; ++i ) {
      vals.setJetEta(jet_eta[i]);
      vals.setJetPt(jet_pt[i]*(1.-jet_rawcorr[i]));
      vals.setJetA(jet_area[i]);
      vals.setRho(rho);
      const auto corr = m_jetCorrector->getCorrection(vals);
      if ( corr > 0. ) {
        const double newc = (1.-jet_rawcorr[i])*corr;
        pt_nom[i]   *= newc;
        mass_nom[i] *= newc;
      }
    }
#ifdef BAMBOO_JME_DEBUG
    LogDebug << "JME:: with reapplied JEC: ";
    for ( std::size_t i{0}; i != nJets; ++i ) {
      LogDebug << "(PT=" << pt_nom[i] << ", ETA=" << jet_eta[i] << ", PHI=" << jet_phi[i] << ", M=" << mass_nom[i] << ") ";
    }
    LogDebug << std::endl;
#endif
  } else {
    LogDebug << "JME:: Not reapplying JEC" << std::endl;
  }
  // smearing and JER
  std::size_t iVar = 1; // after nominal
  if ( m_doSmearing ) {
    LogDebug << "JME:: Smearing (seed=" << seed << ")" << std::endl;
    auto& rg = getTRandom3(seed);
    p4compv_t pt_jerUp(pt_nom.size(), 0.), mass_jerUp(mass_nom.size(), 0.);
    p4compv_t pt_jerDown(pt_nom.size(), 0.), mass_jerDown(mass_nom.size(), 0.);
    for ( std::size_t i{0}; i != nJets; ++i ) {
      const auto eOrig = ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>(pt_nom[i], jet_phi[i], jet_eta[i], mass_nom[i]).E();
      double smearFactor_nom{1.}, smearFactor_down{1.}, smearFactor_up{1.};
      if ( pt_nom[i] > 0. ) {
        JME::JetParameters jPar{
            {JME::Binning::JetPt , pt_nom[i]},
            {JME::Binning::JetEta, jet_eta[i]},
            {JME::Binning::Rho   , rho} };
        const auto ptRes  = m_jetPtRes.getResolution(jPar);
        LogDebug << "JME:: JetParameters: pt=" << pt_nom[i] << ", eta=" << jet_eta[i] << ", rho=" << rho << "; ptRes=" << ptRes << std::endl;
        LogDebug << "JME:: ";
        float genPt = -1;
        if ( m_smearDoGenMatch ) {
          const auto iGen = findGenMatch(pt_nom[i], jet_eta[i], jet_phi[i], genjet_pt, genjet_eta, genjet_phi, ptRes*pt_nom[i]);
          if ( iGen != genjet_pt.size() ) {
            genPt = genjet_pt[iGen];
            LogDebug << "genPt=" << genPt << " ";
          }
        }
        const auto rand = ( genPt < 0. ) ? rg.Gaus(0, ptRes) : -1.;
        LogDebug << "jet_pt_resolution: " << ptRes << ", rand: " << rand << std::endl;
        smearFactor_nom  = jetESmearFactor(pt_nom[i], eOrig, genPt, ptRes, m_jetEResSF.getScaleFactor(jPar, Variation::NOMINAL), rand);
        smearFactor_down = jetESmearFactor(pt_nom[i], eOrig, genPt, ptRes, m_jetEResSF.getScaleFactor(jPar, Variation::DOWN   ), rand);
        smearFactor_up   = jetESmearFactor(pt_nom[i], eOrig, genPt, ptRes, m_jetEResSF.getScaleFactor(jPar, Variation::UP     ), rand);
        // LogDebug << "  scalefactors are NOMINAL=" << m_jetEResSF.getScaleFactor(jPar, Variation::NOMINAL) << ", DOWN=" << m_jetEResSF.getScaleFactor(jPar, Variation::DOWN) << ", UP=" << m_jetEResSF.getScaleFactor(jPar, Variation::UP) << std::endl;
        // LogDebug << "  smearfactors are NOMINAL=" << smearFactor_nom << ", DOWN=" << smearFactor_down << ", UP=" << smearFactor_up << std::endl;
      }
      pt_jerDown[i]   = pt_nom[i]*smearFactor_down;
      mass_jerDown[i] = mass_nom[i]*smearFactor_down;
      pt_jerUp[i]     = pt_nom[i]*smearFactor_up;
      mass_jerUp[i]   = mass_nom[i]*smearFactor_up;
      pt_nom[i]       *= smearFactor_nom;
      mass_nom[i]     *= smearFactor_nom;
    }
    if ( m_splitJER ) {
      ROOT::VecOps::RVec<int> jerBin(pt_nom.size(), -1);
      for ( std::size_t j{0}; j != nJets; ++j ) {
        jerBin[j] = jerSplitID(pt_nom[j], jet_eta[j]);
      }
      for ( int i{0}; i != 6; ++i ) {
        p4compv_t pt_jeriUp{pt_nom}, mass_jeriUp{mass_nom};
        p4compv_t pt_jeriDown{pt_nom}, mass_jeriDown{mass_nom};
        for ( std::size_t j{0}; j != nJets; ++j ) {
          if ( jerBin[j] == i ) {
            pt_jeriUp[j] = pt_jerUp[j];
            pt_jeriDown[j] = pt_jerDown[j];
            mass_jeriUp[j] = mass_jerUp[j];
            mass_jeriDown[j] = mass_jerDown[j];
          }
        }
        out.set(iVar++, std::move(pt_jeriUp)  , std::move(mass_jeriUp)  );
        out.set(iVar++, std::move(pt_jeriDown), std::move(mass_jeriDown));
      }
    } else {
      out.set(iVar++, std::move(pt_jerUp)  , std::move(mass_jerUp)  );
      out.set(iVar++, std::move(pt_jerDown), std::move(mass_jerDown));
    }
    LogDebug << "JME:: Done with smearing" << std::endl;
  } else {
    LogDebug << "JME:: No smearing" << std::endl;
  }
  out.set(0, pt_nom, mass_nom);

  // HEM issue 2018, see https://hypernews.cern.ch/HyperNews/CMS/get/JetMET/2000.html
  if ( m_addHEM2018Issue ) {
    p4compv_t pt_down(pt_nom.size(), 0.), mass_down(mass_nom.size(), 0.);
    for ( std::size_t j{0}; j != nJets; ++j ) {
      const auto delta = deltaHEM2018Issue(pt_nom[j], jet_jetId[j], jet_phi[j], jet_eta[j]);
      pt_down[j] = pt_nom[j]*delta;
      mass_down[j] = mass_nom[j]*delta;
    }
    out.set(iVar++, pt_nom, mass_nom);
    out.set(iVar++, std::move(pt_down), std::move(mass_down));
  }
  // JES uncertainties
  for ( auto& jesUnc : m_jesUncSources ) {
    LogDebug << "JME:: evaluating JES uncertainty: " << jesUnc.first << std::endl;
    p4compv_t pt_jesDown(pt_nom.size(), 0.), mass_jesDown(mass_nom.size(), 0.);
    p4compv_t pt_jesUp(pt_nom.size(), 0.), mass_jesUp(mass_nom.size(), 0.);
    for ( std::size_t i{0}; i != nJets; ++i ) {
      const auto delta = getUncertainty(jesUnc.second, { {JME::Binning::JetPt, pt_nom[i]}, {JME::Binning::JetEta, jet_eta[i]} }, true);
      pt_jesDown[i]   = pt_nom[i]*(1.-delta);
      mass_jesDown[i] = mass_nom[i]*(1.-delta);
      pt_jesUp[i]     = pt_nom[i]*(1.+delta);
      mass_jesUp[i]   = mass_nom[i]*(1.+delta);
    }
    out.set(iVar++, std::move(pt_jesUp)  , std::move(mass_jesUp)  );
    out.set(iVar++, std::move(pt_jesDown), std::move(mass_jesDown));
  }

#ifdef BAMBOO_JME_DEBUG
  assert(iVar == out.size());
  LogDebug << "JME:: returning " << out.size() << " modified jet collections" << std::endl;
  const auto varNames = available();
  assert(varNames.size() == nVariations);
  for ( std::size_t i{0}; i != nVariations; ++i ) {
    LogDebug << "JME:: Jet_" << varNames[i] << ": ";
    for ( std::size_t j{0}; j != nJets; ++j ) {
      LogDebug << "(PT=" << out.pt(i)[j] << ", ETA=" << jet_eta[j] << ", PHI=" << jet_phi[j] << ", M=" << out.mass(i)[j] << ") ";
    }
    LogDebug << std::endl;
  }
#endif
  return out;
}

std::vector<std::string> JetVariationsCalculator::available(const std::string&) const
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
  return products;
}
