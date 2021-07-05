#ifndef UTILITY_H
#define UTILITY_H
#include <map>
#include <ROOT/RVec.hxx>
#include "TRandom3.h"
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

namespace {
  TRandom3& getTRandom3(uint32_t seed) {
    static thread_local TRandom3 rg{};
    rg.SetSeed(seed);
    return rg;
  }

  double jetESmearFactor( double pt, double eOrig, float genPt, float ptRes, float sfUncert, double rand )
  {
    double smear = 1.;
    if ( genPt > 0. ) {
      smear = 1. + (sfUncert-1.)*(pt - genPt)/pt;
    } else if ( sfUncert > 1. ) {
      smear = 1. + rand*std::sqrt(sfUncert*sfUncert-1.);
    }
    if ( smear*eOrig < 1.e-2 ) {
      smear = 1.e-2/eOrig;
    }
    return smear;
  }
  
  // because something goes wrong with linking ROOT::Math::VectorUtil::Phi_mpi_pi
  template<typename T>
  T phi_mpi_pi(T angle) {
    if ( angle <= M_PI && angle > -M_PI ) {
      return angle;
    }
    if ( angle > 0 ) {
      const int n = static_cast<int>(.5*(angle*M_1_PI+1.));
      angle -= 2*n*M_PI;
    } else {
      const int n = static_cast<int>(-.5*(angle*M_1_PI-1.));
      angle += 2*n*M_PI;
    }
    return angle;
  }

    std::vector<float> fillVector(const std::vector<std::string>& names, const JME::JetParameters::value_type& jetParams)
  {
    static const std::unordered_map<std::string,JME::Binning> jmeBinningFromString = {
        {"JetEta", JME::Binning::JetEta},
        {"JetPt" , JME::Binning::JetPt},
        // TODO JetPhi, JetEMF, LepPx, LepPy, LepPz                                                                                                                                    
        {"JetE"  , JME::Binning::JetE}
      };
    std::vector<float> result;
    result.reserve(names.size());
    for ( const auto& nm : names ) {
      const auto it_key = jmeBinningFromString.find(nm);
      if ( std::end(jmeBinningFromString) == it_key ) {
        throw std::runtime_error{"Unknown binning variable: "+nm};
      } else {
        const auto it_par = jetParams.find(it_key->second);
        if ( std::end(jetParams) == it_par ) {
          throw std::runtime_error{"Binning variable "+nm+" not found"};
        } else {
          result.push_back(it_par->second);
        }
      }
    }
    return result;
  }

  float getUncertainty(const SimpleJetCorrectionUncertainty& uncert, const JME::JetParameters::value_type& jetParams, bool direction)
  {
    const auto vx = fillVector(uncert.parameters().definitions().binVar(), jetParams);
    const auto vy = fillVector(uncert.parameters().definitions().parVar(), jetParams);
    return uncert.uncertainty(vx, vy[0], direction);
  }

  float deltaHEM2018Issue(float pt_nom, int jetId, float phi, float eta ) {
    float delta = 1.;
    if ( pt_nom > 15. && ( jetId & 0x2 ) && phi > -1.57 && phi < -0.87 ) {
      if ( eta > -2.5 && eta < -1.3 ) {
        delta = 0.8;
      } else if ( eta <= -2.5 && eta > -3. ) {
        delta = 0.65;
      }
    }
    return delta;
  }

    int jerSplitID(float pt, float eta) {
    const auto aEta = std::abs(eta);
    if ( aEta < 1.93 )
      return 0;
    else if ( aEta < 2.5 )
      return 1;
    else if ( aEta < 3. )
      if ( pt < 50. )
        return 2;
      else
        return 3;
    else
      if ( pt < 50. )
        return 4;
      else
        return 5;
  }
}

  

#endif
