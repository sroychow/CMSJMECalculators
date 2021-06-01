import pytest
import os, os.path

testData = os.path.join(os.path.dirname(__file__), "data")

if "CMSSW_BASE" in os.environ:  # CMSSW/scram version
    from UserCode.CMSJMECalculators.CMSJMECalculators.utils import (
            loadJMESystematicsCalculators,
            fatjet_jmsValues,
            fatjet_jmrValues,
            fatjet_gmsValues_tau21DDT,
            fatjet_gmrValues_tau21DDT,
            fatjet_puppi_msd_params
            )
    from UserCode.CMSJMECalculators.CMSJMECalculators.jetdatabasecache import JetDatabaseCache
else:  # pip version
    from CMSJMECalculators.utils import (
            loadJMESystematicsCalculators,
            fatjet_jmsValues,
            fatjet_jmrValues,
            fatjet_gmsValues_tau21DDT,
            fatjet_gmrValues_tau21DDT,
            fatjet_puppi_msd_params
            )
    from CMSJMECalculators.jetdatabasecache import JetDatabaseCache

jecDBCache = JetDatabaseCache("JECDatabase", repository="cms-jet/JECDatabase", cachedir=testData, mayWrite=True)
jrDBCache = JetDatabaseCache("JRDatabase", repository="cms-jet/JRDatabase", cachedir=testData, mayWrite=True)

def toRVecFloat(values):
    import ROOT as gbl
    RVec_float = getattr(gbl, "ROOT::VecOps::RVec<float>")
    res = RVec_float(len(values), 0.)
    for i,val in enumerate(values):
        res[i] = val
    return res

def toRVecInt(values):
    import ROOT as gbl
    RVec_int = getattr(gbl, "ROOT::VecOps::RVec<int>")
    res = RVec_int(len(values), 0)
    for i,val in enumerate(values):
        res[i] = val
    return res

def getEventWith(f, condition=lambda ev : True, treeName="Events"):
    tup = f.Get(treeName)
    tup.GetEntry(0)
    i = 0
    while not condition(tup):
        i += 1
        tup.GetEntry(i)
    yield tup

def getJetMETArgs(tup, isMC=True, forMET=False, isMETFixEE2017=False, addHEM2018Issue=False):
    args = [
        toRVecFloat(tup.Jet_pt),
        toRVecFloat(tup.Jet_eta),
        toRVecFloat(tup.Jet_phi),
        toRVecFloat(tup.Jet_mass),
        toRVecFloat(tup.Jet_rawFactor),
        toRVecFloat(tup.Jet_area),
        ]
    if forMET:
        args += [
            toRVecFloat(tup.Jet_muonSubtrFactor),
            toRVecFloat(tup.Jet_neEmEF),
            toRVecFloat(tup.Jet_chEmEF),
            ]
    if not ( forMET and isMETFixEE2017 ):
        args.append(toRVecInt(tup.Jet_jetId if addHEM2018Issue else []))
    args.append(tup.fixedGridRhoFastjetAll)
    if isMC:
        args += [
            (tup.run<<20) + (tup.luminosityBlock<<10) + tup.event + 1 + ( int(tup.Jet_eta[0]/.01) if tup.nJet != 0 else 0),
            toRVecFloat(tup.GenJet_pt),
            toRVecFloat(tup.GenJet_eta),
            toRVecFloat(tup.GenJet_phi),
            toRVecFloat(tup.GenJet_mass)
            ]
    else:
        args += [ 0, toRVecFloat([]), toRVecFloat([]), toRVecFloat([]), toRVecFloat([]) ]
    if forMET:
        args += [ tup.RawMET_phi, tup.RawMET_pt ]
        if not isMETFixEE2017:
            args += [ tup.MET_MetUnclustEnUpDeltaX, tup.MET_MetUnclustEnUpDeltaY ]
        else:
            args += [ tup.METFixEE2017_MetUnclustEnUpDeltaX, tup.METFixEE2017_MetUnclustEnUpDeltaY ]
        args += [ toRVecFloat(getattr(tup, f"CorrT1METJet_{varNm}"))
                for varNm in ("rawPt", "eta", "phi", "area", "muonSubtrFactor") ]
        args += [ toRVecFloat([]), toRVecFloat([]) ]
        if isMETFixEE2017:
            args += [ tup.MET_phi, tup.MET_pt, tup.METFixEE2017_phi, tup.METFixEE2017_pt ]
    return args

def getFatJetArgs(tup, isMC=True, addHEM2018Issue=False):
    args = [
        toRVecFloat(tup.FatJet_pt),
        toRVecFloat(tup.FatJet_eta),
        toRVecFloat(tup.FatJet_phi),
        toRVecFloat(tup.FatJet_mass),
        toRVecFloat(tup.FatJet_rawFactor),
        toRVecFloat(tup.FatJet_area),
        toRVecFloat(tup.FatJet_msoftdrop),
        toRVecInt(tup.FatJet_subJetIdx1),
        toRVecInt(tup.FatJet_subJetIdx2),
        toRVecFloat(tup.SubJet_pt),
        toRVecFloat(tup.SubJet_eta),
        toRVecFloat(tup.SubJet_phi),
        toRVecFloat(tup.SubJet_mass),
        toRVecInt(tup.FatJet_jetId if addHEM2018Issue else []),
        tup.fixedGridRhoFastjetAll
        ]
    if isMC:
        args += [
            (tup.run<<20) + (tup.luminosityBlock<<10) + tup.event + 1 + ( int(tup.Jet_eta[0]/.01) if tup.nJet != 0 else 0),
            toRVecFloat(tup.GenJetAK8_pt),
            toRVecFloat(tup.GenJetAK8_eta),
            toRVecFloat(tup.GenJetAK8_phi),
            toRVecFloat(tup.GenJetAK8_mass),
            toRVecFloat(tup.SubGenJetAK8_pt),
            toRVecFloat(tup.SubGenJetAK8_eta),
            toRVecFloat(tup.SubGenJetAK8_phi),
            toRVecFloat(tup.SubGenJetAK8_mass)
            ]
    else:
        args += [ 0, toRVecFloat([]), toRVecFloat([]), toRVecFloat([]), toRVecFloat([]), toRVecFloat([]), toRVecFloat([]), toRVecFloat([]), toRVecFloat([]) ]
    return args

def configureCalc(calc, jecTag=None, jerTag=None, jetType="AK4PFchs", levels=None, levels_l1=None, splitJER=False, uncSources=None):
    import ROOT as gbl
    if jecTag:
        jcp = {}
        if levels:
            jecParams = getattr(gbl, "std::vector<JetCorrectorParameters>")()
            for iLev in levels:
                if iLev not in jcp:
                    jcp[iLev] = gbl.JetCorrectorParameters(jecDBCache.getPayload(jecTag, iLev, jetType))
                jecParams.push_back(jcp[iLev])
            calc.setJEC(jecParams)
        if uncSources:
            for jus in uncSources:
                param = gbl.JetCorrectorParameters(jecDBCache.getPayload(jecTag, "UncertaintySources", jetType), jus)
                calc.addJESUncertainty(jus, param)
        if levels_l1: ## for MET
            jecParams_l1 = getattr(gbl, "std::vector<JetCorrectorParameters>")()
            for iLev in levels_l1:
                if iLev not in jcp:
                    jcp[iLev] = gbl.JetCorrectorParameters(jecDBCache.getPayload(jecTag, iLev, jetType))
                jecParams_l1.push_back(jcp[iLev])
            calc.setL1JEC(jecParams_l1)
    if jerTag:
        calc.setSmearing(jrDBCache.getPayload(jerTag, "PtResolution", jetType), jrDBCache.getPayload(jerTag, "SF", jetType), splitJER, True, 0.2, 3.)

def configureFatJetCalc(calc, isMC=False, year=None, doSmearing=False, jmr=None, jms=None, isTau21DDT=False, gmr=None, gms=None, puppiGen=None, puppiRecoCorrCen=None, puppiRecoCorrFwd=None, puppiResolCen=None, puppiResolFwd=None):
    ## defaults
    if isMC: # mass scale
        if jms is None:
            jms = fatjet_jmsValues[year]
        if gms is None:
            if isTau21DDT:
                gms = fatjet_gmsValues_tau21DDT[year]
            else:
                gms = jms
    if jms is not None:
        print(f"Setting JMS values {jms}")
        if hasattr(jms, "__iter__"):
            calc.setJMSValues(*jms)
        else: # just nominal
            calc.setJMSValues(jms)
    if gms is not None:
        print(f"Setting GMS values {gms}")
        if hasattr(gms, "__iter__"):
            calc.setGMSValues(*jms)
        else: # just nominal
            calc.setGMSValues(jms)
    if doSmearing: # mass resolutions
        if jmr is None:
            jmr = fatjet_jmrValues[year]
        if gmr is None:
            if isTau21DDT:
                gmr = fatjet_gmrValues_tau21DDT[year]
            else:
                gmr = jmr
        print(f"Setting JMR values: {jmr}")
        calc.setJMRValues(*jmr)
        print(f"Setting GMR values: {gmr}")
        calc.setGMRValues(*gmr)
    ## PUPPI: always
    if puppiGen is None:
        puppiGen = fatjet_puppi_msd_params["gen"]
    if puppiRecoCorrCen is None:
        puppiRecoCorrCen = fatjet_puppi_msd_params["reco_cen"]
    if puppiRecoCorrFwd is None:
        puppiRecoCorrFwd = fatjet_puppi_msd_params["reco_fwd"]
    if puppiResolCen is None:
        puppiResolCen = fatjet_puppi_msd_params["resol_cen"]
    if puppiResolFwd is None:
        puppiResolFwd = fatjet_puppi_msd_params["resol_fwd"]
    from cppyy import gbl
    arrType = getattr(gbl, "std::array<double, 6>")
    def toArr(*args):
        assert len(args) == 6
        arr = arrType()
        for i,v in enumerate(args):
            arr[i] = v
        return arr
    print(f"Setting PUPPI gen formula to {puppiGen}")
    calc.setPuppiCorrections(puppiGen, *(toArr(*arg) for arg in (puppiRecoCorrCen, puppiRecoCorrFwd, puppiResolCen, puppiResolFwd)))

@pytest.fixture(scope="module")
def nanojetargsMC16():
    import ROOT as gbl
    f = gbl.TFile.Open(os.path.join(testData, "DY_M50_2016.root"))
    for tup in getEventWith(f, (lambda tup : tup.nJet >= 5)):
        yield getJetMETArgs(tup, isMC=True, forMET=False)

@pytest.fixture(scope="module")
def nanojetargsMC16_postvalues():
    import ROOT as gbl
    f = gbl.TFile.Open(os.path.join(testData, "DY_M50_2016postproc_JMEKin_bTagShape_puWeight.root"))
    tup = f.Get("Events")
    res = []
    for i in range(tup.GetEntries()):
        tup.GetEntry(i)
        jet_vars = {
            "nominal" : (toRVecFloat(tup.Jet_pt_nom    ), toRVecFloat(tup.Jet_mass_nom    )),
            }
        from itertools import chain
        jet_vars.update(dict(chain.from_iterable(
            { f"jer{i:d}{vdir.lower()}" : tuple(toRVecFloat(getattr(tup, f"Jet_{ivar}_jer{i:d}{vdir}")) for ivar in ("pt", "mass"))
                for vdir in ("Up", "Down") }.items() for i in range(6))))
        jet_vars.update(dict(chain.from_iterable(
            { f"jes{src}{vdir.lower()}" : tuple(toRVecFloat(getattr(tup, f"Jet_{ivar}_jes{src}{vdir}".format(ivar, src))) for ivar in ("pt", "mass"))
                for vdir in ("Up", "Down") }.items() for src in ("AbsoluteStat", "AbsoluteScale")
            )))
        res.append((getJetMETArgs(tup, isMC=True, forMET=False), jet_vars))
    yield res

@pytest.fixture(scope="module")
def nanoMETargsMC16_postvalues():
    import ROOT as gbl
    RVec_float = getattr(gbl, "ROOT::VecOps::RVec<float>")
    f = gbl.TFile.Open(os.path.join(testData, "DY_M50_2016postproc_JMEKin_bTagShape_puWeight.root"))
    tup = f.Get("Events")
    res = []
    for i in range(tup.GetEntries()):
        tup.GetEntry(i)
        met_vars = {
            "nominal" : (tup.MET_T1_pt , tup.MET_T1_phi ),
            "unclustEnup"   : (tup.MET_T1_pt_unclustEnUp  , tup.MET_T1_phi_unclustEnUp  ),
            "unclustEndown" : (tup.MET_T1_pt_unclustEnDown, tup.MET_T1_phi_unclustEnDown)
            }
        from itertools import chain
        met_vars.update(dict(
            ("{0}{1}".format(nm, var.lower()), (getattr(tup, "MET_T1_pt_{0}{1}".format(nm, var)), getattr(tup, "MET_T1_phi_{0}{1}".format(nm, var))))
            for var in ("Up", "Down") for nm in [f"jer{i:d}" for i in range(6)]+[ "jes{0}".format(jsnm) for jsnm in ("AbsoluteScale", "AbsoluteStat") ]
            ))
        res.append((tuple(getJetMETArgs(tup, isMC=True, forMET=True)), met_vars))
    yield res

@pytest.fixture(scope="module")
def nanoSmearMETargsMC16_postvalues():
    import ROOT as gbl
    RVec_float = getattr(gbl, "ROOT::VecOps::RVec<float>")
    f = gbl.TFile.Open(os.path.join(testData, "DY_M50_2016postproc_JMEKin_bTagShape_puWeight.root"))
    tup = f.Get("Events")
    res = []
    for i in range(tup.GetEntries()):
        tup.GetEntry(i)
        met_vars = {
            "nominal" : (tup.MET_T1Smear_pt , tup.MET_T1Smear_phi ),
            "unclustEnup"   : (tup.MET_T1Smear_pt_unclustEnUp  , tup.MET_T1Smear_phi_unclustEnUp  ),
            "unclustEndown" : (tup.MET_T1Smear_pt_unclustEnDown, tup.MET_T1Smear_phi_unclustEnDown)
            }
        from itertools import chain
        met_vars.update(dict(
            ("{0}{1}".format(nm, var.lower()), (getattr(tup, "MET_T1Smear_pt_{0}{1}".format(nm, var)), getattr(tup, "MET_T1Smear_phi_{0}{1}".format(nm, var))))
            for var in ("Up", "Down") for nm in [f"jer{i:d}" for i in range(6)]+[ "jes{0}".format(jsnm) for jsnm in ("AbsoluteScale", "AbsoluteStat") ]
            ))
        res.append((tuple(getJetMETArgs(tup, isMC=True, forMET=True)), met_vars))
    yield res

@pytest.fixture(scope="module")
def nanoMETFixEE2017argsMC17_postvalues():
    import ROOT as gbl
    RVec_float = getattr(gbl, "ROOT::VecOps::RVec<float>")
    f = gbl.TFile.Open(os.path.join(testData, "DY_M50_2017postproc_JMEKin_METFixEE2017.root"))
    tup = f.Get("Events")
    res = []
    for i in range(tup.GetEntries()):
        tup.GetEntry(i)
        met_vars = {
            "nominal" : (tup.METFixEE2017_T1_pt, tup.METFixEE2017_T1_phi),
            "unclustEnup"   : (tup.METFixEE2017_pt_unclustEnUp  , tup.METFixEE2017_phi_unclustEnUp  ),
            "unclustEndown" : (tup.METFixEE2017_pt_unclustEnDown, tup.METFixEE2017_phi_unclustEnDown)
            }
        from itertools import chain
        met_vars.update(dict(
            ("{0}{1}".format(nm, var.lower()), (getattr(tup, "METFixEE2017_T1_pt_{0}{1}".format(nm, var)), getattr(tup, "METFixEE2017_T1_phi_{0}{1}".format(nm, var))))
            for var in ("Up", "Down") for nm in ["jer"]+[ "jes{0}".format(jsnm) for jsnm in ("AbsoluteScale", "AbsoluteStat") ]
            ))
        res.append((tuple(getJetMETArgs(tup, isMC=True, forMET=True, isMETFixEE2017=True)), met_vars))
    yield res

@pytest.fixture(scope="module")
def nanojetargsMC18_postvalues_hem():
    import ROOT as gbl
    f = gbl.TFile.Open(os.path.join(testData, "DY_M50_2018postproc_hemfix.root"))
    tup = f.Get("Events")
    res = []
    for i in range(tup.GetEntries()):
        tup.GetEntry(i)
        jet_vars = {
            "nominal" : (toRVecFloat(tup.Jet_pt_nom    ), toRVecFloat(tup.Jet_mass_nom    )),
            }
        from itertools import chain
        jet_vars.update(dict(chain.from_iterable(
            { f"{src}{vdir.lower()}" : tuple(toRVecFloat(getattr(tup, f"Jet_{ivar}_{src}{vdir}".format(ivar, src))) for ivar in ("pt", "mass"))
                for vdir in ("Up", "Down") }.items() for src in ("jer", "jesAbsoluteStat", "jesAbsoluteScale", "jesHEMIssue")
            )))
        res.append((getJetMETArgs(tup, isMC=True, forMET=False, addHEM2018Issue=True), jet_vars))
    yield res

@pytest.fixture(scope="module")
def nanofatjetargsMC16_postvalues():
    import ROOT as gbl
    f = gbl.TFile.Open(os.path.join(testData, "DY_M50_2016postproc_FatJetKin.root"))
    tup = f.Get("Events")
    res = []
    for i in range(tup.GetEntries()):
        tup.GetEntry(i)
        jet_vars = {
            "nominal" : (toRVecFloat(tup.FatJet_pt_nom), toRVecFloat(tup.FatJet_mass_nom), toRVecFloat(tup.FatJet_msoftdrop_nom)),
            }
        from itertools import chain
        jet_vars.update(dict(chain.from_iterable(
            { f"jer{i:d}{vdir.lower()}" : tuple(toRVecFloat(getattr(tup, f"FatJet_{ivar}_jer{i:d}{vdir}")) for ivar in ("pt", "mass", "msoftdrop"))
                for vdir in ("Up", "Down") }.items() for i in range(6))))
        jet_vars.update(dict(chain.from_iterable(
            { f"jes{src}{vdir.lower()}" : tuple(toRVecFloat(getattr(tup, f"FatJet_{ivar}_jes{src}{vdir}".format(ivar, src))) for ivar in ("pt", "mass", "msoftdrop"))
                for vdir in ("Up", "Down") }.items() for src in ("AbsoluteStat", "AbsoluteScale")
            )))
        jet_vars.update(dict(chain.from_iterable(
            { f"{src}{vdir.lower()}" : tuple(toRVecFloat(getattr(tup, f"FatJet_{ivar}_{src}{vdir}".format(ivar, src))) for ivar in ("mass", "msoftdrop"))
                for vdir in ("Up", "Down") }.items() for src in ("jmr", "jms")
            )))
        res.append((getFatJetArgs(tup, isMC=True), jet_vars))
    yield res

@pytest.fixture(scope="module")
def nanofatjetargsMC16():
    import ROOT as gbl
    f = gbl.TFile.Open(os.path.join(testData, "DY_M50_2016postproc_FatJetKin.root"))
    for tup in getEventWith(f):
        yield getFatJetArgs(tup, isMC=True, addHEM2018Issue=False)

@pytest.fixture(scope="module")
def jetvarcalc_empty():
    import ROOT as gbl
    loadJMESystematicsCalculators()
    calc = gbl.JetVariationsCalculator()
    yield calc

@pytest.fixture(scope="module")
def jetvarcalcMC16_smear():
    import ROOT as gbl
    loadJMESystematicsCalculators()
    calc = gbl.JetVariationsCalculator()
    configureCalc(calc, jerTag="Summer16_25nsV1_MC", splitJER=True)
    yield calc

@pytest.fixture(scope="module")
def jetvarcalcMC16_jec():
    import ROOT as gbl
    loadJMESystematicsCalculators()
    calc = gbl.JetVariationsCalculator()
    configureCalc(calc, jerTag="Summer16_25nsV1_MC", splitJER=True,
            jecTag="Summer16_07Aug2017_V11_MC", levels=["L1FastJet", "L2Relative"])
    yield calc

@pytest.fixture(scope="module")
def jetvarcalcMC16_jesunc():
    import ROOT as gbl
    loadJMESystematicsCalculators()
    calc = gbl.JetVariationsCalculator()
    configureCalc(calc, jerTag="Summer16_25nsV1_MC", splitJER=True,
            jecTag="Summer16_07Aug2017_V11_MC", levels=["L1FastJet", "L2Relative"],
            uncSources=["AbsoluteStat", "AbsoluteScale"])
    yield calc

@pytest.fixture(scope="module")
def metvarcalcMC16_jesunc():
    import ROOT as gbl
    loadJMESystematicsCalculators()
    calc = gbl.Type1METVariationsCalculator()
    calc.setUnclusteredEnergyTreshold(15.)
    configureCalc(calc, jerTag="Summer16_25nsV1_MC", splitJER=True,
            jecTag = "Summer16_07Aug2017_V11_MC", levels=["L1FastJet", "L2Relative"],
            uncSources=["AbsoluteStat", "AbsoluteScale"], levels_l1=["L1FastJet"])
    yield calc

@pytest.fixture(scope="module")
def metvarcalcMC17_FixEE():
    ## a better test even would be to have a different JEC than production
    import ROOT as gbl
    loadJMESystematicsCalculators()
    calc = gbl.FixEE2017Type1METVariationsCalculator()
    calc.setUnclusteredEnergyTreshold(15.)
    configureCalc(calc, jerTag="Fall17_V3_MC", splitJER=False,
            jecTag="Fall17_17Nov2017_V32_MC", levels=["L1FastJet", "L2Relative"],
            uncSources=["AbsoluteStat", "AbsoluteScale"], levels_l1=["L1FastJet"])
    yield calc

@pytest.fixture(scope="module")
def jetvarcalcMC18_hem():
    import ROOT as gbl
    loadJMESystematicsCalculators()
    calc = gbl.JetVariationsCalculator()
    configureCalc(calc, jerTag="Autumn18_V7b_MC", splitJER=False,
            jecTag="Autumn18_V19_MC", levels=["L1FastJet", "L2Relative"],
            uncSources=["AbsoluteStat", "AbsoluteScale"])
    calc.setAddHEM2018Issue(True)
    yield calc

@pytest.fixture(scope="module")
def fatjetvarcalcMC16():
    import ROOT as gbl
    loadJMESystematicsCalculators()
    tp = gbl.FatJetVariationsCalculator
    calc = gbl.FatJetVariationsCalculator()
    configureCalc(calc, jetType="AK8PFPuppi", jerTag="Summer16_25nsV1_MC", splitJER=True,
            jecTag="Summer16_07Aug2017_V11_MC", levels=["L1FastJet", "L2Relative"],
            uncSources=["AbsoluteStat", "AbsoluteScale"])
    configureFatJetCalc(calc, isMC=True, year="2016", doSmearing=True, isTau21DDT=False)
    yield calc

def test_jetvarcalc_empty(jetvarcalc_empty):
    assert jetvarcalc_empty

def test_jetvarcalcMC16_smear(jetvarcalcMC16_smear):
    assert jetvarcalcMC16_smear

def test_jetvarcalcMC16_nano_smear(jetvarcalcMC16_smear, nanojetargsMC16):
    res = jetvarcalcMC16_smear.produce(*nanojetargsMC16)
    assert res

def test_jetvarcalcMC16_nano_jec(jetvarcalcMC16_jec, nanojetargsMC16):
    res = jetvarcalcMC16_jec.produce(*nanojetargsMC16)
    assert res

def test_jetvarcalcMC16_nano_jesunc(jetvarcalcMC16_jesunc, nanojetargsMC16):
    res = jetvarcalcMC16_jesunc.produce(*nanojetargsMC16)
    assert res

def test_fatjetvarcalcMC16_nano(fatjetvarcalcMC16, nanofatjetargsMC16):
    res = fatjetvarcalcMC16.produce(*nanofatjetargsMC16)
    assert res

import math
def isclose_float(a, b, tol=1.):
    import ROOT as gbl
    return math.isclose(a, b, rel_tol=tol*getattr(gbl, "std::numeric_limits<float>").epsilon())

def compareJets(names, calcRes, postValues, tol=1.):
    hasDiff = False
    for ky,(post_pt, post_mass) in postValues.items():
        idx = names.index(ky)
        print(ky, "pt", calcRes.pt(idx), post_pt)
        print(ky, "m ", calcRes.mass(idx), post_mass)
        if not ( all(isclose_float(a,b, tol=tol) for a,b in zip(post_pt, calcRes.pt(idx))) and all(isclose_float(a,b, tol=tol) for a,b in zip(post_mass, calcRes.mass(idx))) ):
            print(f"FAIL: Difference for {ky}")
            hasDiff = True
    return not hasDiff

def compareFatJets(namesAll, namesM, calcRes, postValues, tol=1.):
    hasDiff = False
    for ky,post in postValues.items():
        if ky in namesAll:
            idx = namesAll.index(ky)
            (post_pt, post_mass, post_msd) = post
            print(ky, "pt ", calcRes.pt(idx), post_pt)
            print(ky, "m  ", calcRes.mass(idx), post_mass)
            print(ky, "msd", calcRes.msoftdrop(idx), post_msd)
            eq_pt  = all(isclose_float(a,b, tol=tol) for a,b in zip(post_pt, calcRes.pt(idx)))
            eq_m   = all(isclose_float(a,b, tol=tol) for a,b in zip(post_mass, calcRes.mass(idx)))
            eq_msd = all(isclose_float(a,b, tol=tol) for a,b in zip(post_msd, calcRes.msoftdrop(idx)))
            if not ( eq_pt and eq_m and eq_msd ):
                what = ", ".join((["pt"] if not eq_pt else [])+(["mass"] if not eq_m else [])+(["msd"] if not eq_msd else []))
                print(f"FAIL: Difference for {ky} (in {what})")
                hasDiff = True
        else:
            idx = len(namesAll)+namesM.index(ky)
            (post_mass, post_msd) = post
            print(ky, "m  ", calcRes.mass(idx), post_mass)
            print(ky, "msd", calcRes.msoftdrop(idx), post_msd)
            eq_m   = all(isclose_float(a,b, tol=tol) for a,b in zip(post_mass, calcRes.mass(idx)))
            eq_msd = all(isclose_float(a,b, tol=tol) for a,b in zip(post_msd, calcRes.msoftdrop(idx)))
            if not ( eq_m and eq_msd ):
                what = ", ".join((["mass"] if not eq_m else [])+(["msd"] if not eq_msd else []))
                print(f"FAIL: Difference for {ky} (in {what})")
                hasDiff = True
    return not hasDiff

def compareMET(names, calcRes, postValues, reltol_pt=1.e-6, reltol_phi=1.e-6, abstol_phi=1.e-6):
    hasDiff = False
    for ky,(post_pt, post_phi) in postValues.items():
        idx = names.index(ky)
        print(ky, "pt ", calcRes.pt(idx), post_pt)
        print(ky, "phi", calcRes.phi(idx), post_phi)
        if not ( math.isclose(calcRes.pt(idx), post_pt, rel_tol=1.e-6) and math.isclose(calcRes.phi(idx), post_phi, rel_tol=1.e-6, abs_tol=1.e-6) ):
            print(f"FAIL: Difference for {ky}")
            hasDiff = True
    return not hasDiff

def test_jetvarcalc_nanopost_jesunc(jetvarcalcMC16_jesunc, nanojetargsMC16_postvalues):
    for nanojetargsMC16, postValues in nanojetargsMC16_postvalues:
        assert compareJets([ str(nm) for nm in jetvarcalcMC16_jesunc.available() ],
                jetvarcalcMC16_jesunc.produce(*nanojetargsMC16),
                postValues)

def test_metvarcalc_nanopost_jesunc(metvarcalcMC16_jesunc, nanoMETargsMC16_postvalues):
    metvarcalcMC16_jesunc.setIsT1SmearedMET(False)
    for nanoMETargsMC16, postValues in nanoMETargsMC16_postvalues:
        assert compareMET([ str(nm) for nm in metvarcalcMC16_jesunc.available() ],
                metvarcalcMC16_jesunc.produce(*nanoMETargsMC16),
                postValues)

def test_metvarcalc_nanopost_jesunc_T1Smear(metvarcalcMC16_jesunc, nanoSmearMETargsMC16_postvalues):
    metvarcalcMC16_jesunc.setIsT1SmearedMET(True)
    for nanoMETargsMC16, postValues in nanoSmearMETargsMC16_postvalues:
        assert compareMET([ str(nm) for nm in metvarcalcMC16_jesunc.available() ],
                metvarcalcMC16_jesunc.produce(*nanoMETargsMC16), postValues)

def test_metvarcalc_nanopost_jesunc_MCFixEE2017(metvarcalcMC17_FixEE, nanoMETFixEE2017argsMC17_postvalues):
    metvarcalcMC17_FixEE.setIsT1SmearedMET(False)
    for nanoMETargsMC17FixEE, postValues in nanoMETFixEE2017argsMC17_postvalues:
        assert compareMET([ str(nm) for nm in metvarcalcMC17_FixEE.available() ],
                metvarcalcMC17_FixEE.produce(*nanoMETargsMC17FixEE), postValues,
                reltol_pt=2.e-6, reltol_phi=2.e-6, abstol_phi=2.e-6)

def test_jetvarcalc_nanopost_jesunc_HEM(jetvarcalcMC18_hem, nanojetargsMC18_postvalues_hem):
    for nanojetargsMC18, postValues in nanojetargsMC18_postvalues_hem:
        assert compareJets([ str(nm) for nm in jetvarcalcMC18_hem.available() ],
                jetvarcalcMC18_hem.produce(*nanojetargsMC18),
                postValues, tol=2.)

def test_fatjetvarcalc_nanopost_jesunc(fatjetvarcalcMC16, nanofatjetargsMC16_postvalues):
    for fatjetargs, postValues in nanofatjetargsMC16_postvalues:
        avlAll = [ str(nm) for nm in fatjetvarcalcMC16.available() ]
        avlM = [ str(nm) for nm in fatjetvarcalcMC16.available("mass") if nm not in avlAll ]
        assert compareFatJets(avlAll, avlM,
                fatjetvarcalcMC16.produce(*fatjetargs),
                postValues)
