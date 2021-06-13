import logging
logger = logging.getLogger(__name__)

import ROOT as gbl

RVec_float = getattr(gbl, "ROOT::VecOps::RVec<float>")
RVec_int = getattr(gbl, "ROOT::VecOps::RVec<int>")

def toRVecFloat(values):
    """ Convert values to an RVec<float> """
    res = RVec_float(len(values), 0.)
    for i,val in enumerate(values):
        res[i] = val
    return res

def toRVecInt(values):
    """ Convert values to an RVec<int> """
    res = RVec_int(len(values), 0)
    for i,val in enumerate(values):
        res[i] = val
    return res

def getJetMETArgs(tup, isMC=True, forMET=False, isMETFixEE2017=False, addHEM2018Issue=False):
    """ Get the input values for the jet/met variations calculator from a tree (PyROOT-style) """
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
        args += [ toRVecFloat(getattr(tup, "CorrT1METJet_{0}".format(varNm)))
                for varNm in ("rawPt", "eta", "phi", "area", "muonSubtrFactor") ]
        args += [ toRVecFloat([]), toRVecFloat([]) ]
        if isMETFixEE2017:
            args += [ tup.MET_phi, tup.MET_pt, tup.METFixEE2017_phi, tup.METFixEE2017_pt ]
    return args

def getFatJetArgs(tup, isMC=True, addHEM2018Issue=False):
    """ Get the input values for the fatjet variations calculator from a tree (PyROOT-style) """
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

# W-tagging PUPPI softdrop JMS values: https://twiki.cern.ch/twiki/bin/view/CMS/JetWtagging
# 2016 values
fatjet_jmsValues = {
    "2016": [1.00 , 1.0094, 0.9906],  # nominal, up, down
    "2017": [0.982, 0.986 , 0.978 ],
    # Use 2017 values for 2018 until 2018 are released
    "2018": [0.982, 0.986 , 0.978 ],
    }
# jet mass resolution: https://twiki.cern.ch/twiki/bin/view/CMS/JetWtagging
#nominal, up, down
fatjet_jmrValues = {
    "2016": [1.0 , 1.2 , 0.8 ],
    "2017": [1.09, 1.14, 1.04],
    # Use 2017 values for 2018 until 2018 are released
    "2018": [1.09, 1.14, 1.04],
    }
# JMS&JMR SD corr in tau21DDT region: https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetWtagging#tau21DDT_0_43
fatjet_gmsValues_tau21DDT = {
    "2016": [ 1.014, 1.021, 1.007 ],
    "2017": [ 0.983, 0.99 , 0.976 ],
    "2018": [ 1.000, 1.010, 0.990 ] # tau21DDT < 0.43 WP
    }
fatjet_gmrValues_tau21DDT = {
    "2016": [1.086, 1.176, 0.996],
    "2017": [1.080, 1.161, 0.999],
    "2018": [1.124, 1.208, 1.040]
    }
# PUPPI scale and resolutions
fatjet_puppi_msd_params = {
    "gen": "1.0062610283313527+(-1.061605139842829*((x*0.07999000770091785)^-1.2045376937033758))",
    "reco_cen" : (1.0930197734452352, -0.00015006788774298745, 3.4486611503791434e-07, -2.681003093847412e-10, 8.674402290776817e-14, -1.0011358570698617e-17),
    "reco_fwd" : (1.2721151537214315, -0.0005716403627542301, 8.372894123074334e-07, -5.204332049858346e-10, 1.4537520981877012e-13, -1.5038869243803616e-17),
    "resol_cen": (1.092735080341856, 4.142622682579229e-05, -1.3736805733597026e-07, 1.2295818250625584e-10, -4.197075395161288e-14, 4.923792745086088e-18),
    "resol_fwd": (1.1649278339374347, -0.00012678902807057208, 1.0594037344842974e-07, 6.072087019460118e-12, -1.992427482862693e-14, 3.644006510237158e-18)
    }

def configureCalc(calc, jecTag=None, jerTag=None, jetType="AK4PFchs", levels=None, levels_l1=None, splitJER=False, uncSources=None, jecDBCache=None, jrDBCache=None):
    """ Configure a AK4 jet and MET variations calculator """
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

def configureFatJetCalc(calc, isMC=False, year=None, doSmearing=False, jmr=None, jms=None, isTau21DDT=False, gmr=None, gms=None, puppiGen=None, puppiRecoCorrCen=None, puppiRecoCorrFwd=None, puppiResolCen=None, puppiResolFwd=None, jecDBCache=None, jrDBCache=None):
    """ Configure a AK8 jet variations calculator """
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
        logger.debug("Setting JMS values {0}".format(jms))
        if hasattr(jms, "__iter__"):
            calc.setJMSValues(*jms)
        else: # just nominal
            calc.setJMSValues(jms)
    if gms is not None:
        logger.debug("Setting GMS values {0}".format(gms))
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
        logger.debug("Setting JMR values: {0}".format(jmr))
        calc.setJMRValues(*jmr)
        logger.debug("Setting GMR values: {0}".format(gmr))
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
    arrType = getattr(gbl, "std::array<double, 6>")
    def toArr(*args):
        assert len(args) == 6
        arr = arrType()
        for i,v in enumerate(args):
            arr[i] = v
        return arr
    logger.debug("Setting PUPPI gen formula to {0}".format(puppiGen))
    calc.setPuppiCorrections(puppiGen, *(toArr(*arg) for arg in (puppiRecoCorrCen, puppiRecoCorrFwd, puppiResolCen, puppiResolFwd)))
