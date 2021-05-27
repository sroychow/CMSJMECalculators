def loadJMESystematicsCalculators():
    import pkg_resources
    incDir = pkg_resources.resource_filename("CMSJMECalculators", "include")
    libDir = pkg_resources.resource_filename("CMSJMECalculators", "lib")
    libName = "libCMSJMECalculators"
    import os.path
    import ROOT as gbl
    st = gbl.gSystem.Load(os.path.join(libDir, libName))
    if st == -1:
        raise RuntimeError("Library {0} could not be found".format(libName))
    elif st == -2:
        raise RuntimeError("Version match for library {0}".format(libName))
    gbl.gInterpreter.AddIncludePath(incDir)
    gbl.gROOT.ProcessLine('#include "JMESystematicsCalculators.h"')
    getattr(gbl, "JetVariationsCalculator::result_t")  # trigger dictionary generation

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
