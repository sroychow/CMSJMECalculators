# CMSJMECalculators

This packages provides an efficient
[ROOT::RDataFrame](https://root.cern/doc/master/classROOT_1_1RDataFrame.html)-friendly
implementation of the recipes for jet and MET variations for the CMS experiment,
for use with samples in the NanoAOD format.
The code was adopted from the [bamboo](https://gitlab.cern.ch/cp3-cms/bamboo)
analysis framework.

## Installation

For using these helpers from python, the recommended solution is to install
the package (in a
[virtual](https://packaging.python.org/tutorials/installing-packages/#creating-virtual-environments)
or [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html)
environment) with
```python
pip install git+https://github.com/pieterdavid/CMSJMECalculators.git
```
[scikit-build](https://scikit-build.readthedocs.io/en/latest/) is used to
compile the C++ components against the available ROOT distribution.

From C++ there are two options: [CMake](https://cmake.org/) and (inside
CMSSW)
[scram](https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideScram).

A standalone CMake build can be done using the standard commands
(after cloning the repository):
```bash
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=<your-prefix> [other-options] <source-clone>
make
make install
```
Please note that this will only install the C++ components, not the python
helpers (yet).

TODO add python installation there

Building with scram inside CMSSW is also straightforward:
```bash
cd $CMSSW_BASE/src
git clone -o upstream https://github.com/pieterdavid/CMSJMECalculators.git UserCode/CMSJMECalculators
scram b
```

## Usage

When installed as a python package, the necessary components can be loaded
with:
```python
from CMSJMECalculators import loadJMESystematicsCalculators
loadJMESystematicsCalculators()
```
Note that this will load the shared library and headers in
[cling](https://root.cern/cling/), the [ROOT](https://root.cern/) interpreter,
so they can from then on also be used in JITted code, e.g. from
[RDataFrame](https://root.cern/doc/master/classROOT_1_1RDataFrame.html).

When installed inside a CMSSW environment, the import should be modified to
```python
from UserCode.CMSJMECalculators.CMSJMECalculators import loadJMESystematicsCalculators
```

The variations are calculated by the C++ classes ``JetVariationsCalculator`` and
``FatJetVariationsCalculator`` for the AK4 and AK8 jet JER and JES variations, and
``Type1METVariationsCalculator`` and ``FixEE2017Type1METVariationsCalculator``
for the Type-1 MET variations, using the standard procedure or with the special
recipe for 2017 (Type-1 smeared or standard MET is a configuration option).
To use these, an instance should be created (with the C++ interpreter, to make it
available from JITted code), and additional configuration passed by calling
setter methods, e.g. in PyROOT:
```python
import ROOT as gbl
calc = gbl.JetVariationsCalculator()
calc = getattr(gbl, "myJetVarCalc")
calc = gbl.JetVariationsCalculator()
# redo JEC, push_back corrector parameters for different levels
jecParams = getattr(gbl, "std::vector<JetCorrectorParameters>")()
jecParams.push_back(gbl.JetCorrectorParameters(textfilepath))
calc.setJEC(jecParams)
# calculate JES uncertainties (repeat for all sources)
jcp_unc = gbl.JetCorrectorParameters(textfilepath_UncertaintySources)
calc.addJESUncertainty("Total", jcp_unc)
# Smear jets, with JER uncertainty
calc.setSmearing(textfilepath_PtResolution, textfilepath_SF,
    splitJER,       # decorrelate for different regions
    True, 0.2, 3.)  # use hybrid recipe, matching parameters
```
The varied jet pt's and masses can be obtained by calling the ``produce`` method
with the per-event quantities, converted to
[`ROOT::VecOps::RVec`](https://root.cern/doc/master/classROOT_1_1VecOps_1_1RVec.html):
```python
from CMSJMECalculators.utils import toRVecFloat, toRVecInt
jetVars = calc.produce(toRVecFloat(tree.Jet_pt), toRVecFloat(tree.Jet_eta), ...)
```
since the full list of arguments can be long, and depends on a few parameters
(for data the MC branches are not there, and not needed, and MET needs a few
additional inputs), a helper function is provided, which can be used as follows:
```python
from CMSJMECalculators.utils import getJetMETArgs
jetVars = calc.produce(*getJetMETargs(tree, isMC=True, forMET=False))
```
This will return an object that contains all the variations, e.g.
`jetVars.pt(0)` will return the `RVec` with new nominal jet PTs.
The corresponding names of the variations, which depend on the configuration,
can be retrieved from the calculator by calling its `available()` method.

### From (JITted) RDataFrame

When constructing the RDataFrame graph from python, the calculator needs to be
constructed directly from the cling interpreter, such that it is available in
the global C++ namespace for JITted code:
```python
gbl.gROOT.ProcessLine("JetVariationsCalculator myJetVarCalc{};")
calc = getattr(gbl, "myJetVarCalc")
```
the second line retrieves a reference from PyROOT, such that the configuration
methods can be called as above.

Inside the RDataFrame graph the varied jet pt's and masses can be defined as
a new column:
```python
df.Define("ak4JetVars", "myJetVarcalc.produce(Jet_pt, Jet_eta, Jet_phi, ...)")
```
(the full set of arguments is not reproduced here, but can be found from the
`utils.getJetMETargs` method; since RDataFrame uses `RVec` internally
no conversion is needed).

### From C++

The PyROOT example above relies on the automatically generated bindings, so
the C++ equivalent is almost identical, and straigthforward to obtain.
When calling the `produce` method outside RDataFrame, most of the arguments
may need to be converted to `RVec`, which fortunately supports all common
kinds of array interfaces.

TODO expand C++ examples

### Caching the text files

Since the JEC and JER parameter text files need to be downloaded from the
corresponding repositories, which are quite big, a helper is provided that
downloads only the files that are used, and caches them locally.
It can be used like this (see the tests for more examples):
```python
from CMSJMECalculators.jetdatabasecache import JetDatabaseCache
jecDBCache = JetDatabaseCache("JECDatabase", repository="cms-jet/JECDatabase")
jrDBCache = JetDatabaseCache("JRDatabase", repository="cms-jet/JRDatabase")
# usage example, returns the local path
pl = jecDBCache.getPayload("Summer16_07Aug2017_V11_MC", "L1FastJet", "AK4PFchs")
```
The cache can also be checked and updated with the `checkCMSJMEDatabaseCaches`
script, which has an interactive mode (`-i` flag) that will start an IPython
shell after constructing the two database cache helpers.

## Testing and development

A set of [pytest](https://docs.pytest.org/en/6.2.x/)-based tests are included,
to make sure the implementation stays consistent with the POG-provided python
version in [nanoAOD-tools](https://github.com/cms-nanoAOD/nanoAOD-tools).
The tests compare the contents of the pt and mass branches for all variations.
They can be run with
```python
pytest tests
```
or, inside a CMSSW environment where python2 is the default
```python
python3 -m pytest tests
```

TODO make tests python2-compatible, expand, scripts for larger tests samples?
