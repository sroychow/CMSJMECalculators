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
through with:
```python
from CMSJMECalculators.utils import loadJMESystematicsCalculators
loadJMESystematicsCalculators()
```
Note that this will load the shared library and headers in
[cling](https://root.cern/cling/), the [ROOT](https://root.cern/) interpreter,
so they can from then on also be used in JITted code, e.g. from
[RDataFrame](https://root.cern/doc/master/classROOT_1_1RDataFrame.html).

When installed inside a CMSSW environment, the import should be modified to
```python
from UserCode.CMSJMECalculators.CMSJMECalculators.utils import loadJMESystematicsCalculators
```

TODO more detail on how to call (improve helpers from tests), and a C++ example

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
