from setuptools import find_packages
from setuptools_scm import get_version
from skbuild import setup

setup(
    packages=find_packages(where="python"),
    package_dir={"": "python"},
    cmake_install_dir="python",
    cmake_args=[f"-DCMSJMECALCULATOR_VERSION:STRING={get_version()}"],
)
