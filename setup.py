from setuptools import setup, find_packages
import sys, os

requirements = []
with open('requirements.txt', 'r') as fh:
    for line in fh:
        requirements.append(line.strip())

# PYTHON_VERSION = (3,10)

# if sys.version_info < PYTHON_VERSION:
#    sys.exit(f"Python >= {PYTHON_VERSION[0]}.{PYTHON_VERSION[1]} required.")

readme = open('README.md').read()

setup(
    name="seismic-graph",
    version='0.0.1',
    description="seismic-graph is a plotting library for SEISMIC",
    long_description=readme,
    author="Silvi Rouskin Lab",
    author_email="silvi@hms.harvard.edu",
    url="https://github.com/rouskinlab/seismic-graph",
    packages=find_packages(),
    package_dir={'seismic-graph': 'seismic-graph'},
    include_package_data=True,
    package_data={'seismic_graph': ['resources/*.feather']},
    install_requires=requirements,
    python_requires=">=3.10",
)
