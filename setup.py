from setuptools import setup, find_packages


requirements = []
with open('requirements.txt', 'r') as fh:
    for line in fh:
        requirements.append(line.strip())

# PYTHON_VERSION = (3,10)

# if sys.version_info < PYTHON_VERSION:
#    sys.exit(f"Python >= {PYTHON_VERSION[0]}.{PYTHON_VERSION[1]} required.")

readme = open('README.md').read()

setup(
    name="seismograph",
    version='0.0.0',
    description="seismograph is a plotting library for SEISMIC",
    long_description=readme,
    author="Silvi Rouskin Lab",
    author_email="silvi@hms.harvard.edu",
    url="https://github.com/rouskinlab/seismograph",
    packages=find_packages(),
    package_dir={'seismograph': 'seismograph'},
    include_package_data=True,
    install_requires=requirements,
    python_requires=">=3.10",
)
