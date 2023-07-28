from setuptools import setup, find_packages


requirements = []
with open('draw/requirements.txt', 'r') as fh:
    for line in fh:
        requirements.append(line.strip())

# PYTHON_VERSION = (3,10)

# if sys.version_info < PYTHON_VERSION:
#    sys.exit(f"Python >= {PYTHON_VERSION[0]}.{PYTHON_VERSION[1]} required.")

readme = open('README.md').read()

setup(
    name="draw",
    version='0.0.0',
    description="Draw is a plotting library for DREEM and SEISMIC",
    long_description=readme,
    author="Silvi Rouskin Lab",
    author_email="silvi@hms.harvard.edu",
    url="https://github.com/rouskinlab/draw",
    packages=find_packages(),
    package_dir={'draw': 'draw'},
    include_package_data=True,
    install_requires=requirements,
    python_requires=">=3.10",
)
