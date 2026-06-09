#!/usr/bin/env python
from setuptools import setup, find_packages
import os

# Read the version from version.py without importing the package
version_file = os.path.join(os.path.dirname(__file__), 'src', 'enzywizard_dock', 'version.py')
with open(version_file) as f:
    exec(f.read())  # defines __version__

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="enzywizard-dock",
    version=__version__,                     # dynamically read from version.py (1.0.1)
    author="bioinfbrad",
    description=(
        "Perform molecular docking of one or multiple substrates with a cleaned protein structure "
        "and generate a detailed JSON report."
    ),
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/bioinfbrad/enzywizard-dock",
    package_dir={"": "src"},
    packages=find_packages(where="src"),
    python_requires=">=3.10",
    install_requires=[
        "rdkit>=2026.03.1",          # for cheminformatics, structure generation
        "numpy>=1.23.5",             # for numerical operations
        "biopython>=1.86",           # for protein structure handling
        "vina>=1.2.6",               # AutoDock Vina Python bindings
        "meeko>=0.7.1",              # ligand preparation for docking
        "bio-pyvol>=1.7.8",          # pocket detection (PyVOL)
    ],
    entry_points={
        "console_scripts": [
            "enzywizard-dock = enzywizard_dock.cli:main",
        ],
    },
    include_package_data=True,
    license="MIT",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3.10",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Chemistry",
    ],
)
