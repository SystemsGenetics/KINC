from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="kinc",
    version="4.0.0a",
    author="KINC Development Team",
    description="Knowledge Independent Network Construction in Python",
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
    ],
    python_requires=">=3.8",
    install_requires=[
        "numpy>=1.20.0",
        "scipy>=1.7.0",
        "pandas>=1.3.0",
        "h5py>=3.0.0",
        "scikit-learn>=1.0.0",
        "tqdm>=4.60.0",
    ],
    extras_require={
        "gpu": [
            "cupy-cuda12x>=12.0.0",
        ],
        "mpi": [
            "mpi4py>=3.1.0",
        ],
        "dask": [
            "dask[distributed]>=2021.10.0",
        ],
        "all": [
            "cupy-cuda12x>=12.0.0",
            "mpi4py>=3.1.0",
            "dask[distributed]>=2021.10.0",
        ],
        "dev": [
            "pytest>=7.0.0",
            "pytest-cov>=3.0.0",
            "black>=22.0.0",
            "flake8>=4.0.0",
        ],
    },
    entry_points={
        "console_scripts": [
            "kinc=kinc.cli:main",
        ],
    },
)
