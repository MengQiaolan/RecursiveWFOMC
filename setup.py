import os
from setuptools import find_packages, setup

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name='wfomc_fo2',
    version='0.1',
    license='MIT',
    packages=find_packages(include=['wfomc_fo2', 'wfomc_fo2.*']),
    install_requires=["symengine",
                      "sympy",
                      "lark",
                      "numpy",
                      "networkx",
                      "contexttimer",
                      "logzero",
                      "pandas",
                      "pysat",
                      "tqdm",
                      "pynauty",
                      "dataclasses",
                      "PrettyPrintTree"],
    python_requires='>=3.8',
    long_description=read('README.md'),
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Programming Language :: Python :: 3.8",
        "Topic :: Scientific/Engineering :: Artificial Intelligence",
        "License :: OSI Approved :: MIT License",
    ],
)
