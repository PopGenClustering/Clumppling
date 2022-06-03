import os
import sys
import shutil
import subprocess
from setuptools import setup, find_packages


if sys.version_info.major == 3:
    with open("README.md", "r", encoding="utf-8") as fh:
        long_description = fh.read()
else:
    # for python2
    with open("README.md", "r") as fh:
        long_description = fh.read()


setup(
    name='clumppling',
    version='0.0.1',
    author='Xiran Liu et al',
    zip_safe=False,
    packages=find_packages(),
    install_requires=["matplotlib","networkx","numpy","pandas","python_louvain>=0.16","scipy","cvxopt","cvxpy"],
    tests_require=[],
    include_package_data=True,
    long_description=long_description,
    long_description_content_type='text/markdown'
)   #yapf: disable