#!/usr/bin/env python
"""
Setup script for lammps_pyace package
"""

import os
import sys
from pathlib import Path
from setuptools import setup, Extension, find_packages
from setuptools.command.build_ext import build_ext

# Force a fixed build directory
build_temp_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "build"))


class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=''):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)


class CMakeBuild(build_ext):

    def build_extension(self, ext):
        import subprocess

        extdir = os.path.abspath(self.build_lib)
        python_out = os.path.join(extdir, "lammps_pyace")
        os.makedirs(python_out, exist_ok=True)

        cmake_args = [
            f"-DCMAKE_LIBRARY_OUTPUT_DIRECTORY={python_out}",
            f"-DPYTHON_EXECUTABLE={sys.executable}",
            "-DCMAKE_BUILD_TYPE=Debug",  # Changed to Debug for lldb debugging
        ]

        build_args = []
        
        # Set CMAKE_BUILD_PARALLEL_LEVEL to control the parallel build level
        # across all generators.
        if "CMAKE_BUILD_PARALLEL_LEVEL" not in os.environ:
            # Build in parallel if not specified
            build_args += ['-j4']

        build_temp = os.path.join(self.build_temp, ext.name)
        if not os.path.exists(build_temp):
            os.makedirs(build_temp)

        # Configure
        subprocess.check_call(['cmake', ext.sourcedir] + cmake_args, cwd=build_temp)
        
        # Build
        subprocess.check_call(['cmake', '--build', '.'] + build_args, cwd=build_temp)


setup(
    name='lammps_pyace',
    version='0.3.0',
    author='alphataubio',
    description='Python bindings for ML-PACE (Atomic Cluster Expansion)',
    long_description=open('README.md').read() if os.path.exists('README.md') else '',
    long_description_content_type='text/markdown',
    packages=find_packages(),
    package_data={
        'lammps_pyace': ['*.pckl'],
    },
    ext_modules=[CMakeExtension('lammps_pyace.basis', sourcedir='.')],
    cmdclass={'build_ext': CMakeBuild},
    python_requires='>=3.8',
    install_requires=[
        'numpy>=1.19.0',
    ],
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Physics',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
    ],
)
