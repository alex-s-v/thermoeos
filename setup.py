# -*- coding: utf-8 -*-
"""Utilities for modeling vapor-liquid equilibrium of multicomponent mixtures.
Copyright (C) 2019 Alexandr Vasilyev <alexandr.s.vasilyev@gmail.com>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE."""

from setuptools import setup

setup(
    name="thermoeos",
    packages=["thermoeos"],
    license="MIT",
    version="0.1.0",
    description=[
        "Utilities for modeling vapor-liquid equilibrium"
        "of multicomponent mixtures."
    ],
    author="Alexandr Vasilyev",
    install_requires=["scipy", "numpy"],
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    platforms=["Windows", "Linux", "Mac OS", "Unix"],
    author_email="alexandr.s.vasilyev@gmail.com",
    url="https://github.com/alex-s-v/thermoeos",
    keywords=["chemical engineering", "chemistry",
              "thermodynamics", "cheminformatics", "engineering",
              "vapor pressure", "equation of state"],
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Developers",
        "Intended Audience :: Education",
        "Intended Audience :: Manufacturing",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Operating System :: MacOS",
        "Operating System :: Microsoft :: Windows",
        "Operating System :: POSIX",
        "Operating System :: POSIX :: BSD",
        "Operating System :: POSIX :: Linux",
        "Operating System :: Unix",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Topic :: Education",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Scientific/Engineering :: Physics",
    ],
    include_package_data=True
)
