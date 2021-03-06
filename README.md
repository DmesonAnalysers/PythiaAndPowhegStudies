[![](https://img.shields.io/github/license/DmesonAnalysers/DmesonAnalysis?color=blue)](https://github.com/DmesonAnalysers/DmesonAnalysis/blob/master/LICENSE)
![](https://img.shields.io/github/languages/count/DmesonAnalysers/DmesonAnalysis?color=green)
![](https://img.shields.io/github/last-commit/DmesonAnalysers/DmesonAnalysis?color=red)

# PythiaAndPowhegStudies
Repository with scripts to perform studies with Pythia8 and Powheg

# Instructions for Pythia8 installation as ROOT plugin
- Download and installation of Pythia8

  Download and untar a specific Pythia8 release from [http://home.thep.lu.se/~torbjorn/Pythia.html](http://home.thep.lu.se/~torbjorn/Pythia.html). For example for version 8.243:
  
  ```bash
  wget http://home.thep.lu.se/~torbjorn/pythia8/pythia8243.tgz
  tar xvfz pythia8243.tgz
  cd pythia8243
  ```
  
  Configuration for 64bit systems:
  
  ```bash
  ./configure --enable-64bit --enable-shared
  ```
  
  All configuration options can be listed with
  
  ```bash
  ./configure --help
  ```
  
  They include also a variety of external packages that can be enabled via
  
  ```bash
  --with-PACKAGE[=DIR]
  ```
  
  They include ```root``` for the usage of ROOT trees and histograms with pythia, ```evtgen``` for the particle decays with the EvtGen package, ```powheg``` for the hard process production with POWHEGBOX matrix element executables, ```python``` for an interface to use pythia in python. 
  These options are however not necessary for the installation of Pythia8 as ROOT plugin.
  
  Compile Pythia8:
  ```bash
  make -J N
  ```
  
  On MacOs, copy the dynamic shared library:
  ```bash
  cd lib
  cp libpythia8.dylib libpythia8.so
  ```
  
  Setup the environment including in the ```.bashrc``` (```.bash_profile```)
  
  ```bash
  export PYTHIA8=$HOME/pythia8243
  export PYTHIA8DATA=$HOME/pythia8243/share/Pythia8/xmldoc
  DYLD_LIBRARY_PATH=$PYTHIA8/lib:$DYLD_LIBRARY_PATH
  ```
  
- Compilation of ROOT

  ROOT has to be recompiled with the following configurations:
  
  ```bash
  ./configure --all --with-gsl-incdir="/usr/local/include/gsl" --with-gsl-libdir="/usr/local/lib" --enable-pythia8 --with-pythia8-incdir=$PYTHIA8/include --with-pythia8-libdir=$PYTHIA8/lib
  ```
  
  If ROOT in the [AliPhysics](https://github.com/alisw/AliPhysics) installation is used, the ```root.sh``` recipy in ```alidist``` has to be modified, adding in the configurations:
  ```bash
  -Dpythia8=ON                                                                     \
  -DPYTHIA8_DIR=$PYTHIA8                                                           \
  -DPYTHIA8_INCLUDE_DIR=$PYTHIA8/include                                           \
  -DPYTHIA8_LIBRARY=$PYTHIA8/lib/libpythia8.so                                     \
  ```
  
  ROOT has then to be recompiled
  
- Test that installation was successful
  
  Open root session and execute Pythia8 tutorial:

  ```bash
  root -l
  .x $ROOTSYS/tutorials/pythia/pithia8.C
  ```

# Instructions for POWHEG-BOX installation
- Download and install LHAPDF

  For the installation follow the instructions in [https://lhapdf.hepforge.org/install.html](https://lhapdf.hepforge.org/install.html)
  
  Setup the environment including in the ```.bashrc``` (```.bash_profile```)
  
  ```bash
  export LHAPDF=/path/for/installation
  export PATH=$PATH:$LHAPDF/bin
  LD_LIBRARY_PATH=$LHAPDF/lib:$LD_LIBRARY_PATH
  export LHAPDF_DATA_PATH=$LHAPDF/share/LHAPDF
  ```
  
  Copy sets of pdfs in ```$LHAPDF_DATA_PATH``` from CERN CVMFS ```/cvmfs/sft.cern.ch/lcg/external/lhapdfsets/current/```
  
- Download and install POWHEG-BOX
  
  Download the package with
  ```bash
  svn checkout --username anonymous --password anonymoussvn://powhegbox.mib.infn.it/trunk/POWHEG-BOX  
  ```
  
  Install the desired process (heavy quarks ```hvq``` in this case)
  ```bash
  cd POHEG-BOX/hvq
  make pwhg_main
  make main-PYTHIA-lhef
  ```
