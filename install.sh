#!/bin/bash

################################################################################
## Hongkee Yoon Hongkeeyoon@kaist.ac.kr
## 2019.05
## https://kaist-elst.github.io/DFTforge.jl/
################################################################################


echo 'This is JX installer. If Julia is already installed simply type julia install.jl'

echo ' To install julia see: https://julialang.org/downloads/'
echo ' Platform specfic detail see: https://julialang.org/downloads/platform.html'

echo ' Please install hdf5 lib before the install (e.g. in Ubuntu `apt-get install hdf5-tools`)'
case "$OSTYPE" in
  #solaris*) echo "SOLARIS" ;;
  darwin*)  echo "OSX"
  brew cask install julia
  ;;
  linux*)
  echo "Linux"
  mkdir -p ~/bin
  echo 'export PATH=~/bin:$PATH' >>~/.profile
  echo 'export PATH=~/bin:$PATH' >>~/.bashrc
  JULIA_INSTALL=~/bin  bash -ci "$(curl -fsSL https://raw.githubusercontent.com/abelsiqueira/jill/master/jill.sh)"
  ;;
  bsd*)     echo "BSD"
  echo 'Visit https://julialang.org/downloads/'
  ;;
  *)        echo "unknown: $OSTYPE"
  echo 'Visit https://julialang.org/downloads/'
  ;;
esac

echo 're-open the shell'
echo 'julia install.jl'
# install DFTforge
# julia install.jl
