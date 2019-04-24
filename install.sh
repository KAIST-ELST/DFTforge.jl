#!/bin/bash
echo 'This is JX installer. If Julia is already installed simply type julia install.jl'

echo ' To install julia see: https://julialang.org/downloads/'
echo ' Platform specfic detail see: https://julialang.org/downloads/platform.html'

case "$OSTYPE" in
  #solaris*) echo "SOLARIS" ;;
  darwin*)  echo "OSX"
  brew cask install julia
  ;;
  linux*)
  echo "Linux"
  mkdir -p ~/opt/bin
  echo 'export PATH=~/opt/bin:$PATH' >>~/.profile
  echo 'export PATH=~/opt/bin:$PATH' >>~/.bashrc
  JULIA_INSTALL=~/opt/bin  bash -ci "$(curl -fsSL https://raw.githubusercontent.com/abelsiqueira/jill/master/jill.sh)"
  ;;
  bsd*)     echo "BSD"
  echo 'Visit https://julialang.org/downloads/'
  ;;
  *)        echo "unknown: $OSTYPE"
  echo 'Visit https://julialang.org/downloads/'
  ;;
esac


echo 'julia install.jl'
# install DFTforge
julia install.jl
