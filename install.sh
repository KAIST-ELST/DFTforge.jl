#!/bin/bash

JULIA_INSTALL=~/opt/bin  bash -ci "$(curl -fsSL https://raw.githubusercontent.com/abelsiqueira/jill/master/jill.sh)"

# install DFTforge
julia install.jl
