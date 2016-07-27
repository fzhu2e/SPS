#!/bin/bash

mkdir Publish
cp sp_*.f90 Makefile README.markdown Publish
rm -r Publish/sp_*_bak.f90
mkdir Publish/debug
mkdir Publish/output
mkdir Publish/scripts
cp scripts/*.ncl Publish/scripts

exit 0
