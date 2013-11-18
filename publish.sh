#!/bin/bash

mkdir Publish
cp sp_*.f90 Makefile README.markdown Publish
mkdir Publish/debug
mkdir Publish/output
mkdir Publish/scripts
cp scripts/*.ncl Publish/scripts

exit 0
