#!/bin/bash

SageDir="/Applications/sage-5.4"
SageDevelDir="$SageDir/devel/sage"
SageLocalIncludeDir="$SageDir/local/include"
SageCLibIncludeDir="$SageDevelDir/c_lib/include"
SageExtIncludeDir="$SageDevelDir/sage/ext"
PythonIncludeDir="$SageLocalIncludeDir/python2.7"

cython \
	-I $PythonIncludeDir \
	-I $SageDevelDir \
	-I $SageLocalIncludeDir \
	-I $SageCLibIncludeDir \
	-I $SageExtIncludeDir \
	--cplus \
	algo_cython.pyx

c++ \
	-std=gnu++11 -stdlib=libc++ \
	-ftrapv \
	-I $PythonIncludeDir \
	-I $SageDevelDir \
	-I $SageLocalIncludeDir \
	-I $SageCLibIncludeDir \
	-I $SageExtIncludeDir \
	-c \
	algo_cython.cpp

# on linux:
#["ld"] +
#["-L/usr/local/lib"] +
#infiles +
#options +
#["-lc"] +
#["-shared", "-o", outfile]

# mac:
libtool -dynamic \
	algo_cython.o \
	-undefined dynamic_lookup \
	-o algo_cython.so
