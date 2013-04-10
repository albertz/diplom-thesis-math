#!/bin/bash

SageDir="/Applications/sage-5.4"
SageLocalIncludeDir="$SageDir/local/include"
SageIncludeDir="$SageDir/devel/sage/c_lib/include"
SageExtIncludeDir="$SageDir/devel/sage/sage/ext"
PythonIncludeDir="$SageLocalIncludeDir/python2.7"

cython \
	-I $PythonIncludeDir \
	-I $SageIncludeDir \
	-I $SageLocalIncludeDir \
	-I $SageExtIncludeDir \
	--cplus \
	algo_cython.pyx

c++ \
	-std=gnu++11 -stdlib=libc++ \
	-I $PythonIncludeDir \
	-I $SageIncludeDir \
	-I $SageLocalIncludeDir \
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
