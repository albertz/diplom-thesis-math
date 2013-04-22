#!/bin/bash

SageDir="/Applications/sage-5.4"
[ \! -d $SageDir ] && SageDir=~/sage-5.8

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

CppOpts="-std=gnu++11 -ftrapv"
[ "$(uname)" == "Darwin" ] && CppOpts="$CppOpts -stdlib=libc++"]
[ "$(uname)" == "Linux" ] && CppOpts="$CppOpts -fPIC"

c++ \
	$CppOpts \
	-I $PythonIncludeDir \
	-I $SageDevelDir \
	-I $SageLocalIncludeDir \
	-I $SageCLibIncludeDir \
	-I $SageExtIncludeDir \
	-c \
	algo_cython.cpp

if [ "$(uname)" == "Linux" ]; then
ld  -dynamic-linker algo_cython.o -L/usr/local/lib -lc -shared -o algo_cython.so

elif [ "$(uname)" == "Darwin" ]; then
libtool -dynamic \
	algo_cython.o \
	-undefined dynamic_lookup \
	-o algo_cython.so

else
echo "no linker"
exit 1
fi

