#!/bin/bash

SageDir="/Applications/sage-5.4"
SageLocalIncludeDir="$SageDir/local/include"
SageIncludeDir="$SageDir/devel/sage/c_lib/include"
SageExtIncludeDir="$SageDir/devel/sage/sage/ext"

cython \
	-I $SageIncludeDir \
	-I $SageLocalIncludeDir \
	-I $SageExtIncludeDir \
	--cplus \
	algo.pyx

