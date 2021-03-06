#!/bin/bash

# Hermitian modular forms, https://github.com/albertz/diplom-thesis-math
# Copyright (c) 2013, Albert Zeyer, www.az2000.de
# This code is under the GPL v3 or later, see License.txt in the root directory of this project.

[ "$SageDir" = "" ] && {
	SageDir="---"
	[ "$SAGEDIR" != "" ] && SageDir=$SAGEDIR
}
[ \! -d $SageDir ] && SageDir="/Applications/sage-5.9"
[ \! -d $SageDir ] && SageDir=~/sage-5.9
[ \! -d $SageDir ] && SageDir=/usr/lib/sagemath
[ \! -d $SageDir ] && SageDir=/usr/local/share/sage
[ \! -d $SageDir ] && { echo "sagedir not found!"; exit 1; }

SageDevelDir="$SageDir/devel/sage"
SageLocalIncludeDir="$SageDir/local/include"
SageCLibDir="$SageDevelDir/c_lib"
SageCLibIncludeDir="$SageCLibDir/include"
SageExtIncludeDir="$SageDevelDir/sage/ext"
PythonIncludeDir="$SageLocalIncludeDir/python2.7"

set -x

cython \
	-I $PythonIncludeDir \
	-I $SageDevelDir \
	-I $SageLocalIncludeDir \
	-I $SageCLibIncludeDir \
	-I $SageExtIncludeDir \
	--line-directives \
	--cplus \
	algo_cython.pyx || exit 1

# Note that you either need a recent Clang or at least GCC 4.7.
Cpp="c++"
CppOpts="-ftrapv"
[ "$(uname)" == "Darwin" ] && CppOpts="$CppOpts -std=gnu++11 -stdlib=libc++"
[ "$(uname)" == "Linux" ] && CppOpts="$CppOpts -std=c++0x -fPIC -DOLDGCC"

$Cpp \
	$CppOpts \
	-I $PythonIncludeDir \
	-I $SageDevelDir \
	-I $SageLocalIncludeDir \
	-I $SageCLibIncludeDir \
	-I $SageExtIncludeDir \
	-ggdb \
	-c \
	algo_cython.cpp || exit 1

if [ "$(uname)" == "Linux" ]; then
$Cpp -shared algo_cython.o -ggdb -lc -L $SageCLibDir -lcsage -o algo_cython.so || exit 1

elif [ "$(uname)" == "Darwin" ]; then
libtool -dynamic \
	algo_cython.o \
	-undefined dynamic_lookup \
	-o algo_cython.so || exit 1

else
echo "no linker"
exit 1
fi

