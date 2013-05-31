#!sage

# See here:
# http://ask.sagemath.org/question/2628/run-python-file-from-command-line-in-sage

import sys
from sage.all import *
from sage.calculus.predefined import x
from sage.misc.preparser import preparse

from sage.matrix.matrix2 import Matrix

print "Hello world"
