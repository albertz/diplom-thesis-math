# Important: Keep these imports at the top so that `sage utils.py` works.
# See: http://ask.sagemath.org/question/2628/run-python-file-from-command-line-in-sage
import sys
import sage.all

import utils

if __name__ == "__main__":
	utils.ExecingProcess.checkExec()

	print "this Python file has no functionality (except to be run as a worker process)"
