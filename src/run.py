# -*- coding: utf-8 -*-
#!sage

# Hermitian modular forms, https://github.com/albertz/diplom-thesis-math
# Copyright (c) 2013, Albert Zeyer, www.az2000.de
# All rights reserved.
# This code is under the GPL v3 or later, see License.txt in the root directory of this project.

# Important: Keep these imports at the top so that `sage run.py` works.
# See: http://ask.sagemath.org/question/2628/run-python-file-from-command-line-in-sage
import sys
import sage.all


def calc_herm_modforms():
	from sage.parallel.ncpus import ncpus
	kwargs = {
		"D": -3,
		"HermWeight": 6,
		"B_cF": 7,
		"task_limit": ncpus(),
	}
	print "Calculating Hermitian modular forms with %r" % kwargs
	import algo
	modforms = algo.herm_modform_space__parallel(**kwargs)
	assert modforms
	import utils
	utils.save(
		modforms,
		"herm_modforms__D%i_k%i_B%i__%i.sobj" %
		(kwargs["D"], kwargs["HermWeight"], kwargs["B_cF"], modforms.degree())
	)
	print modforms
	sys.exit(0)


if __name__ == "__main__":
	# Check if we are a worker process. In that case, this would not return.
	import utils
	utils.ExecingProcess.checkExec()

	# Normal start up.
	# Just do some calculation.
	calc_herm_modforms()

