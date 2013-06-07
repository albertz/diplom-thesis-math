# -*- coding: utf-8 -*-
# Hermitian modular forms, https://github.com/albertz/diplom-thesis-math
# Copyright (c) 2013, Albert Zeyer, www.az2000.de
# This code is under the GPL v3 or later, see License.txt in the root directory of this project.

from sage.misc.prandom import randrange
from sage.rings.number_field.number_field import QQ
import algo_cython as C
from sage.matrix.constructor import matrix
from utils import *
from helpers import *


def fast_fail_test_D3_k6(B_cF=5):
	from checks import check_eisenstein_series_D3_weight6

	D = -3
	HermWeight = 6
	calc = C.Calc()
	calc.init(D = D, HermWeight = HermWeight, B_cF=B_cF)
	calc.calcReducedCurlF()

	while True:
		calc.curlS_clearMatrices()
		S = calc.getNextS()
		if S == matrix(2, 2, [3, 0, 0, 1]): break

	l = S.det()
	l = toInt(l)

	M_S = calcRestrictMatrix(calc) # matrix over integer ring
	M_S = M_S.matrix_over_field() # matrix over rational field

	precLimit = M_S.nrows()
	assert calcPrecisionDimension(B_cF=B_cF, S=S) == precLimit

	verbose("get elliptic modform space with precision %i ..." % precLimit)
	ell_dim, fe_expansion_matrix_l = getElliptModFormsBasisMatrix(l, 2*HermWeight, precLimit)
	assert ell_dim == fe_expansion_matrix_l.rank()
	ell_modform_fe_expansions_l = fe_expansion_matrix_l.row_module()
	verbose("dim of elliptic modform space: %i" % ell_modform_fe_expansions_l.dimension())

	verbose("calc M_S_module...")
	M_S_module = M_S.column_module()
	verbose("dimension of M_S column module: %i" % M_S_module.dimension())
	restriction_fe_expansions = ell_modform_fe_expansions_l.intersection( M_S_module )
	verbose("dimension of restriction_fe_expansions: %i" % restriction_fe_expansions.dimension())
	herm_modform_fe_expannsion_S = M_S.solve_right( restriction_fe_expansions.basis_matrix().transpose() )
	herm_modform_fe_expannsion_S_module = herm_modform_fe_expannsion_S.column_module()
	verbose("dimension of herm column module: %i" % herm_modform_fe_expannsion_S_module.dimension())
	verbose("calc M_S_right_kernel...")
	M_S_right_kernel = M_S.right_kernel()
	verbose("dimension of M_S right kernel: %i" % M_S_right_kernel.dimension())
	herm_modform_fe_expannsion_S_module += M_S_right_kernel

	try:
		check_eisenstein_series_D3_weight6(
			vs=herm_modform_fe_expannsion_S_module,
			B_cF=B_cF
		)
	except Exception:
		print "restriction_fe_expansions =", restriction_fe_expansions
		print "M_S_right_kernel =", M_S_right_kernel
		print "herm_modform_fe_expannsion_S_module =", herm_modform_fe_expannsion_S_module
		raise


def _fork_test_func(iterator=None):
	if not iterator:
		import itertools
		iterator = itertools.count()
	x = None
	for i in iterator:
		m = matrix(QQ, 100, [randrange(-100,100) for i in range(100*100)])
		x = m.kernel()
		print x
	return x

def fork_test():
	_fork_test_func(range(10))
	import os
	pid = os.fork()
	if pid != 0:
		print "parent, child: %i" % pid
		os.waitpid(pid, 0)
	else:
		print "child"
		try:
			_fork_test_func()
		finally:
			os._exit(0)

def fork_test2():
	# Also see here: http://www.sagemath.org/doc/reference/sage/parallel/decorate.html
	from sage.parallel.decorate import fork
	test_ = fork(_fork_test_func, verbose=True)
	test_()

def fork_test3(mustExec=False):
	_fork_test_func(range(10))
	import utils
	utils.asyncCall(func=_fork_test_func)


def parall_test(task_limit=1):
	import utils
	parallelizaton = utils.Parallelization(task_limit=task_limit)
	def task_iter_func():
		while True:
			yield lambda: _fork_test_func(range(10))
	parallelizaton.task_iter = task_iter_func()

	print parallelizaton.get_next_result()


def calcCurlSiter_serialization_test():
	D = -3
	HermWeight = 6
	B_cF = 7

	import algo_cython as C
	calc = C.Calc()
	calc.init(D = D, HermWeight = HermWeight, B_cF=B_cF)
	calc.calcReducedCurlF()

	for i in range(10):
		calc.getNextS()

	calc_state = pickle_dumps(calc)
	calc2 = pickle_loads(calc_state)

	for i in range(10):
		calc.getNextS()
		calc2.getNextS()

	assert calc.curlS == calc2.curlS


def test_calcPrecisionDimension(D=-3, HermWeight=6, B_cF=7):
	from helpers import calcPrecisionDimension, CurlO
	space = CurlO(D)

	import algo_cython as C
	calc = C.Calc()
	calc.init(D=D, HermWeight=HermWeight, B_cF=B_cF)
	calc.calcReducedCurlF()

	for i in range(1000):
		calc.curlS_clearMatrices()
		S = calc.getNextS()
		M = calc.calcMatrix()
		S_repr = (S[0,0], space.as_tuple_b(S[0,1]), S[1,1])
		assert M.nrows() == calcPrecisionDimension(B_cF=B_cF, S=S), "fail at i = %i, S = %r, M.nrows = %i" % (i, S_repr, M.nrows())


def test_herm_modform_indexset(D=-3, B_cF=7):
	# This also tests curlF iteration and reduceGL.
	from helpers import herm_modform_indexset_py
	from algo import herm_modform_indexset
	indexset = herm_modform_indexset(D=D, B_cF=B_cF)
	indexset_py = herm_modform_indexset_py(D=D, B_cF=B_cF)
	assert indexset == indexset_py


def test_calcMatrix(D=-3, HermWeight=6, B_cF=7):
	import algo_cython as C
	calc = C.Calc()
	calc.init(D=D, HermWeight=HermWeight, B_cF=B_cF)
	calc.calcReducedCurlF()

	from helpers import calcRestrictMatrix_py

	for i in range(1000):
		calc.curlS_clearMatrices()
		S = calc.getNextS()
		if i % 100 != 0: continue

		M_cpp = calc.calcMatrix()
		M_py = calcRestrictMatrix_py(D=D, HermWeight=HermWeight, B_cF=B_cF, S=S)

		assert M_cpp == M_py


def test_solveR():
	space = CurlO(-3)
	a,b,c,d = 2,1,1,1
	s,t,u = 5,space.Droot,1
	M = matrix(2, 2, [a,b,c,d])
	S = matrix(2, 2, [s,t,t.conjugate(),u])
	gamma,R,tM = solveR(M, S, space)

	a,b,c,d = 0,-1,1,0
	s,t,u = 1,0,2
	M = matrix(2, 2, [a,b,c,d])
	S = matrix(2, 2, [s,t,t.conjugate(),u])
	gamma,R,tM = solveR(M, S, space)

	a,b,c,d = 1,0,2,1
	s,t,u = 1,0,4
	M = matrix(2, 2, [a,b,c,d])
	S = matrix(2, 2, [s,t,t.conjugate(),u])
	gamma,R,tM = solveR(M, S, space)

	a,b,c,d = 1,0,3,1
	s,t,u = 1,2,16
	M = matrix(2, 2, [a,b,c,d])
	S = matrix(2, 2, [s,t,t.conjugate(),u])
	gamma,R,tM = solveR(M, S, space)

	a,b,c,d = 0,-1,1,0
	s,t,u = 2, QQ(0.5) * space.Droot - QQ(0.5), 2
	M = matrix(2, 2, [a,b,c,d])
	S = matrix(2, 2, [s,t,t.conjugate(),u])
	gamma,R,tM = solveR(M, S, space)

	a,b,c,d = 1,0,3,1
	s,t,u = 3, QQ(0.5) * space.Droot - QQ(1.5), 3
	M = matrix(2, 2, [a,b,c,d])
	S = matrix(2, 2, [s,t,t.conjugate(),u])
	gamma,R,tM = solveR(M, S, space)

	a,b,c,d = 1,0,2,1
	s,t,u = 3, QQ(0.5) * space.Droot - QQ(1.5), 3
	M = matrix(2, 2, [a,b,c,d])
	S = matrix(2, 2, [s,t,t.conjugate(),u])
	gamma,R,tM = solveR(M, S, space)

	a,b,c,d = 2,-1,3,-1
	s,t,u = 4, QQ(-0.5) * space.Droot + QQ(5/2.0), 4
	M = matrix(2, 2, [a,b,c,d])
	S = matrix(2, 2, [s,t,t.conjugate(),u])
	gamma,R,tM = solveR(M, S, space)

	a,b,c,d = 1,0,2,1
	s,t,u = 9, QQ(5)/2 * space.Droot - QQ(11)/2, 7
	M = matrix(2, 2, [a,b,c,d])
	S = matrix(2, 2, [s,t,t.conjugate(),u])
	gamma,R,tM = solveR(M, S, space)

	return gamma,R,tM


def test_calcMatrix_S3001(usePyImpl=False):
	D = -3
	HermWeight = 6
	B_cF = 7
	S = matrix(ZZ, 2, [3,0,0,1])
	l = S.det()
	if usePyImpl:	from helpers import calcRestrictMatrix_py as calcMatrix
	else:			from helpers import calcRestrictMatrix_any as calcMatrix
	M_S = calcMatrix(D=D, HermWeight=HermWeight, B_cF=B_cF, S=S)
	M_S = M_S.matrix_over_field() # matrix over rational field

	precLimit = M_S.nrows()
	assert precLimit == calcPrecisionDimension(B_cF=B_cF, S=S)

	# These are the Elliptic modular forms with weight 2*HermWeight to \Gamma_0(l).
	ell_dim, fe_expansion_matrix_l = getElliptModFormsBasisMatrix(l, 2*HermWeight, precLimit)
	assert fe_expansion_matrix_l.rank() == ell_dim
	ell_modform_fe_expansions_l = fe_expansion_matrix_l.row_module()
	assert ell_modform_fe_expansions_l.dimension() == ell_dim

	M_S_module = M_S.column_module()
	restriction_fe_expansions = ell_modform_fe_expansions_l.intersection( M_S_module )

	print "intersection:"
	print restriction_fe_expansions
	assert restriction_fe_expansions.dimension() > 0


def test_calcMatrix_S2a1(usePyImpl=False, B_cF=7):
	D = -3
	HermWeight = 6
	K = QuadraticField(D)
	a, b, c = 2, QQ(1)/2*K.gen() + QQ(1)/2, 1
	S = matrix(K, 2, [a, b, b.conjugate(), c])
	l = S.det()
	if usePyImpl:	from helpers import calcRestrictMatrix_py as calcMatrix
	else:			from helpers import calcRestrictMatrix_any as calcMatrix
	M_S = calcMatrix(D=D, HermWeight=HermWeight, B_cF=B_cF, S=S)
	M_S = M_S.matrix_over_field() # matrix over rational field

	precLimit = M_S.nrows()
	assert precLimit == calcPrecisionDimension(B_cF=B_cF, S=S)

	# These are the Elliptic modular forms with weight 2*HermWeight to \Gamma_0(l).
	ell_dim, fe_expansion_matrix_l = getElliptModFormsBasisMatrix(l, 2*HermWeight, precLimit)
	assert fe_expansion_matrix_l.rank() == ell_dim
	ell_modform_fe_expansions_l = fe_expansion_matrix_l.row_module()
	assert ell_modform_fe_expansions_l.dimension() == ell_dim

	M_S_module = M_S.column_module()
	restriction_fe_expansions = ell_modform_fe_expansions_l.intersection( M_S_module )
	assert restriction_fe_expansions.dimension() > 0

	herm_modform_fe_expannsion_S = M_S.solve_right( restriction_fe_expansions.basis_matrix().transpose() )
	herm_modform_fe_expannsion_S_module = herm_modform_fe_expannsion_S.column_module()
	M_S_right_kernel = M_S.right_kernel()
	herm_modform_fe_expannsion_S_module += M_S_right_kernel

	import checks
	checks.check_eisenstein_series_D3_weight6(
		vs=herm_modform_fe_expannsion_S_module,
		B_cF=B_cF
	)


