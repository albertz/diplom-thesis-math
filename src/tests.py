from sage.misc.prandom import randrange
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
	level = 17
	weight = 12
	n = 50
	for i in iterator:
		mf = ModularForms(Gamma0(level), weight)
		fe_expansion_matrix_l = matrix(QQ, [b.qexp(n).padded_list(n) for b in mf.basis()])
		fe_expansion_matrix_l.echelonize()
		fe_expansion_matrix_l.right_kernel()

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
