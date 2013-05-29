from src import algo_cython as C
from sage.matrix.constructor import matrix
from src.utils import _toInt, verbose
from src.helpers import calcRestrictMatrix, calcPrecisionDimension, getElliptModFormsBasisMatrix
from src.algo import _extra_check_on_herm_superspace


def _fast_fail_test_D3_k6(B_cF=5):
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
	l = _toInt(l)

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
		_extra_check_on_herm_superspace(
			vs=herm_modform_fe_expannsion_S_module,
			D=D, B_cF=B_cF, HermWeight=HermWeight
		)
	except Exception:
		print "restriction_fe_expansions =", restriction_fe_expansions
		print "M_S_right_kernel =", M_S_right_kernel
		print "herm_modform_fe_expannsion_S_module =", herm_modform_fe_expannsion_S_module
		raise