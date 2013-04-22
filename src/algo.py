
from sage.matrix.constructor import matrix
from sage.modular.congroup import Gamma0
from sage.modular.modform.constructor import ModularForms
from sage.modules.free_module import FreeModule
from sage.rings.number_field.number_field import QQ, ZZ
from sage.rings.power_series_ring import PowerSeriesRing
import sys
import algo_cython as C

def reloadC():
	"""
	This is just for testing to reload the C (Cython) module
	after it was recompiled.
	Note that this code is very unsafe! It will likely crash
	when you have other references on the C module.
	This is only for debugging and development!
	"""
	global C
	import ctypes
	try:
		libdl = ctypes.CDLL("libdl.so")
	except Exception:
		# MacOSX:
		libdl = ctypes.CDLL("libdl.dylib")
	libdl.dlclose.argtypes = [ctypes.c_void_p]
	so = ctypes.PyDLL(C.__file__)
	assert(libdl.dlclose(so._handle) == 0)
	reload(C)

def test_algo_calcMatrix():
	calc = C.Calc()
	calc.init(D = -4, HermWeight = 10)

	calc.getNextS()
	calc.getNextS()

	calc.calcMatrix()
	return calc.getMatrix()

Verbose = True

def modform(D, HermWeight, B_cF=10):
	"Main algo"

	calc = C.Calc()
	calc.init(D = D, HermWeight = HermWeight, B_cF=B_cF)
	calc.calcReducedCurlF()
	reducedCurlFSize = calc.matrixColumnCount

	# Step 1. Iterate through square-free numbers l, starting at 1.
	# Init curlS = {}.
	curlS = []

	herm_modform_fe_expannsion = FreeModule(QQ, reducedCurlFSize)

	# Calculate the dimension of Hermitian modular form space.
	if D == -3:
		# dern2003graded, Thm 7
		R = PowerSeriesRing(ZZ, name="t", default_prec = HermWeight + 1)
		t = R.gen()
		dims = (1 + t**45) / (1 - t**4 ) / ( 1 - t**6 ) / ( 1 - t**9 ) / ( 1 - t**10 ) / ( 1 - t**12 )
		dim = dims[HermWeight]
	#elif D == -4:
		# dern2003graded, Corollary 9 and Lemma 3
		# TODO...
		#R = PowerSeriesRing(ZZ, name="t", default_prec = HermWeight + 1)
		#t = R.an_element()
	else:
		raise NotImplementedError, "dimension calculation of Hermitian modular form with D = %i not implemented" % D

	print "current dimension:", herm_modform_fe_expannsion.dimension(), "wanted:", dim

	while True:
		# Step 3. Iterate S \in Mat_2^T(Z). Add to curlS. iterate by denominator.
		# S_11 and S_22 (diagonal entries) are positive.
		# S positive definite.
		# S can be changed arbitrarily by GL(2, \ZZ).
		calc.curlS_clearMatrices()
		S = calc.getNextS()
		curlS += [S]
		if Verbose: print "trying S=", S, "det=", S.det()

		# Step 4. Calculate restriction matrix.
		if Verbose: sys.stdout.write("calc restriction matrix..."); sys.stdout.flush()
		calc.calcMatrix()
		M_S = calc.getMatrix()
		if Verbose: print "done:", M_S

		precLimit = M_S.nrows() # \cF(S)

		# These are the Elliptic modular forms with weight 2*HermWeight to \Gamma_0(l).
		l = S.det()
		mf = ModularForms(Gamma0(l), 2 * HermWeight)
		fe_expansion_matrix_l = matrix(QQ, [b.qexp(precLimit).padded_list(precLimit) for b in mf.basis()])
		fe_expansion_matrix_l.echelonize()

		# or:  fe_expansion_matrix[:n2,:].row_module()
		ell_modform_fe_expansions_l = fe_expansion_matrix_l.row_module()

		if Verbose: sys.stdout.write("calc M_S_module..."); sys.stdout.flush()
		M_S_module = M_S.column_module()
		if Verbose: print "done"
		restriction_fe_expansions = ell_modform_fe_expansions_l.intersection( M_S_module )
		herm_modform_fe_expannsion_S = M_S.solve_right( restriction_fe_expansions.basis_matrix().transpose() )
		herm_modform_fe_expannsion_S += M_S.right_kernel()
		# TODO: or row_space?
		herm_modform_fe_expannsion_S_module = herm_modform_fe_expannsion_S.column_module()

		herm_modform_fe_expannsion = herm_modform_fe_expannsion.intersection( herm_modform_fe_expannsion_S_module )
		current_dimension = herm_modform_fe_expannsion.dimension()
		print "current dimension:", current_dimension, "wanted:", dim

		# Step 5. dimension check
		if dim == current_dimension:
			break


	# TODO:
	# Otherwise, reconstruct fourier expansion.
	return herm_modform_fe_expannsion
