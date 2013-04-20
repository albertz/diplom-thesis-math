from sage.modules.free_module import FreeModule
from sage.rings.number_field.number_field import QQ
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

def modform(D, HermWeight):
	"Main algo"

	calc = C.Calc()
	calc.init(D = D, HermWeight = HermWeight)

	# Step 1. Iterate through square-free numbers l, starting at 1.
	# Init curlS = {}.

	# Step 3. Iterate S \in Mat_2^T(Z). Add to curlS. iterate by denominator.
	# S_11 and S_22 (diagonal entries) are positive.
	# S positive definite.
	# S can be changed arbitrarily by GL(2, \ZZ).

	while True:
		S = calc.getNextS()
		if Verbose: print "trying S=", S, "det=", S.det()

		# Step 3a. Choose B>0 as limit for precision curlF.
		# TODO: how? dependent on S?
		# ... calc.setPrecisionLimit(B)

		# Step 4. Calculate restriction matrix.
		calc.calcMatrix()
		M_S = calc.getMatrix()



		reducedCurlFSize = M_S.column_size()

		# http://sage.math.washington.edu/tmp/sage-2.8.12.alpha0/doc/ref/module-sage.matrix.matrix2.html
		#ell_modform_fe_expansions = TODO ...
		restriction_fe_expansions = ell_modform_fe_expansions.intersection( M_S.column_module() )
		herm_modform_fe_expannsion_S = M_S.solve_right( restriction_fe_expansions.basis_matrix().transpose() )


		herm_modform_fe_expannsion = FreeModule(QQ, reducedCurlFSize)

		# Step 5. dimension check

		# TODO
		#   This comes last. It's written in Dern.

		# If check fails, goto step 3 and enlarge curlS.

		# Otherwise, reconstruct fourier expansion.
