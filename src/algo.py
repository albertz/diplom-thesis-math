
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

def modform(D, HermWeight):
	"Main algo"
	
	calc = C.Calc()
	calc.init(D = D, HermWeight = HermHeight)

	# Step 1. Iterate through square-free numbers l, starting at 1.
	# Init curlS = {}.

	# Step 3. Iterate S \in Mat_2^T(Z). Add to curlS. iterate by denominator.
	# S_11 and S_22 (diagonal entries) are positive.
	# S positive definite.
	# S can be changed arbitrarily by GL(2, \ZZ).

	while True:
		S = calc.getNextS()

		# Step 3a. Choose B>0 as limit for precision curlF.
		# TODO: how? dependent on S?
		# ... calc.setPrecisionLimit(B)

		# Step 4. Calculate restriction matrix.
		M = calc.calcMatrix()

		# Step 5. dimension check
		
		# TODO
		#   This comes last. It's written in Dern.

		# If check fails, goto step 3 and enlarge curlS.

		# Otherwise, reconstruct fourier expansion.
