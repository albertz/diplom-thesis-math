# -*- coding: utf-8 -*-
# Hermitian modular forms, https://github.com/albertz/diplom-thesis-math
# Copyright (c) 2013, Albert Zeyer, www.az2000.de
# This code is under the GPL v3 or later, see License.txt in the root directory of this project.

# In this file, the main algorithm is implemented which calculates the
# vector space of Fourier expansions of Hermitian modular forms.
# This is the function `herm_modform_space()`.

# At the time of writing, we need Sage 5.9 (5.8 might work) and these Sage patches:
# http://trac.sagemath.org/sage_trac/ticket/14240
# http://trac.sagemath.org/sage_trac/ticket/14497
# Download the patch files, go to $SAGEDIR/devel/sage, do:
#   hg import <patch_i>
#   sage -b

from sage.modular.arithgroup.congroup_sl2z import SL2Z
from sage.modular.congroup import Gamma0
from sage.modules.free_module import FreeModule
from sage.rings.infinity import Infinity
from sage.rings.arith import lcm
from sage.rings.number_field.number_field import QQ, ZZ, CyclotomicField
from sage.rings.power_series_ring import PowerSeriesRing
from sage.matrix.matrix2 import Matrix
from sage.misc.cachefunc import cached_function as sage_cached_function
from sage.rings.integer import Integer
from sage.modules.free_module_element import vector
from sage.structure.sage_object import SageObject
from sage.symbolic.ring import SymbolicRing
from sage.symbolic.expression import Expression
import algo_cython as C
from utils import *
from helpers import *
from checks import *
from pprint import pprint


@sage_cached_function
def herm_modform_indexset(D, B_cF):
	"""
	This is the precision of the index set for the Fourier coefficients of the
	Hermitian modular forms (calculated by `herm_modform_space()`).

	The index set \Lambda is the set of positive definite 2*2 matrices over \curlO^#.
	\curlF \subset \Lambda is a precision such that for [a,b,c] \in \curlF,
	we have 0 \le a,c \le B_cF.

	This function returns all the reduced matrices of \curlF. The reduction
	is given by `reduce_GL()` in `reduceGL.hpp`.

	The order of \curlF is predefined and fixed. It is such that
	[a,b,c] < [s,t,u]  <==>  max(a,c) < max(s,u).
	Thus, for B1 < B2, we have

		L1 = herm_modform_indexset(D, B1)
		L2 = herm_modform_indexset(D, B2)
		assert len(L1) < len(L2)
		assert L1 == L2[:len(L1)]
	"""

	HermWeight = 6 # Just some dummy value. `calc.init()` wants something valid.
	calc = C.Calc()
	calc.init(D=D, HermWeight=HermWeight, B_cF=B_cF)
	calc.calcReducedCurlF()
	return calc.getReducedCurlF()


def herm_modform_space_dim(D, HermWeight):
	"""
	Calculates and returns the dimension of the vector space
	of Hermitian modular forms of weight `HermWeight` over \Gamma,
	where \Gamma = \Sp_2(\curlO) and \curlO is the maximal order
	of \QQ(\sqrt{D}).
	"""

	if D == -3:
		# dern2003graded, Thm 7
		R = PowerSeriesRing(ZZ, name="t", default_prec = HermWeight + 1)
		t = R.gen()
		dims = (1 + t**45) / (1 - t**4 ) / ( 1 - t**6 ) / ( 1 - t**9 ) / ( 1 - t**10 ) / ( 1 - t**12 )
		return dims[HermWeight]
	#elif D == -4:
		# dern2003graded, Corollary 9 and Lemma 3
		# TODO...
		#R = PowerSeriesRing(ZZ, name="t", default_prec = HermWeight + 1)
		#t = R.an_element()
	else:
		raise NotImplementedError, "dimension calculation of Hermitian modular form with D = %i not implemented" % D



@persistent_cache(name="modform_cusp_info")
def modform_cusp_info(calc, S, l, precLimit):
	"""
	This goes through all the cusps and compares the space given by `(f|R)[S]`
	with the space of Elliptic modular forms expansion at those cusps.
	"""

	assert l == S.det()
	assert list(calc.curlS) == [S]

	D = calc.D
	HermWeight = calc.HermWeight
	reducedCurlFSize = calc.matrixColumnCount
	herm_modform_fe_expannsion = FreeModule(QQ, reducedCurlFSize)

	if not Integer(l).is_squarefree():
		# The calculation of the cusp expansion space takes very long here, thus
		# we skip them for now.
		return None

	for cusp in Gamma0(l).cusps():
		if cusp == Infinity: continue
		M = cusp_matrix(cusp)

		try:
			gamma, R, tM = solveR(M, S, space=CurlO(D))
		except Exception:
			print (M, S)
			raise
		R.set_immutable() # for caching, we need it hashable

		herm_modforms = herm_modform_fe_expannsion.echelonized_basis_matrix().transpose()
		ell_R_denom, ell_R_order, M_R = calcMatrixTrans(calc, R)
		CycloDegree_R = CyclotomicField(ell_R_order).degree()
		print "M_R[0] nrows, ell_R_denom, ell_R_order, Cyclo degree:", \
			M_R[0].nrows(), ell_R_denom, ell_R_order, CycloDegree_R

		# The maximum precision we can use is M_R[0].nrows().
		# However, that can be quite huge (e.g. 600).
		ce_prec = min(precLimit, M_R[0].nrows())

		ce = cuspExpansions(level=l, weight=2*HermWeight, prec=ce_prec)
		ell_M_denom, ell_M = ce.expansion_at(SL2Z(M))
		print "ell_M_denom, ell_M nrows:", ell_M_denom, ell_M.nrows()
		ell_M_order = ell_R_order # not sure here. just try the one from R. toCyclPowerBase would fail if this doesn't work
		# CyclotomicField(l / prod(l.prime_divisors())) should also work.

		# Transform to same denom.
		denom_lcm = int(lcm(ell_R_denom, ell_M_denom))
		ell_M = addRows(ell_M, denom_lcm / ell_M_denom)
		M_R = [addRows(M_R_i, denom_lcm / ell_R_denom) for M_R_i in M_R]
		ell_R_denom = ell_M_denom = denom_lcm
		print "new denom:", denom_lcm
		assert ell_R_denom == ell_M_denom

		# ell_M rows are the elliptic FE. M_R[i] columns are the elliptic FE.
		# We expect that M_R gives a higher precision for the ell FE. I'm not sure
		# if this is always true but we expect it here (maybe not needed, though).
		print "precision of M_R[0], ell_M, wanted:", M_R[0].nrows(), ell_M.ncols(), ce_prec
		assert ell_M.ncols() >= ce_prec
		prec = min(M_R[0].nrows(), ell_M.ncols())
		# cut to have same precision
		M_R = [M_R_i[:prec,:] for M_R_i in M_R]
		ell_M = ell_M[:,:prec]
		assert ell_M.ncols() == M_R[0].nrows() == prec

		print "M_R[0] rank, herm rank, mult rank:", \
			M_R[0].rank(), herm_modforms.rank(), (M_R[0] * herm_modforms).rank()
		ell_R = [M_R_i * herm_modforms for M_R_i in M_R]

		# I'm not sure on this. Seems to be true and it simplifies things in the following.
		assert ell_M_order <= ell_R_order, "{0}".format((ell_M_order, ell_R_order))
		assert ell_R_order % ell_M_order == 0, "{0}".format((ell_M_order, ell_R_order))

		# Transform to same Cyclomotic Field in same power base.
		ell_M2 = toCyclPowerBase(ell_M, ell_M_order)
		ell_R2 = toLowerCyclBase(ell_R, ell_R_order, ell_M_order)
		# We must work with the matrix. maybe we should transform hf_M instead to a
		# higher order field instead, if this ever fails (I'm not sure).
		assert ell_R2 is not None
		assert len(ell_M2) == len(ell_R2) # They should have the same power base & same degree now.
		print "ell_M2[0], ell_R2[0] rank with order %i:" % ell_M_order, ell_M2[0].rank(), ell_R2[0].rank()

		assert len(M_R) == len(ell_M2)
		for i in range(len(ell_M2)):
			ell_M_space = ell_M2[i].row_space()
			ell_R_space = ell_R2[i].column_space()
			merged = ell_M_space.intersection(ell_R_space)

			herm_modform_fe_expannsion_Ci = M_R[i].solve_right( merged.basis_matrix().transpose() )
			herm_modform_fe_expannsion_Ci_module = herm_modform_fe_expannsion_Ci.column_module()
			herm_modform_fe_expannsion_Ci_module += M_R[i].right_kernel()

			extra_check_on_herm_superspace(
				vs=herm_modform_fe_expannsion_Ci_module,
				D=D, B_cF=calc.B_cF, HermWeight=HermWeight
			)

			herm_modform_fe_expannsion = herm_modform_fe_expannsion.intersection( herm_modform_fe_expannsion_Ci_module )
			print "power", i, merged.dimension(), herm_modform_fe_expannsion_Ci_module.dimension(), \
				herm_modform_fe_expannsion.dimension()
			current_dimension = herm_modform_fe_expannsion.dimension()

	return herm_modform_fe_expannsion


@persistent_cache(name="modform_restriction_info")
def modform_restriction_info(calc, S, l):
	assert l == S.det()
	assert list(calc.curlS) == [S]

	D = calc.D
	HermWeight = calc.HermWeight
	B_cF = calc.B_cF

	# Step 4. Calculate restriction matrix. Via calc.calcMatrix() (algo_cpp.cpp).
	# Note that calcMatrix() depends on the current internal calc.curlS set
	# and on the internal calc.curlF. curlF only depends on B_cF which is not changed here.
	verbose("calc restriction matrix...")
	M_S = calcRestrictMatrix(calc) # matrix over integer ring
	M_S = M_S.matrix_over_field() # matrix over rational field

	# The maximum precision of Elliptic modular forms is given in
	# the text by \cF(S). This is also the number of rows of M_S.
	precLimit = M_S.nrows()
	assert precLimit == calcPrecisionDimension(B_cF=B_cF, S=S)

	# These are the Elliptic modular forms with weight 2*HermWeight to \Gamma_0(l).
	verbose("get elliptic modform space with precision %i ..." % precLimit)
	ell_dim, fe_expansion_matrix_l = getElliptModFormsBasisMatrix(l, 2*HermWeight, precLimit)
	if fe_expansion_matrix_l.rank() < ell_dim:
		verbose("ignoring ell modforms because matrix is not expressive enough. ell_dim=%i, matr rank=%i" % (ell_dim, fe_expansion_matrix_l.rank()))
		return None

	ell_modform_fe_expansions_l = fe_expansion_matrix_l.row_module()
	verbose("dim of elliptic modform space: %i" % ell_modform_fe_expansions_l.dimension())
	assert ell_modform_fe_expansions_l.dimension() == ell_dim

	verbose("calc M_S_module...")
	M_S_module = M_S.column_module()
	verbose("dimension of M_S column module: %i" % M_S_module.dimension())
	restriction_fe_expansions = ell_modform_fe_expansions_l.intersection( M_S_module )
	verbose("dimension of restriction_fe_expansions: %i" % restriction_fe_expansions.dimension())
	if restriction_fe_expansions.dimension() == 0:
		# This must not happen. Probably for more cases, but at least these:
		if D == -3 and HermWeight == 6:
			assert False
	herm_modform_fe_expannsion_S = M_S.solve_right( restriction_fe_expansions.basis_matrix().transpose() )
	herm_modform_fe_expannsion_S_module = herm_modform_fe_expannsion_S.column_module()
	verbose("dimension of herm column module: %i" % herm_modform_fe_expannsion_S_module.dimension())
	verbose("calc M_S_right_kernel...")
	M_S_right_kernel = M_S.right_kernel()
	verbose("dimension of M_S right kernel: %i" % M_S_right_kernel.dimension())
	herm_modform_fe_expannsion_S_module += M_S_right_kernel

	try:
		extra_check_on_herm_superspace(
			vs=herm_modform_fe_expannsion_S_module,
			D=D, B_cF=B_cF, HermWeight=HermWeight
		)
	except Exception:
		print "restriction_fe_expansions =", restriction_fe_expansions
		print "M_S_right_kernel =", M_S_right_kernel
		print "herm_modform_fe_expannsion_S_module =", herm_modform_fe_expannsion_S_module
		raise

	return herm_modform_fe_expannsion_S_module


class IntersectSpacesTask:
	def __init__(self, basespace, spaces):
		verbose("queued intersection of %i spaces" % len(spaces))
		self.basespace = basespace
		self.spaces = spaces
	def __call__(self):
		herm_modform_fe_expannsion = self.basespace
		for descr, space in self.spaces:
			verbose("intersecting %r..." % descr)
			herm_modform_fe_expannsion = herm_modform_fe_expannsion.intersection( space )
			current_dimension = herm_modform_fe_expannsion.dimension()
			verbose("current dimension: %i" % current_dimension)
		if herm_modform_fe_expannsion.dimension() == self.basespace.dimension():
			# No gain, just ignore this intersection.
			return None
		return herm_modform_fe_expannsion

class CalcTask:
	def __init__(self, func, calc, kwargs):
		# This will call func(calc, **kwargs) later.
		self.func = func
		self.calc_state = calc.__getstate__()
		self.kwargs = kwargs
	def __call__(self):
		calc = C.Calc()
		calc.__setstate__(self.calc_state)
		return self.func(calc, **self.kwargs)


# These are flags for the parameter `reduction_method_flags` of `herm_modform_space`.
Method_Elliptic_reduction = 1 << 0
Method_EllipticCusp_reduction = 1 << 1


hermModformSpaceCache = PersistentCache("herm_modform_space__precalc")
def herm_modform_space(D, HermWeight, B_cF=10, parallelization=None, reduction_method_flags=-1):
	"""
	This calculates the vectorspace of Fourier expansions to
	Hermitian modular forms of weight `HermWeight` over \Gamma,
	where \Gamma = \Sp_2(\curlO) and \curlO is the maximal order
	of \QQ(\sqrt{D}).

	Each Fourier coefficient vector is indexed up to a precision
	\curlF which is given by `B_cF` such that for every
	[a,b,c] \in \curlF \subset \Lambda, we have 0 \le a,c \le B_cF.

	The function `herm_modform_indexset()` returns reduced matrices
	of that precision index set \curlF.
	"""

	if HermWeight % 3 != 0:
		raise TypeError, "the modulform is trivial/zero if HermWeight is not divisible by 3"

	# Transform these into native Python objects. We don't want to have
	# any Sage objects (e.g. Integer) here so that the cache index stays
	# unique.
	D = int(D)
	HermWeight = int(HermWeight)
	B_cF = int(B_cF)
	reduction_method_flags = int(reduction_method_flags)

	calc = C.Calc()
	calc.init(D = D, HermWeight = HermWeight, B_cF=B_cF)
	calc.calcReducedCurlF()
	reducedCurlFSize = calc.matrixColumnCount

	# Calculate the dimension of Hermitian modular form space.
	dim = herm_modform_space_dim(D=D, HermWeight=HermWeight)

	cacheIdx = (D, HermWeight, B_cF)
	if reduction_method_flags != -1:
		cacheIdx += (reduction_method_flags,)
	try:
		herm_modform_fe_expannsion, calc, curlS_denoms, pending_tasks = hermModformSpaceCache[cacheIdx]
		if not isinstance(calc, C.Calc): raise TypeError
		print "Resuming from %s" % hermModformSpaceCache._filename_for_key(cacheIdx)
	except (TypeError, ValueError, KeyError, EOFError): # old format or not cached or cache incomplete
		herm_modform_fe_expannsion = FreeModule(QQ, reducedCurlFSize)
		curlS_denoms = set() # the denominators of the visited matrices S
		pending_tasks = ()

	current_dimension = herm_modform_fe_expannsion.dimension()

	verbose("current dimension: %i, wanted: %i" % (herm_modform_fe_expannsion.dimension(), dim))
	if dim == 0:
		print "dim == 0 -> exit"
		return

	def task_iter_func():
		# Iterate S \in Mat_2^T(\curlO), S > 0.
		while True:
			# Get the next S.
			# If calc.curlS is not empty, this is because we have recovered from a resume.
			if len(calc.curlS) == 0:
				S = calc.getNextS()
			else:
				assert len(calc.curlS) == 1
				S = calc.curlS[0]

			l = S.det()
			l = toInt(l)
			curlS_denoms.add(l)

			verbose("trying S={0}, det={1}".format(S, l))

			if reduction_method_flags & Method_Elliptic_reduction:
				yield CalcTask(
					func=modform_restriction_info, calc=calc, kwargs={"S":S, "l":l})

			if reduction_method_flags & Method_EllipticCusp_reduction:
				precLimit = calcPrecisionDimension(B_cF=B_cF, S=S)
				yield CalcTask(
					func=modform_cusp_info, calc=calc, kwargs={"S":S, "l":l, "precLimit": precLimit})

			calc.curlS_clearMatrices() # In the C++ internal curlS, clear previous matrices.

	task_iter = task_iter_func()
	if parallelization:
		parallelization.task_iter = task_iter

	if parallelization and pending_tasks:
		for func, name in pending_tasks:
			parallelization.exec_task(func=func, name=name)

	step_counter = 0
	while True:
		if parallelization:
			new_task_count = 0

			spaces = []
			for task, exc, newspace in parallelization.get_all_ready_results():
				if exc: raise exc

				if newspace is None:
					verbose("no data from %r" % task)
					continue
				if newspace.dimension() == reducedCurlFSize:
					verbose("no information gain from %r" % task)
					continue

				spacecomment = task
				assert newspace.dimension() >= dim, "%r, %r" % (task, newspace)
				if newspace.dimension() < herm_modform_fe_expannsion.dimension():
					# Swap newspace with herm_modform_fe_expannsion.
					herm_modform_fe_expannsion, newspace = newspace, herm_modform_fe_expannsion
					current_dimension = herm_modform_fe_expannsion.dimension()
					if current_dimension == dim:
						if not isinstance(task, IntersectSpacesTask):
							verbose("warning: we expected IntersectSpacesTask for final dim but got: %r" % task)
					verbose("new dimension: %i, wanted: %i" % (current_dimension, dim))
					if current_dimension <= 20:
						pprint(herm_modform_fe_expannsion.basis())
					spacecomment = "<old base space>"

				spaces += [(spacecomment, newspace)]

			if spaces:
				parallelization.exec_task(IntersectSpacesTask(herm_modform_fe_expannsion, spaces))
				new_task_count += 1

			new_task_count += parallelization.maybe_queue_tasks()
			time.sleep(0.1)

		else: # no parallelization
			new_task_count = 1
			if pending_tasks: # from some resuming
				task,_ = pending_tasks[0]
				pending_tasks = pending_tasks[1:]
			else:
				task = next(task_iter)
			newspace = task()

			if newspace is None:
				verbose("no data from %r" % task)
			if newspace is not None and newspace.dimension() == reducedCurlFSize:
				verbose("no information gain from %r" % task)
				newspace = None

			if newspace is not None:
				new_task_count += 1
				spacecomment = task
				herm_modform_fe_expannsion_new = IntersectSpacesTask(
					herm_modform_fe_expannsion, [(spacecomment, newspace)])()
				if herm_modform_fe_expannsion_new is not None:
					herm_modform_fe_expannsion = herm_modform_fe_expannsion_new
				current_dimension = herm_modform_fe_expannsion.dimension()
				verbose("new dimension: %i, wanted: %i" % (current_dimension, dim))

		if new_task_count > 0:
			step_counter += 1

			if step_counter % 10 == 0:
				verbose("save state after %i steps to %s" % (step_counter, os.path.basename(hermModformSpaceCache._filename_for_key(cacheIdx))))
				if parallelization:
					pending_tasks = parallelization.get_pending_tasks()
				hermModformSpaceCache[cacheIdx] = (herm_modform_fe_expannsion, calc, curlS_denoms, pending_tasks)

		if current_dimension == dim:
			verbose("finished!")
			break

	# Test for some other S with other not-yet-seen denominator.
	check_herm_modform_space(
		calc, herm_modform_space=herm_modform_fe_expannsion,
		used_curlS_denoms=curlS_denoms
		)

	return herm_modform_fe_expannsion


def herm_modform_space__parallel(D, HermWeight, B_cF=10, task_limit=4, **kwargs):
	parallelization = Parallelization(task_limit=task_limit)
	return herm_modform_space(
		D=D, HermWeight=HermWeight, B_cF=B_cF,
		parallelization=parallelization,
		**kwargs)
