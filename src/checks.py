# -*- coding: utf-8 -*-
# Hermitian modular forms, https://github.com/albertz/diplom-thesis-math
# Copyright (c) 2013, Albert Zeyer, www.az2000.de
# This code is under the GPL v3 or later, see License.txt in the root directory of this project.

from sage.matrix.constructor import matrix
from sage.modules.free_module_element import vector
from sage.rings.number_field.number_field import QQ


def check_eisenstein_series_D3_weight6(vs, B_cF):
	from algo import herm_modform_indexset
	from helpers import M2T_Odual
	D = -3
	HermWeight = 6
	indexset = herm_modform_indexset(D=D, B_cF=B_cF)
	jacobi_coeffs_1 = {
		(0, (0, 0)): 1,
		(1, (0, 0)): -240,
		(1, (1, 1)): -45,
		(2, (0, 0)): -3690,
		(2, (1, 1)): -1872,
		(3, (0, 0)): -19680,
		(3, (1, 1)): -11565,
		(4, (0, 0)): -57840,
		(4, (1, 1)): -43920
	}

	from reduceGL import reduce_GL
	def reduce_matrix_jacobi1(key, value):
		# In standard base of \cO^#.
		a, (b1, b2) = key
		m, char = reduce_GL((a, b1, b2, 1), D)
		# Like reduce_character_evalutation::value() in reduceGL.hpp.
		h = 6 # because D == -3
		k = -HermWeight
		detchar = char[1] * k
		if detchar % h == 0: detvalue = 1
		elif detchar % h == h/2: detvalue = -1
		else: assert False, "detchar %i" % detchar
		m = M2T_Odual(m, D)
		return m, value * detvalue
	jacobi_coeffs_1_transformed = dict(
		[reduce_matrix_jacobi1(key,value) for (key,value) in jacobi_coeffs_1.items()])

	indexmap = {}
	for m in jacobi_coeffs_1_transformed.keys():
		indexmap[m] = indexset.index(m)
	reverse_indexlist = sorted([(value,key) for (key,value) in indexmap.items()])

	def reduced_vector(v):
		new_v = [None] * len(indexmap)
		for i,(j,m) in enumerate(reverse_indexlist):
			new_v[i] = v[j]
		return vector(QQ, new_v)
	reduced_vs_basis = map(reduced_vector, vs.basis())
	reduced_vs_basis_matrix = matrix(QQ, reduced_vs_basis)

	vector_jac = [None] * len(indexmap)
	for i,(j,m) in enumerate(reverse_indexlist):
		vector_jac[i] = jacobi_coeffs_1_transformed[m]
	vector_jac = vector(QQ, vector_jac)

	try:
		lincomb = reduced_vs_basis_matrix.solve_left(vector_jac)
	except ValueError as e:
		if "no solutions" in str(e):
			print "reduced_vs_basis_matrix =", reduced_vs_basis_matrix, ", rank =", reduced_vs_basis_matrix.rank()
			print "vector_jac =", vector_jac
			assert False
		raise
	return lincomb


def extra_check_on_herm_superspace(vs, D, HermWeight, B_cF):
	# This check seems to be wrong.
	#if D == -3 and HermWeight == 6:
	#	check_eisenstein_series_D3_weight6(vs=vs, B_cF=B_cF)
	pass


def check_herm_modform_space(calc, herm_modform_space, used_curlS_denoms, checkSCount = 10):
	"""
	It uses the C++ calc structure to search for additional S matrices
	which have other denominators than those in used_curlS_denoms.
	For testScount matrices with unique denominators, it calculates
	the restriction matrix via the C++ calc structure for f \mapsto f[S].
	When mapping herm_modform_space, we must only get Elliptic modular forms.
	We check whether they are in ModularForms(\Gamma_0(l), prec).
	"""

	from helpers import calcRestrictMatrix, getElliptModFormsBasisMatrix, toInt
	from algo import verbose

	HermWeight = calc.HermWeight
	curlS_denoms = set(used_curlS_denoms)

	while checkSCount > 0:
		calc.curlS_clearMatrices()
		S = calc.getNextS()
		l = S.det()
		l = toInt(l)
		if l in curlS_denoms: continue

		curlS_denoms.add(l)
		checkSCount -= 1
		verbose("testing with S={0}, det={1}".format(S, l))

		verbose("calc restriction matrix...")
		M_S = calcRestrictMatrix(calc) # matrix over integer ring
		M_S = M_S.matrix_over_field() # matrix over rational field

		precLimit = M_S.nrows() # \cF(S)

		# These are the Elliptic modular forms with weight 2*HermWeight to \Gamma_0(l).
		verbose("get elliptic modform space with precision %i ..." % precLimit)
		ell_dim, fe_expansion_matrix_l = getElliptModFormsBasisMatrix(l, 2*HermWeight, precLimit)
		if fe_expansion_matrix_l.rank() < ell_dim:
			verbose("ignoring ell modforms because matrix is not expressive enough")
			checkSCount += 1
			continue
		ell_modform_fe_expansions_l = fe_expansion_matrix_l.row_module()

		verbose("calc M_S * herm_modforms ...")
		m = M_S * herm_modform_space.basis_matrix().transpose()

		m_module = m.column_module()
		assert m_module.is_subspace(ell_modform_fe_expansions_l), \
			"%r not subspace of %r" % (m_module, ell_modform_fe_expansions_l)