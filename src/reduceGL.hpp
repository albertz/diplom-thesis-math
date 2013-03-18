
#ifndef __SAGE_HERMITIAN_REDUCEGL_HPP__
#define __SAGE_HERMITIAN_REDUCEGL_HPP__

// https://github.com/martinra/psage/blob/paper_computing_jacobi_forms/psage/modform/hermitianmodularforms/hermitianmodularformd2_fourierexpansion_cython.pyx

#include <assert.h>
#include <math.h>
#include <iostream>
#include <stdlib.h>
#include "structs.hpp"

struct reduce_character_evalutation {
	// det_character is alread mod h, i.e. < h.
	// det = exp(2 pi i det_character / h)
	// where h = 2, or if D = -3, then h = 6, or D = -4, then h = 4.

	int transposition;
	int determinant;
	int nu;
	reduce_character_evalutation(int _trans = 0, int _det = 0, int _nu = 0)
	: transposition(_trans), determinant(_det), nu(_nu) {}
	
	int detValue(const int D) {
		// det = exp(2 pi i det_character / h)
		// where h = 2, or if D = -3, then h = 6, or D = -4, then h = 4
		// we expect that it is a unit in \Z here.
		const int h = (D == -3) ? 6 : (D == -4) ? 4 : 2;
		if(Mod(determinant, h) == 0) return 1;
		if(Mod(determinant, h) == h/2) return -1;
		std::cerr << "reduce_character_evalutation.detValue: determinant = " << determinant << ", D = " << D << std::endl;
		abort();
		return 0;
	}
};
inline std::ostream& operator<< (std::ostream& os, const reduce_character_evalutation& c) {
	return os << "_C(" << c.transposition << ", " << c.determinant << ", " << c.nu << ")";
}
inline bool operator==(const reduce_character_evalutation& c1, const reduce_character_evalutation& c2) {
	return c1.transposition == c2.transposition && c1.determinant == c2.determinant && c1.nu == c2.nu;
}
inline bool operator!=(const reduce_character_evalutation& c1, const reduce_character_evalutation& c2) {
	return !(c1 == c2);
}

struct hermitian_form_with_character_evaluation {
	M2T matrix;
	reduce_character_evalutation character;
	hermitian_form_with_character_evaluation(M2T _m = M2T(), reduce_character_evalutation _char = reduce_character_evalutation())
	: matrix(_m), character(_char) {}
};
inline std::ostream& operator<< (std::ostream& os, const hermitian_form_with_character_evaluation& f) {
	return os << "_H(" << f.matrix << ", " << f.character << ")";
}
inline bool operator==(const hermitian_form_with_character_evaluation& f1, const hermitian_form_with_character_evaluation& f2) {
	return f1.matrix == f2.matrix && f1.character == f2.character;
}
inline bool operator!=(const hermitian_form_with_character_evaluation& f1, const hermitian_form_with_character_evaluation& f2) {
	return !(f1 == f2);
}


inline void reduce_GL(M2T matrix, int D, struct hermitian_form_with_character_evaluation& res) {
	auto a = matrix.a, b1 = matrix.b1, b2 = matrix.b2, c = matrix.c;
	
	/*
	 Reduce the positive semi-definite hermitian quatratic form `\lbrackt a,b,c \rbrackt`
	 with respect to `\GL(2, \mathfrak{o}_{\Q(\sqrt{D})})` and the transposition.
	 Return the evaluation of the characters `det` and `\nu` (see Dern where
	 he calls this `\nu_\mathcal{p}`) as well as the number mod 2 of applied
	 transpositions.
	 
	 INPUT:
	 `s` -- A tuple `(a, b1, b2, c)`; We set `b = b1 / \sqrt{D} + b2 (1 + \sqrt{D})/2`.
	 `D` -- An negative integer; The fundamental discriminant of the underlying;
	 imaginary quadratic number field.
	 
	 OUTPUT:
	 A pair `(s, (\mathrm{trans}, \mathrm{det}, \nu))` of a quadratic form `s`,
	 which is represented by a four-tuple, and character evaluations, with
	 `\nu` and `\mathrm{trans}` either `1` or `-1` and `\mathrm{det}` between
	 `0` and `|\mathfrak{o}_{\Q(\sqrt{D})}^\times|`.
	 
	 TEST::
	 sage: from hermitianmodularforms.hermitianmodularformd2_fourierexpansion_cython import reduce_GL
	 sage: reduce_GL((2,1,1,2), -3)
	 ((2, 1, 1, 2), (1, 0, 1))
	 sage: reduce_GL((1,0,0,-1), -3)
	 Traceback (most recent call last):
	 ...
	 NotImplementedError: only implemented for non-positive discriminants: (1, 0, 0, -1)
	 sage: reduce_GL((-1,0,0,1), -3)
	 Traceback (most recent call last):
	 ...
	 NotImplementedError: only implemented for non-positive discriminants: (-1, 0, 0, 1)
	 sage: reduce_GL((-1,0,0,-1), -3)
	 Traceback (most recent call last):
	 ...
	 NotImplementedError: only implemented for non-positive discriminants: (-1, 0, 0, -1)
	 sage: reduce_GL((1,6,0,1), -3)
	 Traceback (most recent call last):
	 ...
	 NotImplementedError: only implemented for non-positive discriminants: (1, 6, 0, 1)
	 sage: reduce_GL((1,50,0,1), -7)
	 Traceback (most recent call last):
	 ...
	 NotImplementedError: only implemented for non-positive discriminants: (1, 50, 0, 1)
	 */
	
	/*
	Reduce the positive semi-definite hermitian quatratic form
	`\lbrackt a,b,c \rbrackt` with respect to `\mathrm{GL}_2(\mathcal{o}_D)`.
	
	A reduced form `\lbrackt a,b,c \rbrackt` satisfies `a \le c`,
	`|\Im(b)| / \sqrt{-D} \le a / 4`, `|\Re(b)| \le a / 2`. It satisfies `b2 \ge 0`.
	If `D = -3` it also satisfies `-D b2 \ge 2 b1 \ge 2 b2`.
	If `D = -4` it also satisfies `2 b2 \ge b1 \ge b2`.
	In all other cases it also satisfies `2 b1 + D b2 \le 0`.
	*/

	// the discriminant will be -D det(M)
	// FIXME: the discriminant can become too big
	if( D >= 0 ) {
		std::cerr << "reduce_GL: we expect negative D. D = " << D << std::endl;
		abort();
	}
	if( a < 0 || c < 0 || matrix.det4D(D) < 0 ) {
		std::cerr << "reduce_GL: invalid input: " << matrix << ", " << D << std::endl;
		abort();
	}

	int q, r;
	int tmp;
	// for a describtion of det, trans and nu see the output documentation
	// of reduce_GL
	int det = 0;
	int trans = 1;
	int nu = 1;

	
	// K = QQ[\sqrt -1]
	if(D == -4) {
		// We want b to be of from the positive real axis by at most pi/4 in
		// positive direction
	
		// if b is not in the right half plane;
		if(b2 < 0) {
			// apply [[-1, 0], [0, 1]]
			b1 = -b1;
			b2 = -b2;
			
			det = det + 2;
		}
		
		// if b is in the upper quarter of the right half plane
		if(b1 < b2) {
			// apply [[i, 0], [0, 1]]
			tmp = -2 * b1 + 5 * b2;
			b2 = -b1 + 2 * b2;
			b1 = tmp;
			
			det = det + 1;
			
			// apply transposition
			b1 = -b1 + 4 * b2;
			
			trans = -trans;
		}
		
		// if b is in the lower quarter of the right half plane
		else if(b1 > 3 * b2) {
			// apply [[-i, 0], [0, 1]]
			tmp = 2 * b1 - 5 * b2;
			b2 = b1 - 2 * b2;
			b1 = tmp;
			
			det = det + 3;
		}
		
		// b is in the second quarter counted from below of the right half plane
		else if(b1 > 2 * b2) {
			// apply transposition
			b1 = -b1 + 4 * b2;
			
			trans = -trans;
		}
	}
	
	// K = QQ[\sqrt -3]
	else if(D == -3) {
		if(b2 < 0) {
			// apply [[-1, 0], [0, 1]]
			b1 = -b1;
			b2 = -b2;
			
			det = det + 3;
		}
		
		// TODO: Use one transformation for each sixth of the left half plane
		// check whether b is between 0 and pi/6 off from the real line in mathematical positive direction.
		// b is in the upper quadrant
		if(b1 < 0) {
			// apply [[\bar eps, 0], [0, 1]]
			tmp = -b1 + 3*b2;
			b2 = -b1 + 2*b2;
			b1 = tmp;
		
			det = det + 5;
		}
		
		// b is in the last third of the lower quadrant
		else if(b1 > 3 * b2) {
			// apply [[eps**2, 0], [0, 1]]
			tmp = -3 * b2 + b1;
			b2 = b1 - 2 * b2;
			b1 = tmp;
			
			det = det + 2;
		}
		
		// b is in the lower quadrant
		else if(2 * b1 > 3 * b2) {
			// apply [[eps, 0], [0, 1]]
			tmp = 2*b1 - 3*b2;
			b2 = b1 - b2;
			b1 = tmp;
			
			det = det + 1;
		}
		
		// now b is between 0 and pi/3 off from the real line in mathematical positive direction.
		// check whether b is more than pi/6 off from the real in mathematical positive direction.
		if(b1 < b2) {
			// apply the transpose and
			// apply [[eps, 0], [0, 1]]
			tmp = 3 * b2 - 2 * b1;
			b2 = 2 * b2 - b1;
			b1 = tmp;

			trans = -trans;
			det = det + 1;
		}
		
	// ! else if(D == -3)
	} else {
		// if b is not in the right half plane
		if(b2 < 0) {
			// apply [[-1, 0], [0, 1]]
			b1 = -b1;
			b2 = -b2;
			
			det = det + 1;
		}
		
		// if b is not in the upper half plane;
		if(2 * b1 + D * b2 > 0) {
			b1 = - b1 - D * b2;
			
			trans = -trans;
		}
	} //! else D == -4

	// while not GL-reduced
	while (! ((a == 0) || ( a<=c && abs(b2) <= a && abs(-4*b1 - 2 * D * b2) <= -D*a )) ) {
		// apply symmetric step (reflection) if necessary
		if(b1 < 0) {
			// apply [[-1, 0], [0, 1]]
			b1 = -b1;
			b2 = -b2;
			if(D == -4)
				det = det + 2;
			else if(D == -3)
				det = det + 3;
			else
				det = det + 1;
		}
		
		// abs(Im b / sqrt D) <= a/4 <=> -4*D abs(Im b) <= -D * a
		if(abs(-4*b1 - 2 * D * b2) > -D*a) {
			q = Div((-2*b1 - D * b2), (-D * a));
			r = Mod((-2*b1 - D * b2), (-D * a));
			
			if(r > (-D * a) / 2) {
				// r -= -D * a;
				q += 1;
			}
			
			// apply [[1, -q*(D + \sqrt D)/2], [0, 1]]
			c = c + q*q * a * Div((D*D - D), 4) + b1 * q;
			b1 = b1 + q * a * Div((D-1) * D, 2);
			b2 = b2 - q * a * D;
		}
		
		// abs(Re b) <= a/2
		if(abs (b2) > a) {
			q = Div(b2, (2 * a));
			r = Mod(b2, (2 * a));

			if(r > a) {
				r = r - 2 * a;
				q = q + 1;
			}
			
			// apply [[1, -q], [0, 1]]
			c = c - b2 * q + a * q*q;
			b1 = b1 + D * q * a ;
			b2 = r; // = b2 - 2 * q * a;
		}
		
		if(a > c) {
			// apply [[0, 1], [1, 0]]
			tmp = c;
			c = a;
			a = tmp;
			b1 = -b1 - D * b2;
			if(D == -4)
				det = det + 2;
			else if(D == -3)
				det = det + 3;
			else
				det = det + 1;
		}
		
		// K = QQ[\sqrt -1]
		if(D == -4) {
			// We want b to be of from the positive real axis by at most pi/4 in
			// positive direction
		
			// if b is not in the right half plane
			if(b2 < 0) {
				// apply [[-1, 0], [0, 1]]
				b1 = -b1;
				b2 = -b2;
				
				det = det + 2;
			}
			
			// if b is in the upper quarter of the right half plane
			if(b1 < b2) {
				// apply [[i, 0], [0, 1]]
				tmp = -2 * b1 + 5 * b2;
				b2 = -b1 + 2 * b2;
				b1 = tmp;
				
				det = det + 1;
				
				// apply transposition
				b1 = -b1 + 4 * b2;
				
				trans = -trans;
			}
			
			// if b is in the lower quarter of the right half plane
			else if(b1 > 3 * b2) {
				// apply [[-i, 0], [0, 1]]
				tmp = 2 * b1 - 5 * b2;
				b2 = b1 - 2 * b2;
				b1 = tmp;
				
				det = det + 3;
			}
			
			// b is in the second quarter counted from below of the right half plane
			else if(b1 > 2 * b2) {
				// apply transposition
				b1 = -b1 + 4 * b2;
				
				trans = -trans;
			}
		}
		
		// K = QQ[\sqrt -3]
		else if(D == -3) {
			if(b2 < 0) {
				// apply [[-1, 0], [0, 1]]
				b1 = -b1;
				b2 = -b2;
				
				det = det + 3;
			}
	
			// TODO: Use one transformation for each sixth of the left half plane
			// check whether b is between 0 and pi/6 off from the real line in mathematical positive direction.
			// b is in the upper quadrant
			if(b1 < 0) {
				// apply [[\bar eps, 0], [0, 1]]
				tmp = -b1 + 3*b2;
				b2 = -b1 + 2*b2;
				b1 = tmp;
				
				det = det + 5;
			}
			
			// b is in the last third of the lower quadrant
			else if(b1 > 3 * b2) {
				// apply [[eps**2, 0], [0, 1]]
				tmp = -3 * b2 + b1;
				b2 = b1 - 2 * b2;
				b1 = tmp;
				
				det = det + 2;
			}
			
			// b is in the lower quadrant
			else if(2 * b1 > 3 * b2) {
				// apply [[eps, 0], [0, 1]]
				tmp = 2*b1 - 3*b2;
				b2 = b1 - b2;
				b1 = tmp;
				
				det = det + 1;
			}
			
			// now b is between 0 and pi/3 off from the real line in mathematical positive direction.
			// check whether b is more than pi/6 off from the real in mathematical positive direction.
			if(b1 < b2) {
				// apply the transpose and
				// apply [[eps, 0], [0, 1]]
				tmp = 3 * b2 - 2 * b1;
				b2 = 2 * b2 - b1;
				b1 = tmp;
				
				trans = -trans;
				det = det + 1;
			}
		} //! elif D == -3;
		else {
			// if b is not in the right half plane
			if(b2 < 0) {
				// apply [[-1, 0], [0, 1]]
				b1 = -b1;
				b2 = -b2;
				
				det = det + 1;
			}
			
			// if b is not in the upper half plane
			if(2 * b1 + D * b2 > 0) {
				b1 = - b1 - D * b2;
				
				trans = -trans;
			}
		}
		//! else D == -4
	
	} // #! while
	
	res.matrix.a = a;
	res.matrix.b1 = b1;
	res.matrix.b2 = b2;
	res.matrix.c = c;
	res.character.transposition = trans;
	if(D == -4)
		res.character.determinant = Mod(det, 4);
	else if(D == -3)
		res.character.determinant = Mod(det, 6);
	else
		res.character.determinant = Mod(det, 2);
	res.character.nu = nu;
	
}

#endif
