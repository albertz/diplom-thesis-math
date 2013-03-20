#include "reduceGL.hpp"
#include "structs.hpp"
#include <iostream>

typedef reduce_character_evalutation _C;
typedef hermitian_form_with_character_evaluation _H;
using namespace std;

void reduce_GL_assertEqual(M2T_O matrix, int D, _H result) {
	_H calcRes;
	cout << "reduce_GL(" << matrix << ", " << D << ") == " << result << endl;
	reduce_GL(matrix, D, calcRes);
	if(result != calcRes) {
		cerr << "not equal: " << calcRes << endl;
		abort();
	}
}

void test_reduceGL() {
	reduce_GL_assertEqual(M2T_O(2,-2,2,5), -4,
		_H(M2T_O(0, 0, 0, 1), _C(1, 1, 1)));
	reduce_GL_assertEqual(M2T_O(0,0,0,1), -4,
		_H(M2T_O(0, 0, 0, 1), _C(1, 0, 1)));
	reduce_GL_assertEqual(M2T_O(10,1,4,2), -3,
		_H(M2T_O(2, 1, 1, 4), _C(1, 4, 1)));
	reduce_GL_assertEqual(M2T_O(1,2,0,2), -3,
		_H(M2T_O(1, 1, 1, 1), _C(1, 3, 1)));
	reduce_GL_assertEqual(M2T_O(1,0,0,2), -3,
		_H(M2T_O(1, 0, 0, 2), _C(1, 0, 1)));
	reduce_GL_assertEqual(M2T_O(2,0,0,1), -3,
		_H(M2T_O(1, 0, 0, 2), _C(1, 3, 1)));
	reduce_GL_assertEqual(M2T_O(2,1,2,2), -3,
		_H(M2T_O(1, 1, 1, 2), _C(-1, 1, 1)));
	reduce_GL_assertEqual(M2T_O(2,3,1,2), -3,
		_H(M2T_O(2, 3, 2, 2), _C(1, 1, 1)));
	reduce_GL_assertEqual(M2T_O(1,1,0,1), -3,
		_H(M2T_O(1, 1, 1, 1), _C(1, 2, 1)));
	reduce_GL_assertEqual(M2T_O(10,7, 8, 9), -3,
		_H(M2T_O(9, 10, 9, 10), _C(1, 0, 1)));
	reduce_GL_assertEqual(M2T_O(50, 76, 30, 19), -3,
		_H(M2T_O(19, 14, 11, 23), _C(1, 1, 1)));
	reduce_GL_assertEqual(M2T_O(1,0,0,2), -4,
		_H(M2T_O(1, 0, 0, 2), _C(1, 0, 1)));
	reduce_GL_assertEqual(M2T_O(1,2,0,2), -4,
		_H(M2T_O(1, 0, 0, 1), _C(1, 3, 1)));
	reduce_GL_assertEqual(M2T_O(1,3,1,1), -4,
		_H(M2T_O(1, 1, 1, 1), _C(-1, 0, 1)));
	reduce_GL_assertEqual(M2T_O(17,0,5,10), -4,
		_H(M2T_O(10, 15, 10, 17), _C(1, 3, 1)));
	reduce_GL_assertEqual(M2T_O(117,43,27,10), -4,
		_H(M2T_O(10, 11, 9, 99), _C(-1, 1, 1)));
	reduce_GL_assertEqual(M2T_O(1,0,0,2), -7,
		_H(M2T_O(1, 0, 0, 2), _C(1, 0, 1)));
	reduce_GL_assertEqual(M2T_O(1,2,0,2), -7,
		_H(M2T_O(1, 2, 1, 2), _C(-1, 1, 1)));
	reduce_GL_assertEqual(M2T_O(17,0,5,10), -7,
		_H(M2T_O(10, 0, 5, 17), _C(-1, 1, 1)));
	reduce_GL_assertEqual(M2T_O(117,43,27,10), -7,
		_H(M2T_O(10, -6, 3, 65), _C(1, 0, 1)));
}

