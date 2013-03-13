#ifndef __SAGE_HERMITIAN_STRUCTS_HPP__
#define __SAGE_HERMITIAN_STRUCTS_HPP__

#include <iostream>
#include <string.h>
#include <assert.h>

typedef int Int;

struct M2T {
	// [[a,b1 + ib2],[b1 - ib2, c]]
	Int a, b1, b2, c;
	M2T(Int _a = 0, Int _b1 = 0, Int _b2 = 0, Int _c = 0)
	: a(_a), b1(_b1), b2(_b2), c(_c) {}
};
inline std::ostream& operator<<(std::ostream& os, const M2T& m) {
	return os << "M2T(" << m.a << "," << m.b1 << "," << m.b2 << "," << m.c << ")";
}
inline bool operator==(const M2T& m1, const M2T& m2) {
	return m1.a == m2.a && m1.b1 == m2.b1 && m1.b2 == m2.b2 && m1.c == m2.c;
}
inline bool operator!=(const M2T& m1, const M2T& m2) {
	return !(m1 == m2);
}


template<typename T>
T Mod(const T& a, const T& b) {
	T res = a % b;
	if(a < 0) res += b;
	assert(res >= 0);
	return res;
}

template<typename T>
T Div(const T& a, const T& b) {
	T res = a / b;
	if(a < 0) res -= 1;
	return res;
}

template<typename T>
T Pow(const T& a, const T& b) {
	if(b == 0) return 1;
	if(b < 0) return 1 / Pow(a, -b);
	if(b == 1) return a;
	T result = Pow(a, b / 2);
	result *= result;
	if(b % 2) result *= a;
	return result;
}


#endif
