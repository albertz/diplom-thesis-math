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
inline int compare(const M2T& m1, const M2T& m2) {
	if(m1.a != m2.a) return (m1.a < m2.a) ? -1 : 1;
	if(m1.b1 != m2.b1) return (m1.b1 < m2.b1) ? -1 : 1;
	if(m1.b2 != m2.b2) return (m1.b2 < m2.b2) ? -1 : 1;
	if(m1.c != m2.c) return (m1.c < m2.c) ? -1 : 1;
	return 0;
}
inline bool operator==(const M2T& m1, const M2T& m2) {
	return compare(m1, m2) == 0;
}
inline bool operator!=(const M2T& m1, const M2T& m2) {
	return compare(m1, m2) != 0;
}
inline bool operator<(const M2T& m1, const M2T& m2) {
	return compare(m1, m2) < 0;
}

template<typename T>
struct Matrix2 {
	T a,b,c,d; // [[a,b],[c,d]]
	Matrix2(T _a = 0, T _b = 0, T _c = 0, T _d = 0)
	: a(_a), b(_b), c(_c), d(_d) {}
	T det() { return a*c - b*d; }
};



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
