#ifndef __SAGE_HERMITIAN_STRUCTS_HPP__
#define __SAGE_HERMITIAN_STRUCTS_HPP__

#include <iostream>
#include <string.h>
#include <assert.h>

typedef int Int;

struct M2T_O {
	// We set `b = b1 / \sqrt{D} + b2 (1 + \sqrt{D})/2`, where
	// D is a negative integer, the fundamental discriminant of
	// the underlying imaginary quadratic number field.
	// This M2T_O-struct represents the matrix [a,b,c].
	Int a, b1, b2, c;
	M2T_O(Int _a = 0, Int _b1 = 0, Int _b2 = 0, Int _c = 0)
	: a(_a), b1(_b1), b2(_b2), c(_c) {}
	// det4D == -D * 4 * det
	// TODO(?): if D is fundamental, we always have 4|(D*D-D), thus we could just use -D * det.
	Int det4D(const int D) const {
		assert(D < 0);
		// det = a*c - |b|^2
		// Re(b) = 1/2 b2
		// Re(b)^2 = 1/4 b2^2
		// Im(b) = -b1/sqrt{D} + 1/2 \sqrt{D} b2
		// Im(b)^2 = b1^2/(-D) - b1*b2 + 1/4 (-D) b2^2
		// -> 4*(-D)*|b∫^2 = 4*b1^2 - 4*(-D)*b1*b2 + (D^2 - D)*b2^2
		return a*c*4*(-D) - 4*b1*b1 + 4*(-D)*b1*b2 - (D*D - D)*b2*b2;
	}
};
inline std::ostream& operator<<(std::ostream& os, const M2T_O& m) {
	return os << "M2T_O(" << m.a << "," << m.b1 << "," << m.b2 << "," << m.c << ")";
}
inline int compare(const M2T_O& m1, const M2T_O& m2) {
	if(m1.a != m2.a) return (m1.a < m2.a) ? -1 : 1;
	if(m1.b1 != m2.b1) return (m1.b1 < m2.b1) ? -1 : 1;
	if(m1.b2 != m2.b2) return (m1.b2 < m2.b2) ? -1 : 1;
	if(m1.c != m2.c) return (m1.c < m2.c) ? -1 : 1;
	return 0;
}
inline bool operator==(const M2T_O& m1, const M2T_O& m2) {
	return compare(m1, m2) == 0;
}
inline bool operator!=(const M2T_O& m1, const M2T_O& m2) {
	return compare(m1, m2) != 0;
}
inline bool operator<(const M2T_O& m1, const M2T_O& m2) {
	return compare(m1, m2) < 0;
}

struct M2T {
	Int a, b, c;
	// This M2T-struct represents the matrix [a,b,c].
	M2T(Int _a = 0, Int _b = 0, Int _c = 0)
	: a(_a), b(_b), c(_c) {}
	Int det() { return a*c - b*b; }
};

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
	// in C: (1%3,0%3,-1%3,-2%3,-3%3,-4%3) == (1,0,-1,-2,0,-1)
	if(a % b < 0) res += b;
	assert(res >= 0);
	assert(res < b);
	return res;
}

template<typename T>
T Div(const T& a, const T& b) {
	T res = a / b;
	if(a % b < 0) res -= 1;
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


#include <chrono>
#include <string>

struct Timer {
	std::string name;
	std::chrono::steady_clock::time_point start;
	Timer(const std::string& _name) : name(_name) { start = std::chrono::steady_clock::now(); }
	~Timer() {
		auto d = std::chrono::steady_clock::now() - start;
		auto d_ms = std::chrono::duration_cast<std::chrono::milliseconds>(d).count();
		float d_s = d_ms / 1000.0f;
		std::cout << name << " took " << d_s << " secs" << std::endl;
	}
};


#endif
