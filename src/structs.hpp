#ifndef __SAGE_HERMITIAN_STRUCTS_HPP__
#define __SAGE_HERMITIAN_STRUCTS_HPP__

#include <iostream>
#include <string.h>
#include <stdexcept>
#include <string>

#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)
#define DOMAIN_CHECK(x) { if(!(x)) throw std::domain_error(\
	(__FILE__ ":" TOSTRING(__LINE__) ": domain-check failed: " #x)); }
#define LOGIC_CHECK(x) { if(!(x)) throw std::logic_error(\
	(__FILE__ ":" TOSTRING(__LINE__) ": logic-check failed: " #x)); }

typedef int Int;


template<typename T>
T Mod(const T& a, const T& b) {
	T res = a % b;
	// in C: (1%3,0%3,-1%3,-2%3,-3%3,-4%3) == (1,0,-1,-2,0,-1)
	if(a % b < 0) res += b;
	LOGIC_CHECK(res >= 0);
	LOGIC_CHECK(res < b);
	return res;
}

template<typename T>
T Div(const T& a, const T& b) {
	T res = a / b;
	if(a % b < 0) res -= 1;
	return res;
}

template<typename T>
T gcd(const T& a, const T& b) {
	if(a < 0) return gcd(-a, b);
	if(b < 0) return gcd(a, -b);
	if(a == 0) return b;
	if(a == 1) return 1;
	if(a > b) return gcd(b, a);
	return gcd(Mod(b,a), a);
}

template<typename T>
T gcd(const T& a, const T& b, const T& c) {
	return gcd(gcd(a,b), c);
}

template<typename T>
T gcd(const T& a, const T& b, const T& c, const T& d) {
	return gcd(a, b, gcd(c, d));
}

template<typename T>
T squareRootInt(const T& y) {
	DOMAIN_CHECK(y >= 0);
	DOMAIN_CHECK(y*y >= y); // this is an overflow check
	// This uses the Newton algorithm.
	T x = y;
	while(!(x*x <= y && (x+1)*(x+1) >= y)) {
		// y < x^2 ==> y/x < x ==> (x+y/x) < 2x ==> new_x < x
		x = (x + y/x) / 2;
	}
	return x;
}


struct M2T_Odual {
	// This represents always an element in Her_2(\cO^#) from our work.
	// We set `b = b1 / \sqrt{D} + b2 (1 + \sqrt{D})/2`, where
	// D is a negative integer, the fundamental discriminant of
	// the underlying imaginary quadratic number field.
	// This M2T_Odual-struct represents the matrix [a,b,c].
	Int a, b1, b2, c;
	M2T_Odual(Int _a = 0, Int _b1 = 0, Int _b2 = 0, Int _c = 0)
	: a(_a), b1(_b1), b2(_b2), c(_c) {}
	// det4D == -D * 4 * det
	// TODO(?): if D is fundamental, we always have 4|(D*D-D), thus we could just use -D * det.
	Int det4D(const int D) const {
		DOMAIN_CHECK(D < 0);
		// det = a*c - |b|^2
		// Re(b) = 1/2 b2
		// Re(b)^2 = 1/4 b2^2
		// Im(b) = -b1/sqrt{D} + 1/2 \sqrt{D} b2
		// Im(b)^2 = b1^2/(-D) - b1*b2 + 1/4 (-D) b2^2
		// -> 4*(-D)*|b∫^2 = 4*b1^2 - 4*(-D)*b1*b2 + (D^2 - D)*b2^2
		return a*c*4*(-D) - 4*b1*b1 + 4*(-D)*b1*b2 - (D*D - D)*b2*b2;
	}
};
inline std::ostream& operator<<(std::ostream& os, const M2T_Odual& m) {
	return os << "M2T_Odual(" << m.a << "," << m.b1 << "," << m.b2 << "," << m.c << ")";
}
inline int compare(const M2T_Odual& m1, const M2T_Odual& m2) {
	if(m1.a != m2.a) return (m1.a < m2.a) ? -1 : 1;
	if(m1.b1 != m2.b1) return (m1.b1 < m2.b1) ? -1 : 1;
	if(m1.b2 != m2.b2) return (m1.b2 < m2.b2) ? -1 : 1;
	if(m1.c != m2.c) return (m1.c < m2.c) ? -1 : 1;
	return 0;
}
inline bool operator==(const M2T_Odual& m1, const M2T_Odual& m2) {
	return compare(m1, m2) == 0;
}
inline bool operator!=(const M2T_Odual& m1, const M2T_Odual& m2) {
	return compare(m1, m2) != 0;
}
inline bool operator<(const M2T_Odual& m1, const M2T_Odual& m2) {
	return compare(m1, m2) < 0;
}

struct M2T_O {
	// This represents always an element in Her_2(\cO) from our work.
	// We set `b = b1 + b2 (D + \sqrt{D})/2`, where
	// D is a negative integer, the fundamental discriminant of
	// the underlying imaginary quadratic number field.
	// This M2T_O-struct represents the matrix [a,b,c].
	Int a, b1, b2, c;
	M2T_O(Int _a = 0, Int _b1 = 0, Int _b2 = 0, Int _c = 0)
	: a(_a), b1(_b1), b2(_b2), c(_c) {}
	// = |b|^2
	Int absBsquare(const int D) {
		DOMAIN_CHECK(D < 0);
		DOMAIN_CHECK(Mod(D*D - D, 4) == 0);
		// Re(b) = b1 + D b2 1/2
		// Re(b)^2 = b1^2 + D b2 + 1/4 D^2 b2^2
		// Im(b) = sqrt{-D} b2 1/2
		// Im(b)^2 = -D b2^2 1/4
		return b1*b1 + D*b2 + Div(D*D-D, 4) * b2*b2;
	}
	// = Gauss_upper( M*|b| )
	Int absBMupper(const int D, const Int M) {
		Int y = absBsquare(D) * M * M;
		Int x = squareRootInt(y);
		if(x*x == y) return x;
		return x + 1;
	}
	Int det(const int D) {
		return a*c - absBsquare(D);
	}
	Int gcd() const { return ::gcd(a, b1, b2, c); }
};
inline std::ostream& operator<<(std::ostream& os, const M2T_O& m) {
	return os << "M2T_O(" << m.a << "," << m.b1 << "," << m.b2 << "," << m.c << ")";
}


struct M2T {
	Int a, b, c;
	// This M2T-struct represents the matrix [a,b,c].
	M2T(Int _a = 0, Int _b = 0, Int _c = 0)
	: a(_a), b(_b), c(_c) {}
	Int det() { return a*c - b*b; }
};
inline std::ostream& operator<<(std::ostream& os, const M2T& m) {
	return os << "M2T(" << m.a << "," << m.b << "," << m.c << ")";
}


template<typename T>
struct Matrix2 {
	T a,b,c,d; // [[a,b],[c,d]]
	Matrix2(T _a = 0, T _b = 0, T _c = 0, T _d = 0)
	: a(_a), b(_b), c(_c), d(_d) {}
	T det() { return a*c - b*d; }
};


// calculates trace(S * T)
// is always an integer
inline Int trace(M2T_O S, M2T_Odual T) {
	// = (S.a * T.a + S.b * \overline{T.b}) + (\overline{S.b} * T.b + S.b * T.b)
	// = S.a * T.a + 2 * Re(S.b * \overline{T.b}) + S.c * T.c
	// = S.a * T.a + S.c * T.c + (S.b1 * T.b2 - S.b2 * T.b1)
	return S.a * T.a + S.c * T.c + S.b1 * T.b2 - S.b2 * T.b1;
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
