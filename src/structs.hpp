// Hermitian modular forms, https://github.com/albertz/diplom-thesis-math
// Copyright (c) 2013, Albert Zeyer, www.az2000.de
// This code is under the GPL v3 or later, see License.txt in the root directory of this project.

#ifndef __SAGE_HERMITIAN_STRUCTS_HPP__
#define __SAGE_HERMITIAN_STRUCTS_HPP__

#include <iostream>
#include <string.h>
#include <stdexcept>
#include <string>
#include <stdint.h>

#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)
#define DOMAIN_CHECK(x) { if(!(x)) throw std::domain_error(\
	(__FILE__ ":" TOSTRING(__LINE__) ": domain-check failed: " #x)); }
#define LOGIC_CHECK(x) { if(!(x)) throw std::logic_error(\
	(__FILE__ ":" TOSTRING(__LINE__) ": logic-check failed: " #x)); }

typedef int32_t Int;


template<typename T>
std::string int_to_bin(const T& n);

template<typename T>
T bin_to_int(const char* s);

template<>
inline
std::string int_to_bin<uint32_t>(const uint32_t& n) {
	std::string ret(4, '\0');
	// big endian
	ret[0] = (char)(uint8_t)(n & 0xff000000 >> 24);
	ret[1] = (char)(uint8_t)(n & 0x00ff0000 >> 16);
	ret[2] = (char)(uint8_t)(n & 0x0000ff00 >> 8);
	ret[3] = (char)(uint8_t)(n & 0x000000ff);
	return ret;
}

template<>
inline
uint32_t bin_to_int<uint32_t>(const char* s) {
	return
	(((uint32_t)(uint8_t)(s[0])) << 24) |
	(((uint32_t)(uint8_t)(s[1])) << 16) |
	(((uint32_t)(uint8_t)(s[2])) << 8) |
	(((uint32_t)(uint8_t)(s[3])));
}

template<>
inline
std::string int_to_bin<int32_t>(const int32_t& n) {
	return int_to_bin<uint32_t>((uint32_t) n);
}

template<>
inline
int32_t bin_to_int<int32_t>(const char* s) {
	return (int32_t) bin_to_int<uint32_t>(s);
}

inline
std::string string_to_bin(const std::string& s) {
	DOMAIN_CHECK(s.size() <= 0xffffffff);
	return int_to_bin<uint32_t>((uint32_t)s.size()) + s;
}

inline
std::string bin_to_string(const char* s, const char* end, uint32_t& len) {
	DOMAIN_CHECK(s + 4 < end);
	uint32_t strsize = bin_to_int<uint32_t>(s);
	DOMAIN_CHECK(s + strsize < end);
	len = 4 + strsize;
	return std::string(s + 4, strsize);
}


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

template<typename T>
T squareRootIntUpper(const T& y) {
	T x = squareRootInt(y);
	if(x*x == y) return x;
	return x + 1;
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
	// detD == -D * det
	// TODO(?): if D is fundamental, we always have 4|(D*D-D), thus we could just use -D * det.
	Int detD(const int D) const {
		DOMAIN_CHECK(D < 0);
		DOMAIN_CHECK(Mod(D*D - D, 4) == 0);
		// det = a*c - |b|^2
		// Re(b) = 1/2 b2
		// Re(b)^2 = 1/4 b2^2
		// Im(b) = -b1/sqrt{-D} + 1/2 \sqrt{-D} b2
		// Im(b)^2 = b1^2/(-D) - b1*b2 + 1/4 (-D) b2^2
		// -> (-D)*|b∫^2 = b1^2 - (-D)*b1*b2 + (D^2 - D)/4*b2^2
		return a*c*(-D) - b1*b1 + (-D)*b1*b2 - Div(D*D-D,4)*b2*b2;
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
	Int absBsquare(const int D) const {
		DOMAIN_CHECK(D < 0);
		DOMAIN_CHECK(Mod(D*D - D, 4) == 0);
		// Re(b) = b1 + D b2 1/2
		// Re(b)^2 = b1^2 + D b1 b2 + 1/4 D^2 b2^2
		// Im(b) = sqrt{-D} b2 1/2
		// Im(b)^2 = -D b2^2 1/4
		return b1*b1 + D*b1*b2 + Div(D*D-D, 4) * b2*b2;
	}
	// = Gauss_upper( M*|b| )
	Int absBMupper(const int D, const Int M) {
		DOMAIN_CHECK(M >= 0); // not implemented otherwise right now
		if(b2 == 0) return M * abs(b1); // fast path
		Int y = absBsquare(D) * M * M;
		return squareRootIntUpper(y);
	}
	Int det(const int D) const {
		return a*c - absBsquare(D);
	}
	Int gcd() const { return ::gcd(a, b1, b2, c); }
	
	std::string getState() const {
		std::string res = int_to_bin(a) + int_to_bin(b1) + int_to_bin(b2) + int_to_bin(c);
		LOGIC_CHECK(res.size() == 4 * sizeof(Int));
		return res;
	}
	void setState(const std::string& state) {
		DOMAIN_CHECK(state.size() == 4 * sizeof(Int));
		a  = bin_to_int<Int>(&state[sizeof(Int)*0]);
		b1 = bin_to_int<Int>(&state[sizeof(Int)*1]);
		b2 = bin_to_int<Int>(&state[sizeof(Int)*2]);
		c  = bin_to_int<Int>(&state[sizeof(Int)*3]);
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
inline std::ostream& operator<<(std::ostream& os, const M2T& m) {
	return os << "M2T(" << m.a << "," << m.b << "," << m.c << ")";
}


struct ElemOfCurlO {
	// We represent `b = b1 + b2 (D + \sqrt{D})/2`.
	Int b1, b2;
	ElemOfCurlO(Int _b1 = 0, Int _b2 = 0) : b1(_b1), b2(_b2) {}
	ElemOfCurlO conjugate(const int D) const {
		ElemOfCurlO res;
		res.b1 = b1 + b2 * D;
		res.b2 = -b2;
		return res;
	}
	ElemOfCurlO operator+(const ElemOfCurlO& other) const {
		ElemOfCurlO res;
		res.b1 = b1 + other.b1;
		res.b2 = b2 + other.b2;
		return res;
	}
	ElemOfCurlO mul(const ElemOfCurlO& other, const int D) const {
		DOMAIN_CHECK(Mod(D*D-D, 4) == 0);
		ElemOfCurlO res;
		res.b1 = b1 * other.b1 - b2 * other.b2 * Div(D*D - D, 4);
		res.b2 = b1 * other.b2 + b2 * other.b1 + D * b2 * other.b2;
		return res;
	}
};
inline std::ostream& operator<<(std::ostream& os, const ElemOfCurlO& m) {
	if(m.b1) os << m.b1;
	if(m.b1 && m.b2) os << "+";
	if(m.b2) os << m.b2 << "(D+sqrt(D))/2";
	if(!m.b1 && !m.b2) os << "0";
	return os;
}

struct ElemOfCurlOdual {
	// We represent `b = b1 / \sqrt{D} + b2 (1 + \sqrt{D})/2`.
	Int b1, b2;
	ElemOfCurlOdual() : b1(0), b2(0) {}
	ElemOfCurlOdual(Int _b1, Int _b2) : b1(_b1), b2(_b2) {}
	ElemOfCurlOdual conjugate(const int D) const {
		ElemOfCurlOdual res;
		res.b1 = -b1 - b2 * D;
		res.b2 = b2;
		return res;
	}
	static ElemOfCurlOdual fromInt(Int a, const int D) {
		ElemOfCurlOdual res;
		res.b1 = -a * D;
		res.b2 = a * 2;
		return res;
	}
	ElemOfCurlOdual operator+(const ElemOfCurlOdual& other) const {
		ElemOfCurlOdual res;
		res.b1 = b1 + other.b1;
		res.b2 = b2 + other.b2;
		return res;
	}
	ElemOfCurlOdual mul(const ElemOfCurlO& other, const int D) const {
		DOMAIN_CHECK(Mod(D - D*D, 4) == 0);
		ElemOfCurlOdual res;
		res.b1 = b1 * other.b1 + b2 * other.b2 * Div(D - D*D, 4);
		res.b2 = b1 * other.b2 + b2 * other.b1 + b2 * other.b2 * D;
		return res;
	}
	Int asInt(const int D) const {
		DOMAIN_CHECK(2 * b1 == -b2 * D); // must be in \RR
		DOMAIN_CHECK(Mod(b2, 2) == 0); // and in \ZZ
		return Div(b2, 2);
	}
};
inline std::ostream& operator<<(std::ostream& os, const ElemOfCurlOdual& m) {
	if(m.b1) os << m.b1 << "/sqrt(D)";
	if(m.b1 && m.b2) os << "+";
	if(m.b2) os << m.b2 << "(1+sqrt(D))/2";
	if(!m.b1 && !m.b2) os << "0";
	return os;
}




template<typename ElemType>
struct _M2_withD {
	ElemType a,b,c,d; // [[a,b],[c,d]]
	_M2_withD(ElemType _a = ElemType(), ElemType _b = ElemType(), ElemType _c = ElemType(), ElemType _d = ElemType())
	: a(_a), b(_b), c(_c), d(_d) {}
	_M2_withD conjugate_transpose(const int D) const {
		_M2_withD res;
		res.a = a.conjugate(D);
		res.b = c.conjugate(D);
		res.c = b.conjugate(D);
		res.d = d.conjugate(D);
		return res;
	}
	template<typename OtherMatrixType>
	_M2_withD mulMat(const OtherMatrixType& other, const int D) const {
		_M2_withD res;
		res.a = a.mul(other.a, D) + b.mul(other.c, D);
		res.b = a.mul(other.b, D) + b.mul(other.d, D);
		res.c = c.mul(other.a, D) + d.mul(other.c, D);
		res.d = c.mul(other.b, D) + d.mul(other.d, D);
		return res;
	}
	ElemType trace() const {
		return a + d;
	}
};
template<typename ElemType>
std::ostream& operator<<(std::ostream& os, const _M2_withD<ElemType>& m) {
	return os << "M2(" << m.a << "," << m.b << "," << m.c << "," << m.d << ")";
}

typedef _M2_withD<ElemOfCurlO> M2_O;
typedef _M2_withD<ElemOfCurlOdual> M2_Odual;

static inline M2_O M2_O_from_M2T_O(const M2T_O& m, const int D) {
	M2_O res;
	res.a = m.a;
	res.b = ElemOfCurlO(m.b1, m.b2);
	res.c = res.b.conjugate(D);
	res.d = m.c;
	return res;
}

static inline M2_Odual M2_Odual_from_M2T_Odual(const M2T_Odual& m, const int D) {
	M2_Odual res;
	res.a = ElemOfCurlOdual::fromInt(m.a, D);
	res.b = ElemOfCurlOdual(m.b1, m.b2);
	res.c = res.b.conjugate(D);
	res.d = ElemOfCurlOdual::fromInt(m.c, D);
	return res;
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
inline Int trace_old(M2T_O S, M2T_Odual T) {
	// = (S.a * T.a + S.b * \overline{T.b}) + (\overline{S.b} * T.b + S.b * T.b)
	// = S.a * T.a + 2 * Re(S.b * \overline{T.b}) + S.c * T.c
	// = S.a * T.a + S.c * T.c + (S.b1 * T.b2 - S.b2 * T.b1)
	return S.a * T.a + S.c * T.c + S.b1 * T.b2 - S.b2 * T.b1;
}

// calculates trace(S * T)
inline Int trace(M2T_O S, M2T_Odual T, int D) {
	auto _S = M2_O_from_M2T_O(S, D);
	auto _T = M2_Odual_from_M2T_Odual(T, D);
	auto m = _T.mulMat(_S, D); // T * S because other way around is not implemented
	auto t = m.trace().asInt(D);
	LOGIC_CHECK(t == trace_old(S, T)); // compare with old trace-calc
	return t;
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

#ifndef OLDGCC
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
#else
struct Timer { Timer(const char*) {} };
#endif


#endif
