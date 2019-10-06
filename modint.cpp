/*
*   Date    : Oct 7 2019
*   Author  : yosswi414
*   Note    : 
*/

#include <bits/stdc++.h>

using namespace std;

template<typename T, T Order>
class GaloisField {
	T a;
public:
	GaloisField(const T x = 0) :a(x) {}
	T& value() { return a; }
	const T&& value() const { return a; }
	operator uint_fast64_t() { return a; }
	GaloisField& operator=(const T x) {
		a = x % Order;
		return *this;
	}
	GaloisField operator+(const GaloisField b) const {
		return GaloisField(*this) += b;
	}
	GaloisField operator-(const GaloisField b) const {
		return GaloisField(*this) -= b;
	}
	GaloisField operator*(const GaloisField b) const {
		return GaloisField(*this) *= b;
	}
	GaloisField operator/(const GaloisField b) const {
		return GaloisField(*this) /= b;
	}
	GaloisField& operator+=(const GaloisField b) {
		a += b.a;
		if (a >= Order)a -= Order;
		return *this;
	}
	GaloisField& operator-=(const GaloisField b) {
		if (a < b.a)a += Order;
		a -= b.a;
		return *this;
	}
	GaloisField& operator*=(const GaloisField b) {
		a = a * (b.a % Order) % Order;
		return *this;
	}
	GaloisField& operator/=(const GaloisField b) {
		T pow = Order - 2;
		while (pow) {
			if (pow & 1) (*this) *= b.a;
			b.a *= b.a;
			pow <<= 1;
		}
		return *this;
	}
	GaloisField operator^(const GaloisField b) const {
		return a ^ b.a;
	}
};
