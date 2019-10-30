/*
*   Date    : Oct 30 2019
*   Author  : yosswi414
*   Note    : under construction... it's mainly for the purpose of implementation of Discrete Fourier Transformation
*/

#pragma once

#include <iostream>

using R = double;
using ll = long long;
extern const double eps;
const R pi = 3.141592653589793;
const R e = 2.718281828459045;

// 複素数 (double)
class Complex {
	R _real = 0, _imaginary = 0;
public:
	Complex() :_real(0), _imaginary(0) {}
	Complex(R a) : _real(a), _imaginary(0) {}
	Complex(R a, R b) :_real(a), _imaginary(b) {}

	R& Re() { return _real; }
	R& Im() { return _imaginary; }
	const R Re() const { return _real; }
	const R Im() const { return _imaginary; }
	const R norm2() const { return _real * _real + _imaginary * _imaginary; }
	const R norm() const { return sqrt(norm2()); }
	const R arg() { return atan2(Im(), Re()); }
	const Complex conj() const { return Complex(Re(), -Im()); }

	const Complex operator+() const { return *this; }
	const Complex operator-() const { return Complex(-Re(), -Im()); }

	template<class T>
	Complex& operator=(const T a) {
		return *this = Complex(a, 0);
	}

	Complex& operator+=(const Complex w) {
		return *this = Complex(Re() + w.Re(), Im() + w.Im());
	}
	Complex& operator-=(const Complex w) {
		return *this = Complex(Re() - w.Re(), Im() - w.Im());
	}
	Complex& operator*=(const Complex w) {
		return *this = Complex(Re() * w.Re() - Im() * w.Im(), Re() * w.Im() + Im() * w.Re());
	}
	Complex& operator/=(const Complex w) {
		Complex n = *this * w.conj();
		double d = w.norm2();
		return *this = Complex(n.Re() / d, n.Im() / d);
	}

	bool operator==(const Complex w) const {
		return abs(Re() - w.Re()) < eps && abs(Im() - w.Im()) < eps;
	}
	bool operator!=(const Complex w) const {
		return !(*this == w);
	}

	template<class T, class U>
	friend Complex operator+(const T z, const U w) {
		return Complex(z) += Complex(w);
	}
	template<class T, class U>
	friend Complex operator-(const T z, const U w) {
		return Complex(z) -= Complex(w);
	}
	template<class T, class U>
	friend Complex operator*(const T z, const U w) {
		return Complex(z) *= Complex(w);
	}
	template<class T, class U>
	friend Complex operator/(const T z, const U w) {
		return Complex(z) /= Complex(w);
	}
};

const Complex uI = Complex(0, 1);

bool isnan(Complex z) {
	return isnan(z.Re()) || isnan(z.Im());
}

Complex polar(const R r, const R th) {
	if (abs(th) - 2 * pi > eps) {
		R thd = th / pi;
		ll thdc = (ll)thd;
		thd = (double)(thd - (thdc - thdc % 2)) * pi;
		if (th < 0)thd = -thd;
		return Complex(r * cos(thd), r * sin(thd));
	}
	return Complex(r * cos(th), r * sin(th));
}

Complex eix(const R x) {
	return polar(1, x);
}

Complex uniroot(const ll base, const ll m = 1) {
	return polar(1, (R)(m % base) / base * 2.0 * pi);
}

Complex zeta(const ll m, const ll base) {
	return uniroot(base, m);
}

std::ostream& operator<<(std::ostream& os, const Complex z) {
	if (isnan(z))return os << "nan[" << z.Re() << ", " << z.Im() << "]";
	if (abs(z.Re()) > eps) {
		os << (z.Re() > 0 ? "" : "-") << abs(z.Re());
		if (abs(z.Im()) > eps) {
			os << (z.Im() > 0 ? " + " : " - ");
			if (abs(abs(z.Im()) - 1) > eps) os << abs(z.Im());
			os << "i";
		}
	}
	else if (abs(z.Im()) > eps) {
		os << (z.Im() > 0 ? "" : "-");
		if (abs(abs(z.Im()) - 1) > eps) os << abs(z.Im());
		os << "i";
	}
	else os << 0;
	return os;
}

