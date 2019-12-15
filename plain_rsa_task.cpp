/*
*     author: yosswi414
*     date : 15 Dec 2019
*     note : 
*/

#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <algorithm>
#include <random>
#include <cassert>

#pragma warning(push)
#pragma warning( disable : 4018 )
#pragma warning( disable : 4244 )
#pragma warning( disable : 4267 )
#pragma warning( disable : 26451 )
#pragma warning( disable : 26454 )

#ifndef _BIGINTEGER_4618059_
#define _BIGINTEGER_4618059_
struct BigInteger;

std::ostream& operator<<(std::ostream& os, const BigInteger& b);

template<class T> inline bool operator<(const BigInteger& lhs, const T& rhs);
template<class T> inline bool operator<(const T& lhs, const BigInteger& rhs);
inline bool operator<(const BigInteger& lhs, const BigInteger& rhs);
inline bool operator>(const BigInteger& lhs, const BigInteger& rhs);
inline bool operator==(const BigInteger& lhs, const BigInteger& rhs);
inline bool operator!=(const BigInteger& lhs, const BigInteger& rhs);
inline bool operator<=(const BigInteger& lhs, const BigInteger& rhs);
inline bool operator>=(const BigInteger& lhs, const BigInteger& rhs);

BigInteger operator+(const BigInteger& t1, const BigInteger& t2);
BigInteger operator-(const BigInteger& t1, const BigInteger& t2);
BigInteger operator*(const BigInteger& t1, const BigInteger& t2);
BigInteger operator/(const BigInteger& t1, const BigInteger& t2);
BigInteger operator%(const BigInteger& t1, const BigInteger& t2);

template<> constexpr const BigInteger& std::min(const BigInteger& lhs, const BigInteger& rhs);
template<> constexpr const BigInteger& std::max(const BigInteger& lhs, const BigInteger& rhs);

struct BigInteger {
	using ll = long long;
	static const int base = 100000000;
	static const int digit = 8;
	std::vector<ll> m_num;
	bool negative;
	BigInteger() :m_num(1, 0), negative(false) {}
	template<class T>
	BigInteger(const T&& n) {
		T k = n;
		static_assert(std::is_integral<T>::value, "error: floating point is invalid. cast to proper type.");
		if (k == 0) {
			*this = BigInteger();
			return;
		}
		int aux = 0;
		if (std::is_signed<T>::value&& k < 0) {
			negative = true;
			T nmax = 0 - std::numeric_limits<T>::max();
			if (k == std::numeric_limits<T>::min()) {
				k = std::numeric_limits<T>::max();
				if (std::numeric_limits<T>::min() == nmax - 1) aux = 1;
				else aux = 0;
			}
			else k *= -1;
		}
		else negative = false;
		m_num.clear();
		while (k > 0) {
			m_num.push_back(k % base);
			k /= base;
		}
		m_num[0] += aux;
	}
	template<>
	BigInteger(const std::string&& k) {
		std::string l = k;
		if (k[0] == '-')negative = true, l = l.substr(1);
		else negative = false;

		int first_nonzero = 0;
		while (l[first_nonzero] == '0' && first_nonzero < l.length() - 1)++first_nonzero;
		l = l.substr(first_nonzero);
		int pos = l.length() - 1;
		int x = 0;
		m_num.clear();
		do {
			std::string lp = l.substr(std::max(pos + 1 - digit, 0), std::min(pos + 1, digit));
			pos -= digit;
			//std::cerr << "cnstr l: " << l << ", lp: " << lp << ", pos: " << pos << std::endl;
			assert(lp.size() <= digit);
			x = stoi(lp);
			m_num.push_back(x);
		} while (pos >= 0);
	}
	// REMARK : move(n) is equivalent to static_cast<std::remove_reference<T>::type&&>(n)
	template<class T>
	BigInteger(const T& n) { *this = BigInteger(std::move(n)); }
	template<>
	BigInteger(const std::string& k) { *this = BigInteger(std::move(k)); }
	BigInteger(const char* k) { *this = BigInteger(std::string(k)); }

	inline size_t size() const { return m_num.size(); }
	std::string str() const {
		std::string res;
		for (unsigned int i = 0; i + 1 < size(); ++i) {
			std::string s = std::to_string(m_num[i]);
			while (s.length() < digit)s = "0" + s;
			res = s + res;
		}
		res = std::to_string(*(m_num.end() - 1)) + res;
		if (negative)res = "-" + res;
		return res;
	}
	std::string sstr() const {
		std::string res;
		for (unsigned int i = 0; i + 1 < size(); ++i) {
			std::string s = std::to_string(m_num[i]);
			while (s.length() < digit)s = "0" + s;
			res = s + res;
			res = " " + res; // for debug
		}
		res = std::to_string(*(m_num.end() - 1)) + res;
		if (negative)res = "-" + res;
		return res;
	}
	size_t length() const {
		return this->str().length();
	}
	void carry() {
		int psize = size();
		bool isref = true, issub = false;
		while (isref) {
			isref = false;
			for (int i = size() - 1; i >= 0; --i) {
				if (m_num[i] >= base) {
					isref = true;
					ll c = m_num[i] / base;
					m_num[i] %= base;
					if (i == size() - 1) {
						m_num.push_back(c);
						++psize;
					}
					else m_num[i + 1] += c;
				}
				else if (m_num[i] < 0) {
					isref = true;
					ll d = (-m_num[i]) / base + 1;
					if (i == size() - 1)negative = !negative, m_num.push_back(d), ++psize;
					else {
						m_num[i + 1] -= d;
						if (m_num[i + 1] == 0)issub = true;
					}
					m_num[i] += base * d;
				}
			}
		}
		if (psize > 1 && m_num[psize - 1] == 0) {
			issub = true;
			while (psize > 1 && m_num[psize - 1] == 0)--psize;
		}
		if (issub)m_num.resize(psize);
		if (psize == 1 && negative && m_num[0] == 0)negative = false;
	}


	// Cast
	explicit operator bool() const { return size() > 1 || m_num[0]; }
	explicit operator int() const { return (negative ? 1 : -1) * (int)(-m_num[0] - m_num[1] * base); }
	explicit operator long() const { return (negative ? 1 : -1) * (long)(-m_num[0] + (-m_num[1] - m_num[2] * base) * base); }
	explicit operator std::string() const { return this->str(); }

	BigInteger operator+() const { return BigInteger(*this); }
	BigInteger operator-() const {
		BigInteger bi(*this);
		if (size() > 1 || m_num[0] != 0) bi.negative = !bi.negative;
		return bi;
	}
	// Compound Assignment
	template<class T>
	BigInteger& operator+=(const T& rhs) {
		return *this = NaiveSum(*this, rhs);
	}
	template<class T>
	BigInteger& operator-=(const T& rhs) {
		return *this = NaiveSum(*this, -rhs);
	}
	template<class T>
	BigInteger& operator*=(const T& right) {
		BigInteger rhs(right);
		if (std::min(abs(*this), abs(BigInteger(rhs))) <= base)return *this = NaiveProduct(*this, rhs);
		else if (std::max(this->length(), rhs.length()))
			return *this = KaratsubaProduct(*this, rhs);
		else return *this = NaiveProduct(*this, rhs);
	}
	template<class T>
	BigInteger& operator/=(const T& rhs) {
		return *this = NaiveDivision(*this, rhs, false).first;
	}
	template<class T>
	BigInteger& operator%=(const T& rhs) {
		return *this = NaiveDivision(*this, rhs, true).second;
	}
	BigInteger& operator<<=(const size_t& rhs) {
		size_t k(rhs), su(26);
		while (k > 0) {
			for (auto& i : m_num)i <<= std::min(k, su);
			this->carry();
			k -= std::min(k, su);
		}
		return *this;
	}

	friend BigInteger abs(const BigInteger& bi) {
		BigInteger r(bi);
		r.negative = false;
		return r;
	}
	friend BigInteger pow10(size_t k) {
		if (!k)return 1;
		BigInteger res;
		std::vector<ll> n(k / digit, 0);
		k %= digit;
		n.push_back(1);
		for (; k > 0; --k)*(n.end() - 1) *= 10;
		res.m_num = n;
		return res;
	}
	friend BigInteger pow10(const BigInteger& bi, size_t k) {
		if (!k)return BigInteger(bi);
		BigInteger res;
		std::vector<ll> n(k / digit, 0);
		k %= digit;
		for (auto i : bi.m_num)n.push_back(i);
		res.m_num = n;
		res.negative = bi.negative;
		for (; k > 0; --k)res *= 10;
		return res;
	}
	friend BigInteger powb(size_t k) {
		if (!k)return 1;
		BigInteger res;
		std::vector<ll> n(k, 0);
		n.push_back(1);
		res.m_num = n;
		return res;
	}
	friend BigInteger powb(const BigInteger& bi, size_t k) {
		if (!k)return BigInteger(bi);
		BigInteger res;
		std::vector<ll> n(k, 0);
		for (auto i : bi.m_num)n.push_back(i);
		res.m_num = n;
		res.negative = bi.negative;
		return res;
	}
	friend BigInteger NaiveSum(const BigInteger& left, const BigInteger& right) {
		if (abs(left) < abs(right)) {
			return NaiveSum(right, left);
		}
		BigInteger lhs(left), rhs(right);
		int retsign = 0;
		if (lhs.negative) ++retsign, lhs.negative = !lhs.negative, rhs.negative = !rhs.negative;

		BigInteger p;

		p.m_num.resize(lhs.size());

		if (rhs.negative)
			for (int i = 0; i < rhs.size(); ++i)p.m_num[i] = lhs.m_num[i] - rhs.m_num[i];
		else
			for (int i = 0; i < rhs.size(); ++i)p.m_num[i] = lhs.m_num[i] + rhs.m_num[i];
		for (int i = rhs.size(); i < lhs.size(); ++i)p.m_num[i] = lhs.m_num[i];
		p.negative = retsign & 1;

		p.carry();
		return p;
	}
	friend BigInteger NaiveProduct(const BigInteger& left, const BigInteger& right) {
		BigInteger lhs(left), rhs(right);
		int retsign = 0;

		if (lhs.negative && rhs.negative) lhs.negative = !lhs.negative, rhs.negative = !rhs.negative;
		if (lhs.negative)++retsign, lhs.negative = !lhs.negative;
		if (rhs.negative)++retsign, rhs.negative = !rhs.negative;
		if (lhs.size() < rhs.size())std::swap(lhs, rhs);

		BigInteger res(0);
		for (int i = 0; i < lhs.size(); ++i) {
			BigInteger p;
			p.m_num = std::vector<ll>(rhs.size(), 0);
			for (int j = 0; j < rhs.size(); ++j) {
				p.m_num[j] += lhs.m_num[i] * rhs.m_num[j];
			}
			p.carry();
			res.m_num.resize(p.size() + i);
			for (int j = 0; j < p.size(); ++j)res.m_num[i + j] += p.m_num[j];
			res.carry();
		}
		if (retsign & 1)res.negative = !res.negative;
		return res;
	}
	friend BigInteger KaratsubaProduct(const BigInteger& left, const BigInteger& right, int bound = 64) {
		if (bound > 0 && std::min(left.size(), right.size()) <= bound)return NaiveProduct(left, right);
		BigInteger lhs(left), rhs(right);
		if (lhs.size() < rhs.size())std::swap(lhs, rhs);
		int retsign = 0;
		if (lhs.negative)++retsign;
		if (rhs.negative)++retsign;

		int n = std::max(lhs.size(), rhs.size());
		int rsize = n - (n / 2);

		BigInteger a0, a1, b0, b1;

		bool f[4] = { false };
		for (int i = rsize - 1; i >= 0; --i) {
			if (f[0] || (i + rsize < lhs.size() && lhs.m_num[i + rsize]>0)) {
				if (!f[0]) {
					a0.m_num.resize(i + 1);
					f[0] = true;
				}
				a0.m_num[i] = lhs.m_num[i + rsize];
			}
			if (f[1] || (i + rsize < rhs.size() && rhs.m_num[i + rsize]>0)) {
				if (!f[1]) {
					b0.m_num.resize(i + 1);
					f[1] = true;
				}
				b0.m_num[i] = rhs.m_num[i + rsize];
			}
			if (f[2] || lhs.m_num[i] > 0) {
				if (!f[2]) {
					a1.m_num.resize(i + 1);
					f[2] = true;
				}
				a1.m_num[i] = lhs.m_num[i];
			}
			if (f[3] || rhs.m_num[i] > 0) {
				if (!f[3]) {
					b1.m_num.resize(i + 1);
					f[3] = true;
				}
				b1.m_num[i] = rhs.m_num[i];
			}
		}
		BigInteger R2, R0;
		R2 = KaratsubaProduct(a0, b0, bound);
		R0 = KaratsubaProduct(a1, b1, bound);

		BigInteger R1 = R2 + R0 + -KaratsubaProduct(a0 + -a1, b0 + -b1, bound);
		std::cout << "K:\n";
		std::cout << R2 << "\n" << R1 << "\n" << R0 << "\n";

		BigInteger S1, S2;
		S1.m_num.resize(R1.size() + rsize);
		S2.m_num.resize(R2.size() + rsize * 2);
		S1.negative = R1.negative;
		S2.negative = R2.negative;
		for (int i = R1.size() - 1; i >= 0; --i)S1.m_num[i + rsize] = R1.m_num[i];
		for (int i = R2.size() - 1; i >= 0; --i)S2.m_num[i + rsize * 2] = R2.m_num[i];
		S2 += S1 + R0;
		if (retsign & 1)S2.negative = !S2.negative;
		return S2;
	}

	friend std::pair<BigInteger, BigInteger> NaiveDivision(const BigInteger& left, const BigInteger& right, bool isrem = true) {
		std::pair<BigInteger, BigInteger> res(0, 0);
		if (left < right) {
			res.second = left;
			return res;
		}

		BigInteger a(left), b(right);
		int retsign = 0;
		if (a.negative)++retsign, a.negative = false;
		if (b.negative)++retsign, b.negative = false;
		ll s = left.size() - 1, t = right.size() - 1;
		BigInteger& q = res.first;
		BigInteger& r = res.second;
		q.negative = r.negative = retsign % 2;
		ll d = base;
		if (s == t) {
			if (s == 0) {
				q = a.m_num[0] / b.m_num[0];
				r = a.m_num[0] % b.m_num[0];
				return res;
			}
			ll qu = std::min(a.m_num[s] / b.m_num[t], d - 1) + 1;
			ll ql = 1, qm;
			while (ql + 1 < qu) {
				qm = (qu + ql) / 2;
				r = a - qm * b;
				if (r.negative)qu = qm;
				else ql = qm;
			}
			q = ql;
			r = a - q * b;
			return res;
		}
		bool start = false, end = false;
		ll k = d / (b.m_num[t] + 1);
		a *= k;
		b *= k;
		s = a.size() - 1;
		t = b.size() - 1;
		ll u = s - t;
		ll qq;

		bool f = a < powb(b, u);
		while (!end) {
			f = a < powb(b, u);
			u -= f;
			if (u <= 0)end = true;
			if (!start) q.m_num.resize(u + 1), start = true;
			if (!f) {
				q.m_num[u] = 1;
				a -= powb(b, u);
			}
			else {
				qq = (a.m_num[s] * d + a.m_num[s - 1]) / b.m_num[t] + 1;
				a -= powb(b * qq, u);
				while (a.negative) qq -= 1, a += powb(b, u);
				assert(!a.negative);
				q.m_num[u] += qq;
			}
			s = a.size() - 1;
		}
		if (isrem) {
			if (k != 1 && a != 0)r = NaiveDivision(a, k, false).first;
			else r = a;
		}
		return res;
	}
};

template<class T>
inline bool operator<(const BigInteger& lhs, const T& rhs) {
	return lhs < BigInteger(rhs);
}
template<class T>
inline bool operator<(const T& lhs, const BigInteger& rhs) {
	return BigInteger(lhs) < rhs;
}
inline bool operator<(const BigInteger& lhs, const BigInteger& rhs) {
	if (lhs.negative != rhs.negative)return lhs.negative;
	if (lhs.negative) return -rhs < -lhs;
	if (lhs.size() != rhs.size())return lhs.size() < rhs.size();
	for (int i = lhs.size() - 1; i >= 0; --i) {
		if (lhs.m_num[i] != rhs.m_num[i])return lhs.m_num[i] < rhs.m_num[i];
	}
	return false;
}
inline bool operator>(const BigInteger& lhs, const BigInteger& rhs) { return rhs < lhs; }
inline bool operator==(const BigInteger& lhs, const BigInteger& rhs) { return lhs.negative == rhs.negative && lhs.m_num == rhs.m_num; }
inline bool operator!=(const BigInteger& lhs, const BigInteger& rhs) { return !(lhs == rhs); }
inline bool operator<=(const BigInteger& lhs, const BigInteger& rhs) { return lhs < rhs || lhs == rhs; }
inline bool operator>=(const BigInteger& lhs, const BigInteger& rhs) { return lhs > rhs || lhs == rhs; }

// Arithmetic Operation
BigInteger operator+(const BigInteger& t1, const BigInteger& t2) { return BigInteger(t1) += t2; }
BigInteger operator-(const BigInteger& t1, const BigInteger& t2) { return BigInteger(t1) -= t2; }
BigInteger operator*(const BigInteger& t1, const BigInteger& t2) { return BigInteger(t1) *= t2; }
BigInteger operator/(const BigInteger& t1, const BigInteger& t2) {
	if (t1 < t2)return 0;
	return BigInteger(t1) /= t2;
}
BigInteger operator%(const BigInteger& t1, const BigInteger& t2) {
	if (t2 == 2) {
		return t1.m_num[0] & 1;
	}
	if (t1 < t2) {
		return BigInteger(t1);
	}
	return BigInteger(t1) %= t2;
}

std::ostream& operator<<(std::ostream& os, const BigInteger& b) { return os << b.str(); }

// other function
template<> constexpr const BigInteger& std::min(const BigInteger& lhs, const BigInteger& rhs) { return rhs < lhs ? rhs : lhs; }
template<> constexpr const BigInteger& std::max(const BigInteger& lhs, const BigInteger& rhs) { return rhs > lhs ? rhs : lhs; }

template<class T>
std::pair<T, T> divrem(const T& lhs, const T& rhs) { return std::pair<T, T>(lhs / rhs, lhs % rhs); }
template<>
std::pair<BigInteger, BigInteger> divrem(const BigInteger& lhs, const BigInteger& rhs) { return NaiveDivision(lhs, rhs); }

#endif // #ifndef _BIGINTEGER_4618059_

template<class T>
T power(T x, T y, T p) {
	T res = 1;
	x = x % p;
	while (y) {
		if (y % 2 != 0)res = (res * x) % p;
		y /= 2;
		x = (x * x) % p;
	}
	return res;
}

template<class T>
struct PublicKey {
	T p, q;
	T e = (1 << 16) + 1;
};

template<class T>
T generatePublicKey(const PublicKey<T>& pk) {
	T phi = (pk.p - 1) * (pk.q - 1);

	T r0 = pk.e, r1 = phi;

	std::vector<T> k;
	while (r1 != 0) {
		assert(r0 > 0);
		assert(r1 > 0);
		std::pair<T, T> qr = divrem(r0, r1);
		k.push_back(qr.first);
		r0 = r1;
		r1 = qr.second;
	}
	assert(k.size() > 0);
	T x0, y0, u0, v0;
	T x1 = 0, y1 = 1, u1 = 1, v1 = -k[0];
	int h = (int)k.size();
	for (int i = 1; i < h; ++i) {
		u0 = u1;
		v0 = v1;
		x0 = x1;
		y0 = y1;
		x1 = u0;
		y1 = v0;
		u1 = x0 - u0 * k[i];
		v1 = y0 - v0 * k[i];
	}
	while (x1 < 0)x1 += phi;
	return x1;
}

BigInteger uniformBigInteger(const BigInteger& lower, const BigInteger& upper) {
	assert(lower <= upper);
	BigInteger nu = upper - lower;
	std::random_device rnd;
	std::mt19937_64 mt(rnd());
	std::uniform_int_distribution<> randbi(0, BigInteger::base - 1);
	std::uniform_int_distribution<> randmsb(0, nu.m_num[nu.size() - 1]);
	BigInteger rn;
	rn.m_num.resize(nu.size());
	for (int i = 0; i < rn.size() - 1; ++i)rn.m_num[i] = randbi(mt);
	rn.m_num[rn.size() - 1] = randmsb(mt);

	return rn + lower;
}

const int bound_mrpt = 1e8;
bool MillerRabbinPrimalityTest(const BigInteger& n) {
	if (n.negative) return false;
	if (n.size() == 1) {
		if (n.m_num[0] <= 1) return false;
		if (n.m_num[0] == 2) return true;
	}
	if (n % 2 == 0)return false;

	BigInteger d(n - 1), two(2);
	int s = 0;
	std::pair<BigInteger, BigInteger> dr;
	for (;;) {
		dr = divrem(d, two);
		if (dr.second != 0)break;
		d = dr.first;
		++s;
	}
	for (int k = 0; k < bound_mrpt; ++k) {
		BigInteger a(uniformBigInteger(1, n - 1));
		BigInteger ad(power(a, d, n));
		if (ad == 1) continue;
		for (int r = 1; r < s && ad != n - 1; ++r) {
			ad = (ad * ad) % n;
		}
		if (ad == n - 1)continue;
		else return false;
	}
	return true;
}


int main() {
	unsigned long long k = 0;
	using T = BigInteger;
	PublicKey<T> pk;
	pk.p = "12057981209962932719";
	pk.q = "17669034071322594473";

	std::string str;

	std::cout << "###### generate public key ######\n" << std::endl;
	do {
		std::cout << "please input prime number for p." << std::endl;
		std::cout << "enter a value for p ( 0 for default: p=" << pk.p << " )" << std::endl;
		std::cout << ">>> ";
		std::cin >> str;
		if (str == "0") break;
		if (!str.empty())pk.p = str;
	} while (!MillerRabbinPrimalityTest(pk.p));
	do {
		std::cout << "please input prime number for q." << std::endl;
		std::cout << "enter a value for q ( 0 for default: q=" << pk.q << " )" << std::endl;
		std::cout << ">>> ";
		std::cin >> str;
		if (str == "0") break;
		if (!str.empty())pk.q = str;
	} while (!MillerRabbinPrimalityTest(pk.q));

	std::cout << "please input certain value for e ( 0 for default: e=" << pk.e << " )" << std::endl;
	std::cout << ">>> ";
	std::cin >> str;
	if (str != "0" && !str.empty())pk.e = str;

	std::cout << "\n(p, q, e) = ( " << pk.p << ", " << pk.q << ", " << pk.e << " )" << std::endl;

	std::cout << "\n\n###### calculate secret key ######\n" << std::endl;
	T d = generatePublicKey(pk);
	std::cout << "secret key d = " << d << std::endl;

	BigInteger m = "31415926535897932384626433832795028841";
	m %= pk.p * pk.q;
	std::cout << "\n\n###### encrypt plain number ######\n" << std::endl;
	do {
		std::cout << "please input certain value to be encrypted ( lower than n = pq = " << pk.p * pk.q << " )" << std::endl;
		std::cout << "enter a value for m ( 0 for default: m=" << m << " ):" << std::endl;
		std::cout << ">>> ";
		std::cin >> str;
		if (str == "0") break;
		if (!str.empty())m = str;
	} while (m >= pk.p * pk.q);

	BigInteger c;
	std::cout << "\n[ encryption ]\nc = m ^ e mod n\n  = " << m << " ^ " << pk.e << " mod " << pk.p * pk.q << std::endl;
	std::cout << "  = " << (c = power(m, pk.e, pk.p * pk.q)) << std::endl;

	std::cout << "\n\n###### decrypt cryptgram ######\n" << std::endl;
	BigInteger mm;
	std::cout << "[ decryption ]\nm = c ^ d mod n\n  = " << c << " ^ " << d << " mod " << pk.p * pk.q << std::endl;
	std::cout << "  = " << (mm = power(c, d, pk.p * pk.q)) << std::endl;
}

#pragma warning(pop)
