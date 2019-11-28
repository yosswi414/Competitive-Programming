#pragma once

#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <algorithm>
#include <cassert>

struct BigInteger {
	static const int base = 10000;
	static const int digit = 4;
	std::vector<int> m_num;
	bool negative;
	BigInteger() :m_num(1, 0), negative(false) {}
	template<class T>
	BigInteger(const T&& n) {
		T k = n;
		static_assert(std::is_integral<T>::value, "error: floating point is invalid. cast to proper type.");
		int aux = 0;
		if (std::is_signed<T>::value && k < 0) {
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
		while (l[first_nonzero] == '0' && first_nonzero < (int)l.length() - 1)++first_nonzero;
		l = l.substr(first_nonzero);
		assert(isdigit(l[0]));
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
};

std::ostream& operator<<(std::ostream& os, const BigInteger& b) {
	std::string res;
	for (unsigned int i = 0; i + 1 < b.m_num.size(); ++i) {
		std::string s = std::to_string(b.m_num[i]);
		while (s.length() < b.digit)s = "0" + s;
		res = s + res;
	}
	res = std::to_string(*(b.m_num.end() - 1)) + res;
	if (b.negative)res = "-" + res;
	return os << res;
}