/*
*   Date    : Oct 30 2019
*   Author  : yosswi414
*   Note    : I learned a lot about C++ by writing this... basic function/operations like det(A) will be added later
*/

#pragma once

#include <iostream>
#include <vector>
#include <utility>
#include <exception>

extern const double eps;

// 要素の型が T の行列
template<class T>
class Matrix {
	std::vector<std::vector<T>> mat;
public:
	Matrix() {
		mat = std::vector<std::vector<T>>(0);
	};
	Matrix(size_t r, size_t c, T x = static_cast<T>(0)) {
		mat = std::vector<std::vector<T>>(r, std::vector<T>(c, x));
	}
	Matrix(std::pair<size_t, size_t> S, T x = static_cast<T>(0)) {
		mat = std::vector<std::vector<T>>(S.first, std::vector<T>(S.second, x));
	}
	Matrix(std::initializer_list<std::vector<T>> init) {
		size_t newRow = init.size(), newCol = newRow > 0 ? init.begin()->size() : 0;
		*this = Matrix(newRow, newCol);
		size_t r = 0;
		for (const auto v : init) {
			if (newCol != v.size()) throw std::invalid_argument("Matrix Operation Error: Cannot interpret non-matrix argument");
			mat[r++] = v;
		}
	}
	operator std::vector<std::vector<T>>() { return mat; }
	size_t row() const {
		return mat.size();
	}
	size_t column() const {
		return row() > 0 ? mat[0].size() : 0;
	}
	std::vector<T>& operator [] (const size_t i) { return mat[i]; }
	const std::vector<T>& operator [] (const size_t i) const { return mat[i]; }
	std::pair<size_t, size_t> size() const {
		return { row(), column() };
	}
	Matrix operator +(const Matrix B) const {
		const Matrix& A = *this;
		if (A.size() != B.size()) throw std::invalid_argument("Matrix Operation Error: Cannot operate two different size matrice (operator: +)");
		Matrix X(A.size(), static_cast<T>(0));
		for (size_t i = 0; i < row(); ++i) for (size_t j = 0; j < column(); ++j) X[i][j] = A[i][j] + B[i][j];
		return X;
	}
	Matrix& operator +=(const Matrix A) {
		return *this = (*this) + A;
	}
	Matrix operator *(const Matrix B) const {
		const Matrix& A = *this;
		if (A.column() != B.row()) throw std::invalid_argument("Matrix Operation Error: Column of left and row of right matrix must be equal in terms of multiplication");
		Matrix X(A.row(), B.column(), static_cast<T>(0));
		for (size_t i = 0; i < A.row(); ++i)for (size_t j = 0; j < B.column(); ++j) {
			for (size_t k = 0; k < A.column(); ++k) X[i][j] += A[i][k] * B[k][j];
		}
		return X;
	}
	Matrix& operator *=(const Matrix A) {
		return *this = (*this) * A;
	}
	Matrix operator *(const T k) const {
		Matrix X = *this;
		for (size_t i = 0; i < X.row(); ++i)for (size_t j = 0; j < X.column(); ++j)X[i][j] *= k;
		return X;
	}
	Matrix& operator *=(const T k) {
		Matrix& X = *this;
		for (size_t i = 0; i < X.row(); ++i)for (size_t j = 0; j < X.column(); ++j)X[i][j] *= k;
		return X;
	}
	Matrix operator ~() const {
		const Matrix& A = *this;
		Matrix X(column(), row());
		for (size_t j = 0; j < column(); ++j)for (size_t i = 0; i < row(); ++i)X[j][i] = A[i][j];
		return X;
	}
	Matrix operator ^(const Matrix& B) const {
		const Matrix& A = *this;
		if (A.size() != B.size()) throw std::invalid_argument("Matrix Operation Error: Cannot operate two different size matrice (operator: ^)");
		Matrix X = A;
		for (size_t i = 0; i < X.row(); ++i)for (size_t j = 0; j < X.column(); ++j)X[i][j] = A[i][j] ^ B[i][j];
		return X;
	}
	
};

template<class T>
std::ostream& operator<< (std::ostream& os, const Matrix<T> A) {
	if (A.row() < 1)return os << std::endl;

	if (A.row() == 1) {
		os << "( ";
		for (size_t i = 0; i < A.column(); ++i) {
			try {
				os << (fabs(A[0][i]) < eps ? 0.0 : A[0][i]);
			}
			catch (std::exception & excp) {
				os << A[0][i];
			}
			os << " ";
		}
		os << ")\n";
		return os;
	}

	else for (size_t i = 0; i < A.row(); ++i) {
		os << (i ? (i < A.row() - 1 ? "| " : "\\ ") : "/ ");
		for (size_t j = 0; j < A.column(); ++j) {
			try {
				os << (fabs(A[i][j]) < eps ? 0.0 : A[i][j]);
			}
			catch (std::exception & excp) {
				os << A[i][j];
			}
			os << " ";
		}
		os << (i ? (i < A.row() - 1 ? "|" : "/") : "\\") << std::endl;
	}
	return os;
}
