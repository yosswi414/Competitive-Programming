/*
*   Date    : Oct 7 2019
*   Author  : yosswi414
*   Note    : I learned a lot about C++ by writing this... basic function/operations like det(A) will be added later
*/

#include <bits/stdc++.h>

using namespace std;

template<typename T>
class Matrix {
	vector<vector<T>> mat;
public:
	Matrix() {
		mat = vector<vector<T>>(0);
	};
	Matrix(size_t r, size_t c, T x = static_cast<T>(0)) {
		mat = vector<vector<T>>(r, vector<T>(c, x));
	}
	Matrix(pair<size_t, size_t> S, T x = static_cast<T>(0)) {
		mat = vector<vector<T>>(S.first, vector<T>(S.second, x));
	}
	Matrix(initializer_list<vector<T>> init) {
		size_t newRow = init.size(), newCol = newRow > 0 ? init.begin()->size() : 0;
		*this = Matrix(newRow, newCol);
		size_t r = 0;
		for (const auto v : init) {
			if (newCol != v.size()) throw invalid_argument("Matrix Operation Error: Cannot interpret non-matrix argument");
			mat[r++] = v;
		}
	}
	operator vector<vector<T>>() { return mat; }
	size_t row() const {
		return mat.size();
	}
	size_t column() const {
		return row() > 0 ? mat[0].size() : 0;
	}
	vector<T>& operator [] (const size_t i) { return mat[i]; }
	const vector<T>& operator [] (const size_t i) const { return mat[i]; }
	pair<size_t, size_t> size() const {
		return { row(), column() };
	}
	Matrix operator +(const Matrix B) const {
		const Matrix& A = *this;
		if (A.size() != B.size()) throw invalid_argument("Matrix Operation Error: Cannot operate two different size matrice (operator: +)");
		Matrix X(A.size(), static_cast<T>(0));
		for (size_t i = 0; i < row(); ++i) for (size_t j = 0; j < column(); ++j) X[i][j] = A[i][j] + B[i][j];
		return X;
	}
	Matrix& operator +=(const Matrix A) {
		return *this = (*this) + A;
	}
	Matrix operator *(const Matrix B) const {
		const Matrix& A = *this;
		if (A.column() != B.row()) throw invalid_argument("Matrix Operation Error: Column of left and row of right matrix must be equal in terms of multiplication");
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
		if (A.size() != B.size()) throw invalid_argument("Matrix Operation Error: Cannot operate two different size matrice (operator: ^)");
		Matrix X = A;
		for (size_t i = 0; i < X.row(); ++i)for (size_t j = 0; j < X.column(); ++j)X[i][j] = A[i][j] ^ B[i][j];
		return X;
	}
  
  // to be implemented ...
  T det(const Matrix& A) const;
  Matrix inv() const;
  Matrix power(unsigned int t) const;
  // this work is so fun i'll add some more function
};
