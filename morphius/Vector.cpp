#include <iostream>
#include <assert.h>
#include <math.h>

#define MAX(x,y) ((x>y)?x:y)
#define MIN(x,y) ((x<y)?x:y)

template <typename T> class Vector {
  private:
    size_t size;
    T *vec, *ptr;
    void resize(size_t new_size) {
      size_t index = this->ptr-this->vec;
      this->size = new_size;
      T *new_vec = (T*) malloc(sizeof(T)*this->size);
      for (int i = 0; i < this->size; i++) *(new_vec+i) = *(this->vec+i);
      free(this->vec);
      this->vec = new_vec;
      this->ptr = this->vec+index;
    }
    void fill(T value) {
      for (int i = 0; i < this->size; i++) {
        *(this->vec+i) = value;
      }
    }
  public:
    void allocate(size_t N) {
      this->size = N;
      this->vec = (T*) malloc(sizeof(T)*N);
      this->ptr = this->vec;
      this->fill(0);
    }
    void push(T value) {
      if (this->ptr-this->vec == this->size) {
        this->resize(this->size+1);
      }
      *(ptr++) = value;
    }
    // index pass by reference for assignment
    T& operator[](const size_t index) {
      assert(index >= 0 && index < this->size);
      return *(this->vec+index);
    }
    // dot product between two vectors (of equal size)
    friend T operator*(Vector<T> &v1, Vector<T> &v2) {
      assert(v1.size == v2.size);
      T dot = 0;
      for (int i = 0; i < v1.size; i++) dot += v1[i]*v2[i];
      return dot;
    }
    // multiply vector by scalar
    friend Vector<T>& operator*(Vector<T> &v1, T scalar) {
      static Vector<T> tmp; tmp.allocate(v1.size);
      for (int x = 0; x < v1.size; x++) *(tmp.vec+x) = *(v1.vec+x)*scalar;
      return tmp;
    }
    // add two vectors
    friend Vector<T>& operator+(Vector<T> &v1, Vector<T> &v2) {
      assert(v1.size == v2.size);
      static Vector<T> tmp; tmp.allocate(v1.size);
      for (int i = 0; i < v1.size; i++) *(tmp.vec+i) = *(v1.vec+i) + *(v2.vec+i);
      return tmp;
    }
    // subtract two vectors
    friend Vector<T>& operator-(Vector<T> &v1, Vector<T> &v2) {
      assert(v1.size == v2.size);
      static Vector<T> tmp; tmp.allocate(v1.size);
      for (int i = 0; i < v1.size; i++) *(tmp.vec+i) = *(v1.vec+i) - *(v2.vec+i);
      return tmp;
    }
    // print to std::out
    friend std::ostream& operator<<(std::ostream& os, Vector<T> &o) {
      for (int i = 0; i < o.size; i++) os << '[' << o.vec[i] << ']';
      return os;
    }
    // equality
    friend bool operator==(const Vector<T> &v1, const Vector<T> &v2) {
      for (int i = 0; i < MIN(v1.size, v2.size); i++) {
        if (*(v1.vec+i) != *(v2.vec+i)) { return false; }
      }
      return (v1.size == v2.size);
    }
    friend bool operator!=(const Vector<T> &v1, const Vector<T> &v2) {
      return !(v1 == v2);
    }
		/*Vector<T>& operator=(Vector<T>& v1) {
			// assign one vector to another
			if (this == &v1) return *this;
			if (this->size != v1.size) {
				this->allocate(v1.size);
			}
		}*/
};

template <typename T> Vector<T> createVec(size_t N) {
  assert(std::is_arithmetic<T>::value);
  Vector<T> tmp;
  tmp.allocate(N);
  return tmp;
}

template <typename T> class Matrix {
  private:
    size_t rows, cols;
    Vector<T> *mat;
    Vector<T> columnVec(size_t C) {
      assert(C >= 0 && C < this->cols);
      Vector<T> tmp = createVec<T>(this->rows);
      for (int R = 0; R < this->rows; R++) tmp.push(this->get(R, C));
      return tmp;
    }
    // submatrix: don't include rows and column of index.
    Matrix<T>& submatrix(size_t index) {
      static Matrix<T> sub; sub.allocate(this->rows-1, this->cols-1);
      for (int r = 0; r < this->rows; r++) {
        if (r == index) continue;
        for (int c = 0; c < this->cols; c++) {
          if (c == index) continue;
          sub.insert(r-(r>index?1:0), c-(c>index?1:0), this->get(r,c));
        }
      }
      return sub;
    }
  public:
    bool isSquare() {
      return (this->cols == this->rows);
    }
    size_t getRows() { return this->rows; }
    size_t getCols() { return this->cols; }
    void allocate(size_t R, size_t C) {
      this->rows = R; this->cols = C;
      this->mat = (Vector<T>*) malloc(sizeof(Vector<T>)*this->rows);
      for (int i = 0; i < this->rows; i++) {
        *(this->mat+i) = createVec<T>(this->cols);
      }
    }
    // insert to position in the matrix
    void insert(size_t R, size_t C, T value) {
      assert(R >= 0 && R < this->rows && C >= 0 && C < this->cols);
      this->mat[R][C] = value;
    }
    // load 1-d array using (r*rows+cols = size)
    void loadArray(T* arr, size_t size) {
      for (int r = 0; r < this->rows; r++) {
        for (int c = 0; c < this->cols; c++) {
          size_t index = r*this->rows+c;
          if (index < size) {
            this->mat[r][c] = *(arr+index);
          } else { break; }
        }
      }
    }
    T get(size_t R, size_t C) {
      assert(R >= 0 && R < this->rows && C >= 0 && C < this->cols);
      return this->mat[R][C];
    }
    // switching the row and columns to produce the transpose
    Matrix<T>& transpose() {
      static Matrix<T> tmp; tmp.allocate(this->cols, this->rows);
      for (int i = 0; i < this->cols; i++) {
        tmp.mat[i] = this->columnVec(i);
      }
      return tmp;
    }
		/* Matrix must be square
		 * Convert to reduced row echelon form
		 * Determinant = Product of diagonal values
	  */
    T determinant() {
      assert(this->isSquare());
      T det = 0;
      if (this->rows == 2) {
        return (this->mat[0][0]*this->mat[1][1])-(this->mat[0][1]*this->mat[1][0]);
      } else {
        for (int i = 0; i < this->rows; i++) {
          det += this->submatrix(i).determinant();
        }
      }
      return det;
    }
		void swapRows(size_t i, size_t k) {
			Vector<T> tmp = this->mat[k];
			this->mat[k] = this->mat[i];
			this->mat[i] = tmp;
		}
    // passing in lambda function
    template <typename L> void lambdaFill(const L& lfunc) {
      for (int r = 0; r < this->rows; r++) {
        for (int c = 0; c < this->cols; c++) {
          this->mat[r][c] = lfunc(r+c);
        }
      }
    }
    Vector<T>& operator[](const size_t index) {
      assert(index >= 0 && index < this->rows);
      return *(this->mat+index);
    }
    // MATMUL (m1 rows * m2 cols)
    friend Matrix<T>& operator*(Matrix<T> &m1, Matrix<T> &m2) {
      assert(m1.cols == m2.rows);
      static Matrix<T> product; product.allocate(m1.rows, m2.cols);
      for (int r = 0; r < m1.rows; r++) {
        for (int c = 0; c < m2.cols; c++) {
          product[r][c] = m1[r] * m2[c];
        }
      }
      return product;
    }
    friend bool operator==(Matrix<T> &m1, Matrix<T> &m2) {
      for (int r = 0; r < MIN(m1.rows, m2.rows); r++) {
        if (m1[r] != m2[r]) return false;
      }
      return (m1.rows == m2.rows && m1.cols == m2.cols);
    }
    friend std::ostream& operator<<(std::ostream& os, Matrix<T> &m) {
      for (int i = 0; i < m.rows; i++) {
        os << *(m.mat+i) << (i!=m.rows-1?"\n":"");
      }
      return os;
    }

};

template <typename T> Matrix<T> createMatrix(size_t R, size_t C) {
  assert(std::is_arithmetic<T>::value);
  Matrix<T> tmp;
  tmp.allocate(R, C);
  return tmp;
}

/* Identity Matrix
(3x3 I)  [1, 0, 0]
         [0, 1, 0]
         [0, 0, 1]
 */
template <typename T> Matrix<T> createIdentity(size_t N) {
  Matrix<T> I = createMatrix<T>(N, N);
  for (int i = 0; i < N; i++) {
    I[i][i] = 1;
  }
  return I;
}

/* How does inverting a matrix work?
  - GOAL: Find the matrix N such that N*M = I (where I is the identity matrix)
  - Reduce M into Reduced Row Echelon Form (RREF --> M == I)
*/
template <typename T> Matrix<T> invert(Matrix<T>& m) {
  assert(m.isSquare() && m.determinant() != 0);
  Matrix<T> I = createIdentity<T>(m.getRows());
	// O(n^2)...
	for (size_t k = 0; k < m.getRows(); ++k) {
		bool zeroed = false;
		// swap rows if (k,k) == 0
		if ((zeroed = (m[k][k] == 0))) {
			for (size_t i = 0; i < m.getRows(); ++i) {
				if (m[i][k] != 0) { m.swapRows(i,k); I.swapRows(i,k); }
			}
			std::cout << zeroed << "\n";
		}

		// normalize using index at (k,k)
		I[k] = I[k]*(1/m[k][k]);
		m[k] = m[k]*(1/m[k][k]);
		for (size_t i = 0; i < m.getRows(); ++i) {
			// create reduced row echelon form
			if (k == i) continue;
			double s = m[i][k];
			m[i] = m[i]-m[k]*s;
			I[i] = I[i]-I[k]*s;
		}
	}

	return I;
}

int main(void) {
  auto fill = [](const auto x) { return x+1; };
	/*float f[] = {3,0,2,2,0,-2,0,1,1};
	Matrix<float> m = createMatrix<float>(3, 3);
	m.loadArray((float*) f, 9);*/
	Matrix<float> m = createMatrix<float>(4,4);
	float f[] = {4,0,0,0,0,0,2,0,0,1,2,0,1,0,0,1};
	m.loadArray((float*) f, 16);
	std::cout << m << "\n";
	std::cout << m.determinant() << "\n";


	/*
	Matrix<float> inverse = invert(m);
	std::cout << inverse << "\n";*/
}
