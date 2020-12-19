#ifndef SJTU_MATRIX_HPP
#define SJTU_MATRIX_HPP

#include <cstddef>
#include <initializer_list>
#include <utility>
#include <iterator>
#include <string>

using std::size_t;

namespace sjtu
{
	
	template <class T>
	class Matrix
	{
	private:
		// your private member variables here.
		size_t row_max, column_max;
		T *storage;
	public:
		Matrix() = default;
		
		Matrix(size_t n, size_t m, T _init = T())
		{
		    row_max = n; column_max = m;
			storage = new T[row_max * column_max];
            for (int i = 0; i < row_max * column_max; ++i) {
                storage[i] = _init;
            }
		}
		
		explicit Matrix(std::pair<size_t, size_t> sz, T _init = T())
		{
		    row_max = sz.first; column_max = sz.second;
		    storage = new T[row_max * column_max];
            for (int i = 0; i < row_max * column_max ; ++i) {
                storage[i] = _init;
            }
		}
		
		Matrix(const Matrix &o)
		{
			row_max = o.row_max; column_max = o.column_max;
			storage = new T[row_max * column_max];
            for (int i = 0; i < row_max * column_max; ++i) {
                storage[i] = o.storage[i];
            }
		}
		
		template <class U>
		Matrix(const Matrix<U> &o)
		{
			row_max = o.row_max; column_max = o.column_max;
			storage = new T[row_max * column_max];
            for (int i = 0; i < row_max * column_max; ++i) {
                storage[i] = o.storage[i];
            }
		}
		
		Matrix &operator=(const Matrix &o)
		{
			row_max = o.row_max; column_max = o.column_max;
            for (int i = 0; i < row_max * column_max; ++i) {
                storage[i] = o.storage[i];
            }
            return *this;
		}
		
		template <class U>
		Matrix &operator=(const Matrix<U> &o)
		{
			row_max = o.row_max; column_max = o.column_max;
            for (int i = 0; i < row_max * column_max; ++i) {
                storage[i] = o.storage[i];
            }
            return *this;
		}
		
		Matrix(Matrix &&o) noexcept
		{
			row_max = o.row_max; column_max = o.column_max;
			storage = o.storage; o.storage = nullptr;
		}
		
		Matrix &operator=(Matrix &&o) noexcept
		{
		    row_max = o.row_max; column_max = o.row_max;
		    storage = o.storage; o.storage = nullptr;
            return *this;
		}
		
		~Matrix()
		{
            delete [] storage;
		}
		
		Matrix(std::initializer_list<std::initializer_list<T>> il)
		{
			
		}
		
	public:
        size_t rowLength() const { return row_max; }
		
        size_t columnLength() const { return column_max; }
		//搞定
		void resize(size_t _n, size_t _m, T _init = T())
		{
			T *tmp_storage, *swap_storage;
			tmp_storage = new T[_n * _m];
			if (_n * _m <= row_max * column_max)
			{
                for (int i = 0; i < _n * _m; ++i) {
                    tmp_storage[i] = storage[i];
                }
			}
			else
            {
                for (int i = 0; i < row_max * column_max; ++i) {
                    tmp_storage[i] = storage[i];
                }
                for (int i = row_max * column_max; i < _n * _m; ++i) {
                    tmp_storage[i] = _init;
                }
            }
		    row_max = _n; column_max = _m;
            swap_storage = storage;
            storage = tmp_storage;
            delete [] swap_storage;
		}
		//搞定
		void resize(std::pair<size_t, size_t> sz, T _init = T())
		{
		    T *tmp_storage, *swap_storage;
		    tmp_storage = new T[sz.first * sz.second];
		    if (sz.first * sz.second <= row_max * column_max)
            {
                for (int i = 0; i < sz.first * sz.second; ++i) {
                    tmp_storage[i] = storage[i];
                }
            }
		    else
            {
                for (int i = 0; i < row_max * column_max; ++i) {
                    tmp_storage[i] = storage[i];
                }
                for (int i = row_max * column_max; i < sz.first * sz.second; ++i) {
                    tmp_storage[i] = _init;
                }
            }
		    swap_storage = storage;
		    storage = tmp_storage;
            delete [] swap_storage;
		}
		
		std::pair<size_t, size_t> size() const
		{
			std::pair<size_t, size_t> size_tmp;
			size_tmp.first = row_max; size_tmp.second = column_max;
            return size_tmp;
		};
		
		void clear()
		{
            delete [] storage;
		}
		
	public:
		const T &operator()(size_t i, size_t j) const
		{
            return storage[(i - 1) * column_max + j - 1];
		}
		
		T &operator()(size_t i, size_t j)
		{
            return storage[(i - 1) * column_max + j - 1];
		}
		
		Matrix<T> row(size_t i) const
		{
			Matrix<T> row_matrix;
			row_matrix.row_max = 1; row_matrix.column_max = column_max;
			row_matrix.storage = new T[column_max];
            for (int j = (i - 1) * column_max; j < i * column_max; ++j) {
                row_matrix.storage[j - (i - 1) * column_max] = storage[i];
            }
            return row_matrix;
		}
		
		Matrix<T> column(size_t i) const
		{
			Matrix<T> column_matrix;
			column_matrix.row_max = row_max; column_matrix.column_max = 1;
			column_matrix.storage = new T[row_max];
            for (int j = 0; j < row_max; ++j) {
                column_matrix.storage[j] = storage[i + j * column_max - 1];
            }
            return column_matrix;
		}
		
		
	public:
		template <class U>
		bool operator==(const Matrix<U> &o) const
		{
            if (row_max != o.row_max || column_max != o.column_max) return false;
            else
            {
                for (int i = 0; i < row_max * column_max; ++i) {
                    if (storage[i] != o.storage[i]) return false;
                }
            }
            return true;
		}
		
		template <class U>
		bool operator!=(const Matrix<U> &o) const
		{
            if (row_max != o.row_max || column_max != o.column_max) return true;
            else
            {
                for (int i = 0; i < row_max * column_max; ++i) {
                    if (storage[i] != o.storage[i]) return true;
                }
            }
            return false;
		}
		
		Matrix operator-() const
		{
            for (int i = 0; i < row_max * column_max; ++i) {
                storage[i] = (-1) * storage[i];
            }
            return *this;
		}

		///这里的异常抛出？两个数组类型不同的矩阵是否需要强制类型转化？
		template <class U>
		Matrix &operator+=(const Matrix<U> &o)
		{
            if (o.row_max != row_max || o.column_max != column_max) throw invalid_argument("fuck you");
            for (int i = 0; i < row_max * column_max; ++i) {
                storage[i] += o.storage[i];
            }
            return *this;
		}

		template <class U>
		Matrix &operator-=(const Matrix<U> &o)
		{
            for (int i = 0; i < row_max * column_max; ++i) {
                storage[i] -= o.storage[i];
            }
            return *this;
		}
		
		template <class U>
		Matrix &operator*=(const U &x)
		{
            for (int i = 0; i < row_max * column_max; ++i) {
                storage[i] *= x;
            }
		}
		
		Matrix tran() const
		{
			Matrix<T> tmp_matrix;
			size_t i, j;
			tmp_matrix.row_max = column_max; tmp_matrix.column_max = row_max;
			tmp_matrix.storage = new T[tmp_matrix.row_max * tmp_matrix.column_max];
            for (i = 1; i <= row_max; ++i) {
                for (j = 1; j <= column_max; ++j) {
                    tmp_matrix.storage[(j - 1) * tmp_matrix.column_max + i - 1] = storage[(i - 1) * column_max + j - 1];
                }
            }
            return tmp_matrix;
		}
		
	public: // iterator
		class iterator
		{
		public:
			using iterator_category = std::random_access_iterator_tag;
			using value_type        = T;
			using pointer           = T *;
			using reference         = T &;
			using size_type         = size_t;
			using difference_type   = std::ptrdiff_t;
			
			iterator() = default;
			
			iterator(const iterator &) = default;
			
			iterator &operator=(const iterator &) = default;
			
		private:

			
		public:
			difference_type operator-(const iterator &o)
			{
				
			}
			
			iterator &operator+=(difference_type offset)
			{
				
			}
			
			iterator operator+(difference_type offset) const
			{
				
			}
			
			iterator &operator-=(difference_type offset)
			{
				
			}
			
			iterator operator-(difference_type offset) const
			{
				
			}
			
			iterator &operator++()
			{
				
			}
			
			iterator operator++(int)
			{
				
			}
			
			iterator &operator--()
			{
				
			}
			
			iterator operator--(int)
			{
				
			}
			
			reference operator*() const
			{
				
			}
			
			pointer operator->() const
			{
				
			}
			
			bool operator==(const iterator &o) const
			{
				
			}
			
			bool operator!=(const iterator &o) const
			{
				
			}
		};
		
		iterator begin()
		{
			
		}
		
		iterator end()
		{
			
		}
		
		std::pair<iterator, iterator> subMatrix(std::pair<size_t, size_t> l, std::pair<size_t, size_t> r)
		{
			
        }
	};
		
}

//
namespace sjtu
{
	template <class T, class U>
	auto operator*(const Matrix<T> &mat, const U &x)
	{
		
	}
	
	template <class T, class U>
	auto operator*(const U &x, const Matrix<T> &mat)
	{
		
	}
	
	template <class U, class V>
	auto operator*(const Matrix<U> &a, const Matrix<V> &b)
	{
		
	}
	
	template <class U, class V>
	auto operator+(const Matrix<U> &a, const Matrix<V> &b)
	{
		
	}
	
	template <class U, class V>
	auto operator-(const Matrix<U> &a, const Matrix<V> &b)
	{
		
	}
	
}

#endif //SJTU_MATRIX_HPP

