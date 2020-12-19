#ifndef SJTU_MATRIX_HPP
#define SJTU_MATRIX_HPP

#include <cstddef>
#include <initializer_list>
#include <utility>
#include <iterator>
#include <string>
#include <algorithm>
#include <iostream>

using std::size_t;
using std::invalid_argument;

namespace sjtu
{
	
	template <class T>
	class Matrix
	{
    public:
		// your private member variables here.
		size_t row_max, column_max;
		T *storage;
	public:

        Matrix()
        {
            row_max = 0; column_max = 0;
            storage = nullptr;
        }
		
		Matrix(size_t n, size_t m, T _init = T())
		{
		    row_max = n; column_max = m;
			storage = new T[row_max * column_max];
			if (storage == nullptr) throw invalid_argument("something is wrong");
            for (int i = 0; i < row_max * column_max; ++i) {
                storage[i] = _init;
            }
		}
		
		explicit Matrix(std::pair<size_t, size_t> sz, T _init = T())
		{
		    row_max = sz.first; column_max = sz.second;
		    storage = new T[row_max * column_max];
		    if (storage == nullptr) throw invalid_argument("something is wrong");
            for (int i = 0; i < row_max * column_max ; ++i) {
                storage[i] = _init;
            }
		}

		Matrix(const Matrix &o)
		{
			row_max = o.row_max; column_max = o.column_max;
			storage = new T[row_max * column_max];
            if (storage == nullptr) throw invalid_argument("something is wrong");
            for (int i = 0; i < row_max * column_max; ++i) {
                storage[i] = o.storage[i];
            }
		}

		template <class U>
		Matrix(const Matrix<U> &o)
		{
			row_max = o.row_max; column_max = o.column_max;
			storage = new T[row_max * column_max];
            if (storage == nullptr) throw invalid_argument("something is wrong");
            for (int i = 0; i < row_max * column_max; ++i) {
                storage[i] = o.storage[i];
            }
		}

		Matrix &operator=(const Matrix &o)
		{
            if (storage == o.storage) return *this;
            delete [] storage;
			row_max = o.row_max; column_max = o.column_max;
			storage = new T[row_max * column_max];
            if (storage == nullptr) throw invalid_argument("something is wrong");
            for (int i = 0; i < row_max * column_max; ++i) {
                storage[i] = o.storage[i];
            }
            return *this;
		}
		
		template <class U>
		Matrix &operator=(const Matrix<U> &o)
		{
            delete [] storage;
			row_max = o.row_max; column_max = o.column_max;
			storage = new T[row_max * column_max];
            if (storage == nullptr) throw invalid_argument("something is wrong");
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
            if (storage == o.storage) return *this;
            delete [] storage;
		    row_max = o.row_max; column_max = o.column_max;
		    storage = o.storage; o.storage = nullptr;
            return *this;
		}

		~Matrix()
		{
            row_max = 0; column_max = 0;
            delete [] storage;
		}
		
		Matrix(std::initializer_list<std::initializer_list<T>> il)
		{
		    int i = 0;
		    auto it1 = il.begin(); auto it2 = (*it1).begin();
		    row_max = il.size(); column_max = (*it1).size();
		    storage = new T[row_max * column_max];
            if (storage == nullptr) throw invalid_argument("something is wrong");
            for (; it1 != il.end(); ++it1) {
                if ((*it1).size() != column_max ) {
                    delete [] storage;
                    throw invalid_argument("something is wrong");
                }
                it2 = (*it1).begin();
                for (; it2 != (*it1).end(); ++it2){
                    storage[i] = *it2; ++i;
                }
            }
		}
		
	public:
        size_t rowLength() const { return row_max; }
		
        size_t columnLength() const { return column_max; }

		void resize(size_t _n, size_t _m, T _init = T())
		{
            if (row_max * column_max == _n * _m)
            {
                row_max = _n; column_max = _m;
                return;
            }
			T *tmp_storage = new T[_n * _m];
            if (tmp_storage == nullptr) throw invalid_argument("something is wrong");
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
            delete [] storage;
            storage = tmp_storage;
		}

		void resize(std::pair<size_t, size_t> sz, T _init = T())
		{
            if (row_max * column_max == sz.first * sz.second)
            {
                row_max = sz.first; column_max = sz.second;
                return;
            }
		    T *tmp_storage = new T[sz.first * sz.second];
            if (tmp_storage == nullptr) throw invalid_argument("something is wrong");
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
		    row_max = sz.first; column_max = sz.second;
            delete [] storage;
		    storage = tmp_storage;
		}
		
		std::pair<size_t, size_t> size() const
		{
			std::pair<size_t, size_t> size_tmp;
			size_tmp = std::make_pair(row_max, column_max);
            return size_tmp;
		};
		
		void clear()
		{
		    if (storage != nullptr) {
                row_max = 0; column_max = 0;
                delete [] storage;
                storage = nullptr;
		    }
		}
		
	public:
		const T &operator()(size_t i, size_t j) const
		{
		    if (i < 0 || i >= row_max || j < 0 || j >= column_max) throw invalid_argument("something is wrong");
            return storage[i * column_max + j];
		}
		
		T &operator()(size_t i, size_t j)
		{
            if (i < 0 || i >= row_max || j < 0 || j >= column_max) throw invalid_argument("something is wrong");
            return storage[i * column_max + j];
		}
		
		Matrix<T> row(size_t i) const
		{
            if (i >= row_max || i < 0 ) throw invalid_argument("something is wrong");
			Matrix<T> row_matrix;
			row_matrix.row_max = 1; row_matrix.column_max = column_max;
			row_matrix.storage = new T[column_max];
            if (row_matrix.storage == nullptr) throw invalid_argument("something is wrong");
            for (int j = 0; j < column_max; ++j) {
                row_matrix.storage[j] = storage[i * column_max + j];
            }
            return row_matrix;
		}
		
		Matrix<T> column(size_t i) const
		{
            if (i >= column_max || i < 0) throw invalid_argument("something is wrong");
			Matrix<T> column_matrix;
			column_matrix.row_max = row_max; column_matrix.column_max = 1;
			column_matrix.storage = new T[row_max];
            if (column_matrix.storage == nullptr) throw invalid_argument("something is wrong");
            for (int j = 0; j < row_max; ++j) {
                column_matrix.storage[j] = storage[j * column_max + i];
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
		    Matrix<T> tmp_matrix;
		    tmp_matrix.row_max = row_max; tmp_matrix.column_max = column_max;
		    tmp_matrix.storage = new T[tmp_matrix.row_max * tmp_matrix.column_max];
            if (tmp_matrix.storage == nullptr) throw invalid_argument("something is wrong");
            for (int i = 0; i < tmp_matrix.row_max * tmp_matrix.column_max; ++i) {
                tmp_matrix.storage[i] = (-1) * storage[i];
            }
            return tmp_matrix;
		}

		template <class U>
		Matrix &operator+=(const Matrix<U> &o)
		{
            typedef decltype(T() + U()) V;
            if (o.row_max != row_max || o.column_max != column_max) throw invalid_argument("something is wrong");
            V *tmp_storage = new V[row_max * column_max];
            if (tmp_storage == nullptr) throw invalid_argument("something is wrong");
            for (int i = 0; i < row_max * column_max; ++i) {
                tmp_storage[i] = storage[i] + o.storage[i];
            }
            delete [] storage;
            storage = tmp_storage;
            return *this;
		}

		template <class U>
		Matrix &operator-=(const Matrix<U> &o)
		{
            typedef decltype(T() - U()) V;
            if (o.row_max != row_max || o.column_max != column_max) throw invalid_argument("something is wrong");
            V *tmp_storage = new V[row_max * column_max];
            if (tmp_storage == nullptr) throw invalid_argument("something is wrong");
            for (int i = 0; i < row_max * column_max; ++i) {
                tmp_storage[i] = storage[i] - o.storage[i];
            }
            delete [] storage;
            storage = tmp_storage;
            return *this;
		}
		
		template <class U>
		Matrix &operator*=(const U &x)
		{
            typedef decltype(U() * T()) V;
            V *tmp_storage = new V[row_max * column_max];
            if (tmp_storage == nullptr) throw invalid_argument("something is wrong");
            for (int i = 0; i < row_max * column_max; ++i) {
                tmp_storage[i] = x * storage[i];
            }
            delete [] storage;
            storage = tmp_storage;
            return *this;
		}
		
		Matrix tran() const
		{
			Matrix<T> tmp_matrix;
			size_t i, j;
			tmp_matrix.row_max = column_max; tmp_matrix.column_max = row_max;
			tmp_matrix.storage = new T[tmp_matrix.row_max * tmp_matrix.column_max];
            if (tmp_matrix.storage == nullptr) throw invalid_argument("something is wrong");
            for (i = 0; i < row_max; ++i) {
                for (j = 0; j < column_max; ++j) {
                    tmp_matrix.storage[j * tmp_matrix.column_max + i] = storage[i * column_max + j];
                }
            }
            return tmp_matrix;
		}

	public: // iterator
		class iterator
		{
		    friend iterator Matrix<T>::begin();
		    friend iterator Matrix<T>::end();
		    friend std::pair<iterator, iterator> Matrix<T>::subMatrix(std::pair<size_t, size_t> l, std::pair<size_t, size_t> r);
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
		    T *position;
		    size_t column_left, column_right, column_current, one_row, row_current, row_up, row_down;
		public:
			difference_type operator-(const iterator &o)
			{
                return position - o.position;
			}
			
			iterator &operator+=(difference_type offset)
			{
				position += offset;
                return *this;
			}
			
			iterator operator+(difference_type offset) const
			{
				iterator tmp_Itr;
				tmp_Itr.position = position + offset;
                return tmp_Itr;
			}
			
			iterator &operator-=(difference_type offset)
			{
				position -= offset;
                return *this;
			}
			
			iterator operator-(difference_type offset) const
			{
				iterator tmp;
				tmp.position = position - offset;
                return tmp;
			}
			
			iterator &operator++()
			{
                if (column_current < column_right)
                {
                    ++column_current;
                    ++position;
                }
                else
                {
                    if (row_current == row_down)
                    {
                        ++position;
                    }
                    else
                    {
                        column_current = column_left;
                        position += one_row - (column_right - column_left);
                        ++row_current;
                    }
                }
			}
			
			iterator operator++(int)
			{
			    iterator tmp = *this;
                if (column_current < column_right)
                {
                    ++column_current;
                    ++position;
                }
                else
                {
                    if (row_current == row_down)
                    {
                        ++position;
                    }
                    else
                    {
                        column_current = column_left;
                        position += one_row - (column_right - column_left);
                        ++row_current;
                    }
                }
                return tmp;
			}
			
			iterator &operator--()
			{
				if (column_current > column_left)
                {
				    --column_current;
				    --position;
                }
				else
                {
				    if (row_current == row_up)
                    {
				        --position;
                    }
				    else
                    {
                        column_current = column_right;
                        position += column_right - column_left - one_row;
                        --row_current;
                    }
                }
			}
			
			iterator operator--(int)
			{
			    iterator tmp = *this;
                if (column_current > column_left)
                    if (column_current > column_left)
                    {
                        --column_current;
                        --position;
                    }
                    else
                    {
                        if (row_current == row_up)
                        {
                            --position;
                        }
                        else
                        {
                            column_current = column_right;
                            position += column_right - column_left - one_row;
                            --row_current;
                        }
                    }
                return tmp;
			}
			
			reference operator*() const
			{
                return *position;
			}
			
			pointer operator->() const
			{
                return position;
			}
			
			bool operator==(const iterator &o) const
			{
                return position == o.position;
			}
			
			bool operator!=(const iterator &o) const
			{
                return position != o.position;
			}
		};
		
		iterator begin()

		{
			iterator tmp;
			tmp.column_left = 0; tmp.column_right = column_max - 1; tmp.one_row = column_max; tmp.column_current = 0;
            tmp.row_current = 0; tmp.row_up = 0; tmp.row_down = row_max - 1;
			tmp.position = storage;
            return tmp;
		}
		
		iterator end()
		{
			iterator tmp;
			tmp.column_left = 0; tmp.column_right = column_max - 1; tmp.one_row = column_max; tmp.column_current = column_max - 1;
            tmp.row_current = row_max - 1; tmp.row_up = 0;tmp.row_down = row_max - 1;
			tmp.position = storage + row_max * column_max;
            return tmp;
		}
		
		std::pair<iterator, iterator> subMatrix(std::pair<size_t, size_t> l, std::pair<size_t, size_t> r)
		{
			iterator it1, it2;
			std::pair<iterator, iterator> begin_end;
			it1.column_left = l.second; it1.column_right = r.second; it1.column_current = l.second; it1.one_row = column_max;
			it1.row_current = l.first; it1.row_up = l.first; it1.row_down = r.first;
			it2.column_left = l.second; it2.column_right = r.second; it2.column_current = r.second; it2.one_row = column_max;
			it2.row_current = r.first; it2.row_up = l.first; it2.row_down = r.first;
			it1.position = storage + l.first * column_max + l.second;
			it2.position = storage + r.first * column_max + r.second + 1;
			begin_end = std::make_pair(it1, it2);
            return begin_end;
        }

        void print_matrix()
        {
		    size_t i, j;
            for (i = 0; i < row_max; ++i) {
                for (j = 0; j < column_max; ++j) {
                    std::cout << storage[i * column_max + j] << ' ';
                }
                std::cout << std::endl;
            }
        }
	};
		
}

//
namespace sjtu
{
	template <class T, class U>
	auto operator*(const Matrix<T> &mat, const U &x)//->decltype(x * mat.storage[0])
	{
        typedef decltype(x * mat.storage[0]) V;
        Matrix<V> tmp_matrix;
		tmp_matrix.row_max = mat.row_max; tmp_matrix.column_max = mat.column_max;
		tmp_matrix.storage = new V[tmp_matrix.row_max * tmp_matrix.column_max];
        if (tmp_matrix.storage == nullptr) throw invalid_argument("something is wrong");
        for (int i = 0; i < tmp_matrix.row_max * tmp_matrix.column_max; ++i) {
            tmp_matrix.storage[i] = x * mat.storage[i];
        }
        return tmp_matrix;
	}
	
	template <class T, class U>
	auto operator*(const U &x, const Matrix<T> &mat)//->decltype(x * mat.storage[0])
	{
        typedef decltype(x * mat.storage[0]) V;
        Matrix<V> tmp_matrix;
		tmp_matrix.row_max = mat.row_max; tmp_matrix.column_max = mat.column_max;
		tmp_matrix.storage = new V[tmp_matrix.row_max * tmp_matrix.column_max];
        if (tmp_matrix.storage == nullptr) throw invalid_argument("something is wrong");
        for (int i = 0; i < tmp_matrix.row_max * tmp_matrix.column_max; ++i) {
            tmp_matrix.storage[i] = x * mat.storage[i];
        }
        return tmp_matrix;
	}
	
	template <class U, class V>
	auto operator*(const Matrix<U> &a, const Matrix<V> &b)//->decltype(a.storage[0] * b.storage[0])
	{
		if (a.column_max != b.row_max) throw invalid_argument("something is wrong");
        typedef decltype(a.storage[0] * b.storage[0]) W;
        Matrix<W> tmp_matrix;
        size_t i, j;
        tmp_matrix.row_max = a.row_max; tmp_matrix.column_max = b.column_max;
        tmp_matrix.storage = new W[tmp_matrix.row_max * tmp_matrix.column_max];
        if (tmp_matrix.storage == nullptr) throw invalid_argument("something is wrong");
        for (i = 0; i < tmp_matrix.row_max; ++i) {
            for (j = 0; j < tmp_matrix.column_max; ++j) {
                W &element = tmp_matrix.storage[i * tmp_matrix.column_max + j];
                element = 0;
                for (size_t k = 0; k < a.column_max; ++k) {
                    element += a.storage[i * a.column_max + k] * b.storage[k * b.column_max + j];
                }
            }
        }
        return tmp_matrix;
	}

	template <class U, class V>
	auto operator+(const Matrix<U> &a, const Matrix<V> &b)//->decltype(a.storage[0] + b.storage[0])
	{
	    if (a.row_max != b.row_max || a.column_max != b.column_max) throw invalid_argument("something is wrong");
        typedef decltype(a.storage[0] + b.storage[0]) W;
        Matrix<W> tmp_matrix;
		tmp_matrix.row_max = a.row_max; tmp_matrix.column_max = a.column_max;
		tmp_matrix.storage = new W[tmp_matrix.row_max * tmp_matrix.column_max];
        if (tmp_matrix.storage == nullptr) throw invalid_argument("something is wrong");
        for (int i = 0; i < tmp_matrix.row_max * tmp_matrix.column_max; ++i) {
            tmp_matrix.storage[i] = a.storage[i] + b.storage[i];
        }
        return tmp_matrix;
	}
	
	template <class U, class V>
	auto operator-(const Matrix<U> &a, const Matrix<V> &b)//->decltype(a.storage[0] - b.storage[0])
	{
        if (a.row_max != b.row_max || a.column_max != b.column_max) throw invalid_argument("something is wrong");
        typedef decltype(a.storage[0] - b.storage[0]) W;
        Matrix<W> tmp_matrix;
        tmp_matrix.row_max = a.row_max; tmp_matrix.column_max = a.column_max;
        tmp_matrix.storage = new W[tmp_matrix.row_max * tmp_matrix.column_max];
        if (tmp_matrix.storage == nullptr) throw invalid_argument("something is wrong");
        for (int i = 0; i < tmp_matrix.row_max * tmp_matrix.column_max; ++i) {
            tmp_matrix.storage[i] = a.storage[i] - b.storage[i];
        }
        return tmp_matrix;
    }

}

#endif //SJTU_MATRIX_HPP

