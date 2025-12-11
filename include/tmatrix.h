// ННГУ, ИИТММ, Курс "Алгоритмы и структуры данных"
//
// Copyright (c) Сысоев А.В.
//
//

#ifndef __TDynamicMatrix_H__
#define __TDynamicMatrix_H__

#include <iostream>
#include <algorithm>
#include <stdexcept>
#include <cassert>
#include <utility> // for std::move, std::swap
#include <limits>

using namespace std;

const int MAX_VECTOR_SIZE = 100000000;
const int MAX_MATRIX_SIZE = 10000;

// Динамический вектор - 
// шаблонный вектор на динамической памяти
template<typename T>
class TDynamicVector
{
protected:
  size_t sz;
  T* pMem;
public:
    TDynamicVector(int size = 1) : sz(0), pMem(nullptr) {
        if (size <= 0 || size > MAX_VECTOR_SIZE)
            throw out_of_range("Vector size should be in (0, MAX_VECTOR_SIZE]");
        sz = static_cast<size_t>(size);
        pMem = new(std::nothrow) T[sz]();
        if (pMem==nullptr)
            throw std::bad_alloc();
    }

    TDynamicVector(const TDynamicVector& v) : sz(v.sz), pMem(nullptr) {
        pMem = new T[sz];
        std::copy(v.pMem, v.pMem + sz, pMem);
    }

    TDynamicVector(TDynamicVector&& v) noexcept : sz(0), pMem(nullptr) {
        swap(*this, v);
    }

    ~TDynamicVector() {
        delete[] pMem;
        pMem = nullptr;
    }
    TDynamicVector& operator=(const TDynamicVector& v) {
        if (this != &v) {
            T* newMem = new T[v.sz];
            std::copy(v.pMem, v.pMem + v.sz, newMem);
            delete[] pMem;
            pMem = newMem;
            sz = v.sz;
        }
        return *this;
    }

    TDynamicVector& operator=(TDynamicVector&& v) noexcept {
        if (this != &v) {
            delete[] pMem;
            sz = 0;
            pMem = nullptr;
            swap(*this, v);
        }
        return *this;
    }

    size_t size() const noexcept { return sz; }
  // индексация
  T& operator[](size_t ind) 
  { 
      return pMem[ind];
  }
  const T& operator[](size_t ind) const 
  { 
      return pMem[ind];
  }
  
  // индексация с контролем
  T& at(size_t ind)
  {
      if (ind >= size())
          throw out_of_range("Index out of range");
      return pMem[ind];
  }
  const T& at(size_t ind) const
  {
      if (ind >= size())
          throw out_of_range("Index out of range");
      return pMem[ind];
  }

  // сравнение
  bool operator==(const TDynamicVector& v) const noexcept {
      if (sz != v.sz) return false;
      for (size_t i = 0; i < sz; ++i)
          if (!(pMem[i] == v.pMem[i])) return false;
      return true;
  }

  bool operator!=(const TDynamicVector& v) const noexcept {
      return !(*this == v); 
  }

  TDynamicVector operator+(T val) const {
      TDynamicVector res(*this);
      for (size_t i = 0; i < sz; ++i) res.pMem[i] += val;
      return res;
  }

  TDynamicVector operator-(T val) const {
      TDynamicVector res(*this);
      for (size_t i = 0; i < sz; ++i) res.pMem[i] -= val;
      return res;
  }

  TDynamicVector operator*(T val) const {
      TDynamicVector res(*this);
      for (size_t i = 0; i < sz; ++i) res.pMem[i] *= val;
      return res;
  }

  TDynamicVector operator+(const TDynamicVector& v) const {
      if (sz != v.sz) throw length_error("Vectors must have equal size");
      TDynamicVector res(*this);
      for (size_t i = 0; i < sz; ++i) res.pMem[i] += v.pMem[i];
      return res;
  }

  TDynamicVector operator-(const TDynamicVector& v) const {
      if (sz != v.sz) throw length_error("Vectors must have equal size");
      TDynamicVector res(*this);
      for (size_t i = 0; i < sz; ++i) res.pMem[i] -= v.pMem[i];
      return res;
  }

  T operator*(const TDynamicVector& v) const {
      if (sz != v.sz) throw length_error("Vectors must have equal size");
      T acc = T();
      for (size_t i = 0; i < sz; ++i) acc += pMem[i] * v.pMem[i];
      return acc;
  }

  friend void swap(TDynamicVector& lhs, TDynamicVector& rhs) noexcept
  {
      std::swap(lhs.sz, rhs.sz);
      std::swap(lhs.pMem, rhs.pMem);
  }

  // ввод/вывод
  friend istream& operator>>(istream& istr, TDynamicVector& v)
  {
      for (size_t i = 0; i < v.sz; i++)
         istr >> v.pMem[i]; // требуется оператор>> для типа T
      return istr;
  }
  friend ostream& operator<<(ostream& ostr, const TDynamicVector& v)
  {
    for (size_t i = 0; i < v.sz; i++)
      ostr << v.pMem[i] << ' '; // требуется оператор<< для типа T
    return ostr;
  }
};


// Динамическая матрица - 
// шаблонная матрица на динамической памяти
template<typename T>
class TDynamicMatrix : public TDynamicVector<TDynamicVector<T>>
{
private:
    size_t sz;
    TDynamicVector<TDynamicVector<T>> pMem;

public:
    TDynamicMatrix(int s = 1)
    {
        if (s <= 0 || s > MAX_MATRIX_SIZE)
            throw out_of_range("Matrix size should be in (0, MAX_MATRIX_SIZE]");

        sz = s;
        pMem = TDynamicVector<TDynamicVector<T>>(s);
        for (size_t i = 0; i < sz; ++i)
            pMem[i] = TDynamicVector<T>(s);
    }

    TDynamicMatrix(const TDynamicMatrix& m) : sz(m.sz), pMem(m.pMem) {}

    TDynamicMatrix(TDynamicMatrix&& m) noexcept : sz(0), pMem()
    {
        swap(sz, m.sz);
        swap(pMem, m.pMem);
    }

    TDynamicMatrix& operator=(const TDynamicMatrix& m) {
        if (this != &m) {
            sz = m.sz;
            pMem = m.pMem;
        }
        return *this;
    }

    TDynamicMatrix& operator=(TDynamicMatrix&& m) noexcept {
        if (this != &m) swap(*this, m);
        return *this;
    }

    size_t size() const noexcept { return sz; }

    TDynamicVector<T>& operator[](size_t ind) { return pMem[ind]; }
    const TDynamicVector<T>& operator[](size_t ind) const { return pMem[ind]; }

    TDynamicVector<T>& at(size_t ind) {
        if (ind >= sz) throw out_of_range("Index out of range");
        return pMem[ind];
    }

    const TDynamicVector<T>& at(size_t ind) const {
        if (ind >= sz) throw out_of_range("Index out of range");
        return pMem[ind];
    }

    bool operator==(const TDynamicMatrix& m) const noexcept {
        if (sz != m.sz) return false;
        for (size_t i = 0; i < sz; ++i)
            for (size_t j = 0; j < sz; ++j)
                if (!(pMem[i][j] == m.pMem[i][j]))
                    return false;
        return true;
    }

    bool operator!=(const TDynamicMatrix& m) const noexcept { return !(*this == m); }

    TDynamicMatrix operator+(const TDynamicMatrix& m) const {
        if (sz != m.sz) throw length_error("Matrices must have equal size");
        TDynamicMatrix res(*this);
        for (size_t i = 0; i < sz; ++i)
            res.pMem[i] = res.pMem[i] + m.pMem[i];
        return res;
    }

    TDynamicMatrix operator-(const TDynamicMatrix& m) const {
        if (sz != m.sz) throw length_error("Matrices must have equal size");
        TDynamicMatrix res(*this);
        for (size_t i = 0; i < sz; ++i)
            res.pMem[i] = res.pMem[i] - m.pMem[i];
        return res;
    }

    TDynamicMatrix operator*(const TDynamicMatrix& m) const {
        if (sz != m.sz) throw length_error("Matrices must have equal size");
        TDynamicMatrix res(static_cast<int>(sz));
        for (size_t i = 0; i < sz; ++i) {
            for (size_t j = 0; j < sz; ++j) {
                T sum = T();
                for (size_t k = 0; k < sz; ++k)
                    sum += pMem[i][k] * m.pMem[k][j];
                res.pMem[i][j] = sum;
            }
        }
        return res;
    }

    TDynamicVector<T> operator*(const TDynamicVector<T>& v) const {
        if (v.size() != sz) throw length_error("Vector size must equal matrix size");
        TDynamicVector<T> res(static_cast<int>(sz));
        for (size_t i = 0; i < sz; ++i)
            res[i] = pMem[i] * v;
        return res;
    }

    friend void swap(TDynamicMatrix& lhs, TDynamicMatrix& rhs) noexcept {
        std::swap(lhs.sz, rhs.sz);
        std::swap(lhs.pMem, rhs.pMem);
    }

    friend istream& operator>>(istream& istr, TDynamicMatrix& m) {
        for (size_t i = 0; i < m.sz; ++i)
            for (size_t j = 0; j < m.sz; ++j)
                istr >> m.pMem[i][j];
        return istr;
    }

    friend ostream& operator<<(ostream& ostr, const TDynamicMatrix& m) {
        for (size_t i = 0; i < m.sz; ++i) {
            for (size_t j = 0; j < m.sz; ++j)
                ostr << m.pMem[i][j] << ' ';
            if (i + 1 < m.sz) ostr << '\n';
        }
        return ostr;
    }
};

#endif
