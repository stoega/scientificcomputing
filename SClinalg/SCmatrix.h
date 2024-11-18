#pragma once
#include <random>
#include <iostream>
#include <complex>

#include "SCvector.h"

namespace SC
{

    template <typename T = double>
    class LinearOperator
    {
    protected:
        int height;
        int width;

        LinearOperator() : height(0), width(0) {}
        LinearOperator(int h) : height(h), width(h) {}
        LinearOperator(int h, int w) : height(h), width(w) {}

        virtual ~LinearOperator() {}

        LinearOperator(const LinearOperator &) = default;
        LinearOperator &operator=(const LinearOperator &) = default;

    public:
        inline int Height() const { return height;}
        inline int Width() const {return width;}
        bool IsSquare() const {return (width==height);}

        virtual void Apply(const Vector<T> &a, Vector<T> &r, T factor = 1.) const = 0;


        // overload if wanted, otherwise fallback exception called when used
        virtual void ApplyT(const Vector<T> &a, Vector<T> &result, T factor = 1.) const
        {
            throw "LinearOperator::ApplyT called in base class, needs to be overridden in derived class";
        }

        // overload if wanted, otherwise fallback exception called when used
        virtual void ApplyH(const Vector<T> &a, Vector<T> &result, T factor = 1.) const
        {
            throw "LinearOperator::ApplyT called in base class, needs to be overridden in derived class";
        }

        virtual void Print(std::ostream& os) const
        {
            os << "[LinearOperator, size " << height << " x " << width << "]" << std::flush;
        }
    };

    enum class TMatrixStorage
    {
        ROW_MAJOR,
        COL_MAJOR
    };

    template <typename T = double>
    class Matrix : public LinearOperator<T>
    {
    protected:
        T *data;
        int alloc_size;
        TMatrixStorage storage;

    public:
        using LinearOperator<T>::height;
        using LinearOperator<T>::width;

        Matrix() : LinearOperator<T>(), data(nullptr), alloc_size(0), storage(TMatrixStorage::ROW_MAJOR) {}

        Matrix(int h, TMatrixStorage st = TMatrixStorage::ROW_MAJOR)
            : LinearOperator<T>(h), data(nullptr), alloc_size(h*h), storage(st)
        {
            if (h > 0)
            {
                data = new T[h * h];
            }
        }

        Matrix(int h, int w, TMatrixStorage st = TMatrixStorage::ROW_MAJOR)
            : LinearOperator<T>(h, w), data(nullptr), alloc_size(h*w), storage(st)
        {
            if (h > 0 && w > 0)
            {
                data = new T[h * w];
            }
        }

        virtual ~Matrix()
        {
            if (data != nullptr)
                delete[] data;
        }

        Matrix(const Matrix &m) : LinearOperator<T>(m)
        {
            storage = m.storage;
            alloc_size = m.alloc_size;
            if (alloc_size > 0)
            {
                data = new T[alloc_size];

                for (int i = 0; i < height * width; i++)
                    data[i] = m.data[i];
            }
            else
                data = nullptr;
        }

        Matrix &operator=(const Matrix &m)
        {
            if (this == &m)
                return *this;

            storage = m.storage;
            LinearOperator<T>::operator=(m);
            if (alloc_size < m.height*m.width)
                AllocateMemory(m.height*m.width);

            for (int i = 0; i < height * width; i++)
                data[i] = m.data[i];

            return *this;
        }

        inline T &operator()(int i, int j) 
        {
#ifndef NDEBUG
            if (i < 0 || i >= height || j < 0 || j >= width)
                throw "Error in Matrix::operator(): out of range";
#endif
            return data[storage == TMatrixStorage::ROW_MAJOR ? i * width + j : j * height + i]; 
        }
        
        inline const T &operator()(int i, int j) const 
        { 
#ifndef NDEBUG
            if (i < 0 || i >= height || j < 0 || j >= width)
                throw "Error in Matrix::operator(): out of range";
#endif
            return data[storage == TMatrixStorage::ROW_MAJOR ? i * width + j : j * height + i]; 
        }

        void SetAll(T val)
        {
            for (int i = 0; i < height * width; i++)
                data[i] = val;
        }

        void SetRandom(double min=0., double max=1., int seed=0);

        void SetStorage(TMatrixStorage st)
        {
            storage = st;
        }

        void AllocateMemory(int alloc_size_)
        {
            #ifndef NDEBUG
            if (alloc_size < height*width)
                std::cout << "Matrix<T>::AllocateMemory: Warning: alloc_size smaller than height*width";
            #endif

            if (data != nullptr) delete [] data;
            if (alloc_size_ > 0)
                data = new T[alloc_size_];
            else
                data = nullptr;
            alloc_size = alloc_size_;
        }

        void SetSize(int h, int w)
        {
            if (h * w > alloc_size)
                AllocateMemory(h*w);
            height = h;
            width = w;
        }

        void SetDiag(T val)
        {
            for (int i = 0; i < std::min(height, width); i++)
                operator()(i, i) = val;
        }

        Matrix &operator+=(const Matrix &mat2)
        {
#ifndef NDEBUG
            if (this->height != mat2.height || this->width != mat2.width)
                throw "Error in Matrix::operator+=: Matrices must have same size";
#endif

            for (int i = 0; i < height * width; i++)
                data[i] += mat2.data[i];

            return *this;
        }

        Matrix &operator*=(T factor)
        {
            for (int i = 0; i < height * width; i++)
                data[i] *= factor;

            return *this;
        }

        Matrix& Transpose()
        {
            std::swap(height,width);
            if (storage == TMatrixStorage::COL_MAJOR) storage = TMatrixStorage::ROW_MAJOR;
            else storage = TMatrixStorage::COL_MAJOR;
            return *this;
        }

        virtual void Apply(const Vector<T> &x, Vector<T> &r, T factor = 1.) const override
        {
#ifndef NDEBUG
            if (x.Size() != width)
                throw "Matrix::Apply dimensions don't fit";
#endif

            r.SetSize(height);

            for (int i = 0; i < height; ++i)
            {
                r(i) = 0.0;
                for (int j = 0; j < width; ++j)
                {
                    r(i) += operator()(i, j) * x(j);
                }
                r(i) *= factor;
            }
        }

        virtual void ApplyT(const Vector<T> &x, Vector<T> &r, T factor = 1.) const override
        {
#ifndef NDEBUG
            if (x.Size() != height)
                throw "Matrix::ApplyT dimensions don't fit";
#endif

            r.SetSize(width);

            for (int i = 0; i < width; ++i)
            {
                r(i) = 0.0;
                for (int j = 0; j < height; ++j)
                {
                    r(i) += operator()(j, i) * x(j);
                }
                r(i) *= factor;
            }
        }

        void ApplyH(const Vector<T> &x, Vector<T> &r, T factor = 1.) const override
        {
#ifndef NDEBUG
            if (x.Size() != height)
                throw "Matrix::ApplyH dimensions don't fit";
#endif

            r.SetSize(width);

            for (int i = 0; i < width; ++i)
            {
                r(i) = 0.0;
                for (int j = 0; j < height; ++j)
                {
                    r(i) += Conjugate(operator()(j, i)) * x(j);
                }
                r(i) *= factor;
            }
        }



        void Print(std::ostream& os) const override
        {
            if (storage == TMatrixStorage::ROW_MAJOR) os << "RM" ;
            else os << "CM";
            os << "[";
            for (int i=0; i<height; i++)
            {
                for (int j=0; j<width; j++)
                    os << (*this)(i,j) << " ";
                os << std::endl;
            }
            os << "]" << std::flush;
        }
    };

    template <>
    void Matrix<double>::
    SetRandom(double min, double max, int seed)
    {
        std::uniform_real_distribution<double> unif(min, max);
        std::default_random_engine re;
        re.seed(seed);
        for (int i = 0; i < height * width; i++)
            data[i] = unif(re);
    }

    template <>
    void Matrix<std::complex<double>>::
    SetRandom(double min, double max, int seed)
    {
        std::uniform_real_distribution<double> unif(min, max);
        std::default_random_engine re;
        re.seed(seed);
        for (int i = 0; i < height * width; i++)
            data[i] = std::complex<double>(unif(re), unif(re));
    }

    template <>
    void Matrix<int>::
    SetRandom(double min, double max, int seed)
    {
        int imin = min; int imax = max;
        std::uniform_int_distribution<> unif(imin, imax);
        std::default_random_engine re;
        re.seed(seed);
        for (int i = 0; i < height * width; i++)
            data[i] = unif(re);
    }

    template <typename T>
    void Matrix<T>:: 
    SetRandom(double min, double max, int seed)
    {
        std::cout << "Matrix<T>::SetRandom not implemented for template type" << std::endl;
        std::cout << "Setting to zero" << std::endl;
        SetAll(0.);
    }
    
    // matrix-vector product c = a.b
    template <typename T>
    void DotProduct(const Matrix<T>& a, const Vector<T>& b, Vector<T>& c)
    {
#ifndef NDEBUG
        if(a.Width() !=b.Size())
            throw "Error in DotProduct: Sizes don't match";
#endif

        c.SetSize(a.Height());
        c.SetAll(0);
        for (int i=0; i<a.Height(); i++)
            for (int j=0; j<a.Width(); j++)
                c(i) += a(i,j)*b(j);

    }

    // matrix-matrix product c = a.b
    template <typename T>
    void DotProduct(const Matrix<T>& a, const Matrix<T>& b, Matrix<T>& c)
    {
        #ifndef NDEBUG
        if(a.Width() !=b.Height())
            throw "Error in DotProduct: Sizes don't match";
        #endif

        c.SetSize(a.Height(),b.Width());
        c.SetAll(0);
        for (int i=0; i<a.Height(); i++)
            for (int j=0; j<b.Width(); j++)
                for (int k=0; k<a.Width(); k++)
                    c(i,j) += a(i,k)*b(k,j);

    }

}


template <typename T>
inline std::ostream& operator<<(std::ostream& os, const SC::LinearOperator<T>& m)
{
    m.Print(os);
    return os;
}

