#pragma once

#include <iostream>
#include <complex>

#include "SCvector.h"

namespace SC
{


    template <typename T = double>
    class HessenbergMatrix : public LinearOperator<T>
    {
    protected:
        T *data;
        int alloc_size;
        
        
    public:
        using LinearOperator<T>::height;
        using LinearOperator<T>::width;

        int DataSize(int k) const {return ((k+1) * k)/2 + k;}

        HessenbergMatrix() : LinearOperator<T>(), data(nullptr), alloc_size(0) {}

        HessenbergMatrix(int k)
            : LinearOperator<T>(k+1,k), data(nullptr), alloc_size(k)
        {
            if (k > 0)
            {
                data = new T[DataSize(k)];
                SetAll(0);
            }
        }


        virtual ~HessenbergMatrix()
        {
            if (data != nullptr)
                delete[] data;
        }

        HessenbergMatrix(const HessenbergMatrix &m) : LinearOperator<T>(m)
        {
            alloc_size = width;
            if (width > 0)
            {
                data = new T[DataSize(width)];

                for (int i = 0; i < DataSize(width); i++)
                    data[i] = m.data[i];
            }
            else
                data = nullptr;
        }

        HessenbergMatrix &operator=(const HessenbergMatrix &m)
        {
            if (this == &m)
                return *this;

            if (alloc_size != m.alloc_size)
            {
                alloc_size = m.alloc_size;
                if (data != nullptr)
                    delete[] data;
                LinearOperator<T>::operator=(m);
                if (alloc_size > 0)
                    data = new T[ DataSize(alloc_size)];
                else
                    data = nullptr;
            }

            for (int i = 0; i <  DataSize(alloc_size); i++)
                data[i] = m.data[i];

            return *this;
        }

        inline T Get(int i, int j) const
        { 
            if (i > j+1) return 0;
            return data[DataSize(j)+i]; 
        }
        
        inline void Set(int i, int j, T val) 
        { 
            #ifndef NDEBUG
            if (i > j+1) throw "HessenbergMatrix::Set";
            #endif
            data[DataSize(j)+i] = val; 
        }
        inline void AddTo(int i, int j, T val) 
        { 
            #ifndef NDEBUG
            if (i > j+1) throw "HessenbergMatrix::AddTo";
            #endif
            data[DataSize(j)+i] += val; 
        }

        void SetAll(T val)
        {
            for (int i = 0; i < DataSize(width); i++)
                data[i] = val;
        }

        void SetSize(int w)
        {
            if (w > alloc_size)
            {
                alloc_size = w;
                if (data != nullptr)
                    delete[] data;
                if (w > 0)
                    data = new T[DataSize(w)];
                else
                    data = nullptr;
            }
            height = w+1;
            width = w;
        }


        void Apply(const Vector<T> &x, Vector<T> &r, T factor = 1.) const override
        {
            #ifndef NDEBUG
            if (x.Size() != width)
                throw "Matrix::Apply dimensions don't fit";
            #endif

            r.SetSize(height);

            for (int i = 0; i < height; ++i)
            {
                r(i) = 0.0;
                for (int j = std::max(0,i-1); j < width; ++j)
                {
                    r(i) += Get(i, j) * x(j);
                }
                r(i) *= factor;
            }
        }



        void Print(std::ostream& os) const override
        {
            os << "[";
            for (int i=0; i<height; i++)
            {
                for (int j=0; j<width; j++)
                    os << Get(i,j) << " ";
                os << std::endl;
            }
            os << "]" << std::flush;
        }
    };


}


