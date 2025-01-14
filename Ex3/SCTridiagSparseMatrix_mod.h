#pragma once

#include <iostream>
#include <complex>

#include "SCvector.h"
#include "SCmatrix.h"
#include "SCSparseMatrix.h"

namespace SC
{

    template <typename T = double>
    class TridiagSparseMatrix : public SparseMatrix<T>
    {
        using LinearOperator<T>::height;
        using LinearOperator<T>::width;

        using SparseMatrix<T>::data;
        using SparseMatrix<T>::colind;
        using SparseMatrix<T>::rowptr;

    public:
        TridiagSparseMatrix() = delete;

        ~TridiagSparseMatrix()
        {
            ;
        }

        TridiagSparseMatrix(const TridiagSparseMatrix &m) = delete;

        TridiagSparseMatrix(size_t size, Vector<T> &diag, Vector<T> &subdiag, Vector<T> &superdiag) : SparseMatrix<T>(size, size)
        {
            size_t nze = 3 * size - 2;
            data = new T[nze];
            colind = new int[nze];
            rowptr = new int[height + 1];

            size_t ii = 0;
            rowptr[0] = 0;
            data[ii] = diag(0);
            colind[ii++] = 0;
            data[ii] = superdiag(0);
            colind[ii++] = 1;

            rowptr[1] = 2;
            data[ii] = subdiag(0);
            colind[ii++] = 0;
            data[ii] = diag(1);
            colind[ii++] = 1;
            data[ii] = superdiag(1);
            colind[ii++] = 2;
            for (int i = 2; i < height; i++)
            {
                rowptr[i] = rowptr[i - 1] + 3;
                data[ii] = subdiag(i - 1);
                colind[ii++] = i - 1;
                data[ii] = diag(i);
                colind[ii++] = i;
                if (i < height - 1)
                {
                    data[ii] = superdiag(i);
                    colind[ii++] = i + 1;
                }
            }
            rowptr[height] = nze;
        }

        // initialize sparse matrix graph, data is 0.0
        TridiagSparseMatrix(size_t size) : SparseMatrix<T>(size, size)
        {
            size_t nze = 3 * size - 2;
            data = new T[nze];
            colind = new int[nze];
            rowptr = new int[height + 1];

            size_t ii = 0;
            rowptr[0] = 0;
            data[ii] = T(0);
            colind[ii++] = 0;
            data[ii] = T(0);
            colind[ii++] = 1;

            rowptr[1] = 2;
            data[ii] = T(0);
            colind[ii++] = 0;
            data[ii] = T(0);
            colind[ii++] = 1;
            data[ii] = T(0);
            colind[ii++] = 2;
            for (int i = 2; i < height; i++)
            {
                rowptr[i] = rowptr[i - 1] + 3;
                data[ii] = T(0);
                colind[ii++] = i - 1;
                data[ii] = T(0);
                colind[ii++] = i;
                if (i < height - 1)
                {
                    data[ii] = T(0);
                    colind[ii++] = i + 1;
                }
            }
            rowptr[height] = nze;
        }

        // initialize sparse matrix graph, data from diag, subdiag and superdiag vectors
        TridiagSparseMatrix(size_t size, const Vector<T> &diag, const Vector<T> &subdiag, const Vector<T> &superdiag) : SparseMatrix<T>(size, size)
        {
            if (diag.Size() != size || subdiag.Size() + 1 != size || superdiag.Size() + 1 != size)
                throw "TridiagSparseMatrix: initialization with vectors of non-matching size";

            size_t nze = 3 * size - 2;
            data = new T[nze];
            colind = new int[nze];
            rowptr = new int[height + 1];

            size_t ii = 0;
            rowptr[0] = 0;
            data[ii] = diag(0);
            colind[ii++] = 0;
            data[ii] = superdiag(0);
            colind[ii++] = 1;

            rowptr[1] = 2;
            data[ii] = subdiag(0);
            colind[ii++] = 0;
            data[ii] = diag(1);
            colind[ii++] = 1;
            data[ii] = superdiag(1);
            colind[ii++] = 2;
            for (int i = 2; i < height; i++)
            {
                rowptr[i] = rowptr[i - 1] + 3;
                data[ii] = subdiag(i - 1);
                colind[ii++] = i - 1;
                data[ii] = diag(i);
                colind[ii++] = i;
                if (i < height - 1)
                {
                    data[ii] = superdiag(i);
                    colind[ii++] = i + 1;
                }
            }
            rowptr[height] = nze;
        }

        TridiagSparseMatrix &operator=(const TridiagSparseMatrix &m) = delete;

        T &Data(int i) { return data[i]; }
        const T &Data(int i) const { return data[i]; }

        int &ColInd(int i) { return colind[i]; }
        const int &ColInd(int i) const { return colind[i]; }

        int &RowPtr(int i) { return rowptr[i]; }
        const int &RowPtr(int i) const { return rowptr[i]; }

        inline T Get(int i, int j) const
        {
            for (int k = rowptr[i]; k < rowptr[i + 1]; k++)
            {
                if (colind[k] == j)
                    return data[k];
            }
            return 0;
        }

        inline void Set(int i, int j, T value)
        {
            for (int k = rowptr[i]; k < rowptr[i + 1]; k++)
            {
                if (colind[k] == j)
                {
                    data[k] = value;
                    return;
                }
            }
            throw "TridiagSparseMatrix: try to set zero element";
        }

        void SetAll(T val)
        {
            for (int i = 0; i < rowptr[height]; i++)
                data[i] = val;
        }

        void Apply(const Vector<T> &x, Vector<T> &r, T factor = 1.) const override
        {

            if (x.Size() != width)
                throw "TridiagSparseMatrix::Apply dimensions don't fit";
            r.SetSize(height);

            for (int i = 0; i < height; ++i)
            {
                r(i) = 0.0;
                for (int j = rowptr[i]; j < rowptr[i + 1]; ++j)
                {
                    r(i) += data[j] * x(colind[j]);
                }
                r(i) *= factor;
            }
        }

        void ApplyT(const Vector<T> &x, Vector<T> &r, T factor = 1.) const override
        {

            if (x.Size() != height)
                throw "TridiagSparseMatrix::Apply dimensions don't fit";
            r.SetSize(width);
            r.SetAll(0.);

            for (int i = 0; i < height; ++i)
            {
                for (int j = rowptr[i]; j < rowptr[i + 1]; ++j)
                {
                    r(colind[j]) += data[j] * x(i);
                }
            }
            r *= factor;
        }

        void ApplyH(const Vector<T> &x, Vector<T> &r, T factor = 1.) const override
        {

            if (x.Size() != height)
                throw "TridiagSparseMatrix::Apply dimensions don't fit";
            r.SetSize(width);
            r.SetAll(0.);

            for (int i = 0; i < height; ++i)
            {
                for (int j = rowptr[i]; j < rowptr[i + 1]; ++j)
                {
                    r(colind[j]) += Conjugate(data[j]) * x(i);
                }
            }
            r *= factor;
        }

        void Print(std::ostream &os) const override
        {
            os << "[";
            for (int i = 0; i < height; i++)
            {
                for (int j = rowptr[i]; j < rowptr[i + 1]; j++)
                    os << colind[j] << ": " << (data)[j] << " -- ";
                os << std::endl;
            }
            os << "]" << std::flush;
        }
    };

}