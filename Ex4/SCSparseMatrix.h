#pragma once

#include <iostream>
#include <complex>

#include "SCvector.h"
#include "SCmatrix.h"

namespace SC
{

    template <typename T = double>
    class SparseMatrix : public LinearOperator<T>
    {
    protected:
        using LinearOperator<T>::height;
        using LinearOperator<T>::width;

        T *data;
        int *colind;
        int *rowptr;

    public:
        SparseMatrix() = delete;

        ~SparseMatrix()
        {
            if (data != nullptr)
            {
                delete[] data;
                delete[] colind;
                delete[] rowptr;
            }
        }

        SparseMatrix(const SparseMatrix &m) : LinearOperator<T>(m)
        {
            int nze = m.rowptr[height];
            data = new T[nze];
            colind = new int[nze];
            rowptr = new int[height+1];

            for (int i = 0; i < nze; i++)
                data[i] = m.data[i];
            for (int i = 0; i < nze; i++)
                colind[i] = m.colind[i];
            for (int i = 0; i <= height; i++)
                rowptr[i] = m.rowptr[i];
            

        }
        
        protected:
        // for use in derived class TridiagSparseMatrix only
        SparseMatrix(int h, int w) : LinearOperator<T>(h, w)
        {
            data = nullptr;
            colind = nullptr;
            rowptr = nullptr;
        }   
        
        public:
        SparseMatrix(int h, int w, Vector<T>& data_v, Vector<int>& colind_v, Vector<int>& rowptr_v) : LinearOperator<T>(h,w)
        {
            int nze = rowptr_v(height);
            data = new T[nze];
            colind = new int[nze];
            rowptr = new int[height+1];

            for (int i = 0; i < nze; i++)
                data[i] = data_v(i);
            for (int i = 0; i < nze; i++)
                colind[i] = colind_v(i);
            for (int i = 0; i <= height; i++)
                rowptr[i] = rowptr_v(i); 

        }

        SparseMatrix(int h, int w, 
            std::vector<int> row, 
            std::vector<int> col, 
            std::vector<T>& val) : LinearOperator<T>(h,w)
        {
            int nze = row.size();
            data = new T[nze];
            colind = new int[nze];
            rowptr = new int[height+1];

            Vector<int> nze_per_row(height);
            nze_per_row.SetAll(0);

            for (int i=0; i<nze; i++)
            {   
                //initialize all column indices by -1
                colind[i] = -1;
                // count the number of non-zero elements per row
                nze_per_row(row[i])++;
            }

            // Step 1
            // generate the rowprt array, using count for non-zero elements
            rowptr[0] = 0;
            for (int i=0; i<height; i++)
                rowptr[i+1] = rowptr[i] + nze_per_row(i);

            // Step 2
            // generate the colind array
            // unsorted, no check for double entries
            for (int k=0; k<nze; k++)
            {
                // the first available column index
                int* colind_p = &(colind[rowptr[row[k]]]);
                // if -1, it has not been used yet - set to col[k]
                while (*colind_p != -1)
                {
                    colind_p++;
                    if (colind_p >= &(colind[rowptr[row[k]+1]]))
                        throw "SparseMatrix: error, not enough nze in row";
                    if (*colind_p == col[k])
                        throw "SparseMatrix: double entry";
                }
                *colind_p = col[k];
            }

            // Step 3
            // set the data values, using Set
            for (int k=0; k<nze; k++)
            {
                this->Set(row[k], col[k], val[k]);
            }

        }

        SparseMatrix &operator=(const SparseMatrix &m)
        {
            if (this == &m)
                return *this;


            delete[] data;
            delete[] colind;
            delete[] rowptr;



            LinearOperator<T>::operator=(m);
            int nze = m.rowptr[height];
            data = new T[nze];
            colind = new int[nze];
            rowptr = new int[height+1];

            for (int i = 0; i < nze; i++)
                data[i] = m.data[i];
            for (int i = 0; i < nze; i++)
                colind[i] = m.colind[i];
            for (int i = 0; i <= height; i++)
                rowptr[i] = m.rowptr[i];
        
            return *this;
        }

        T& Data(int i) { return data[i]; }
        const T& Data(int i) const { return data[i]; }

        int& ColInd(int i) { return colind[i];}
        const int& ColInd(int i) const { return colind[i]; }

        int& RowPtr(int i) { return rowptr[i]; }
        const int& RowPtr(int i) const { return rowptr[i]; }
        
        inline T Get(int i, int j) const
        { 
            for (int k=rowptr[i]; k<rowptr[i+1]; k++)
                {
                    if (colind[k] == j) return data[k];
                }
            return 0;
        }

        inline void Set(int i, int j, T value) 
        { 
            for (int k=rowptr[i]; k<rowptr[i+1]; k++)
                {
                    if (colind[k] == j) 
                    {
                        data[k] = value;
                        return;
                    }
                }
            throw "SparseMatrix: try to set zero element";
        }

        void SetAll(T val)
        {
            for (int i = 0; i < rowptr[height]; i++)
                data[i] = val;
        }


        void Apply(const Vector<T> &x, Vector<T> &r, T factor = 1.) const override
        {

            if (x.Size() != width)
                throw "Matrix::Apply dimensions don't fit";
            r.SetSize(height);

            for (int i = 0; i < height; ++i)
            {
                r(i) = 0.0;
                for (int j = rowptr[i]; j < rowptr[i+1]; ++j)
                {
                    r(i) += data[j] * x(colind[j]);
                }
                r(i) *= factor;
            }
         }

        void ApplyT(const Vector<T> &x, Vector<T> &r, T factor = 1.) const override
        {

            if (x.Size() != height)
                throw "Matrix::Apply dimensions don't fit";
            r.SetSize(width);
            r.SetAll(0.);

            for (int i = 0; i < height; ++i)
            {
                for (int j = rowptr[i]; j < rowptr[i+1]; ++j)
                {
                    r(colind[j]) += data[j] * x(i);
                }
            }
            r *= factor;
         }

        void ApplyH(const Vector<T> &x, Vector<T> &r, T factor = 1.) const override
        {

            if (x.Size() != height)
                throw "Matrix::Apply dimensions don't fit";
            r.SetSize(width);
            r.SetAll(0.);

            for (int i = 0; i < height; ++i)
            {
                for (int j = rowptr[i]; j < rowptr[i+1]; ++j)
                {
                    r(colind[j]) += Conjugate(data[j]) * x(i);
                }
            }
            r *= factor;
         }

        void Print(std::ostream& os) const override
        {
            os << "[";
            for (int i=0; i<height; i++)
            {
                for (int j=rowptr[i]; j<rowptr[i+1]; j++)
                    os << colind[j] << ": " << (data)[j] << " -- ";
                os << std::endl;
            }
            os << "]" << std::flush;
        }

    };


}