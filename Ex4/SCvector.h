#pragma once

#include <iostream>
#include <complex>
#include <random>

namespace SC
{

// for complex numbers
inline int Conjugate(int x) { return x; }
inline double Conjugate(double x) { return x; }
inline std::complex<double> Conjugate(std::complex<double> x) { return std::conj(x); }


// vector class
template<typename T=double>
class Vector
{
private:
    T* data;
    int size;
    int alloc_size;

public:
    //Vector() { size=0; alloc_size = 0; data = nullptr;}
    Vector() : size(0), alloc_size(0), data(nullptr)
    { ; }

    Vector(int size_) : size(size_), alloc_size(size_), data(nullptr)
    { 
        if (alloc_size > 0) data = new T[alloc_size];
    }

    Vector(int size_, int alloc_size_) : size(size_), alloc_size(alloc_size_), data(nullptr)
    { 
        if (alloc_size < size) alloc_size = size;
        if (alloc_size > 0) data = new T[alloc_size];
    }

    ~Vector() 
    {
        if (data != nullptr) delete[] data; 
    }

    Vector(const Vector& vec) : size(vec.size), alloc_size(vec.alloc_size)
    {
        if (alloc_size > 0)
        {
            data = new T[alloc_size];
            for (int i = 0; i < size; i++)
            {
                data[i] = vec.data[i];
            }
        }
        else
            data = nullptr;
    }

    Vector& operator= (const Vector& vec)
    { 
        if (this == &vec)
            return *this;

        // AllocateMemory(vec.alloc_size);
        SetSize(vec.size);
        
        for (int i=0; i<size; i++)
        {
            data[i] = vec.data[i];
        }
        return *this;
    }

    Vector& operator= (Vector&& vec)
    { 
        std::swap(size,vec.size);
        std::swap(alloc_size,vec.alloc_size);
        std::swap(data,vec.data);
        return *this;
    }


    inline int Size() const { return size;}

    void AllocateMemory(int alloc_size_)
    {
        if (data != nullptr) delete [] data;
        if (alloc_size_ > 0)
            data = new T[alloc_size_];
        else
            data = nullptr;
        alloc_size = alloc_size_;
        size = alloc_size_;
    }


    void SetSize(int size_)
    {
        if (size_ > alloc_size)
            AllocateMemory(size_);
        size = size_;
    }

    // operators () and [], with range check in debug mode
    T& operator()(int i) 
    { 
        #ifndef NDEBUG
        if(i < 0 || i >= size)
            throw "Error in Vector::operator(): out of range";
        #endif
        return data[i]; 
    }    
    
    const T& operator()(int i) const 
    { 
        #ifndef NDEBUG
        if(i < 0 || i >= size)
            throw "Error in Vector::operator(): out of range";
        #endif
        return data[i]; 
    }    


    // 
    void SetAll(T val) { for (int i=0; i<size; i++) data[i] = val; }


    void SetRandom(double min=0., double max=1., int seed=0);

    // *this = (*this) + vec
    void Add(const Vector &vec2)
    {
        #ifndef NDEBUG
        if(this->size != vec2.size)
            throw "Error in Vector::Add: Vectors must have same size";
        #endif
        
        for (int i=0; i<size; i++)
            data[i] += vec2.data[i];
    }

    // *this = factor * (*this)
    void Mult(T factor)
    {
        for (int i=0; i<size; i++)
            data[i] *= factor;
    }

    
    // *this = (*this) + factor * vec2
    void AddMultiple(T factor, const Vector &vec2)
    {
        #ifndef NDEBUG
        if(this->size != vec2.size)
            throw "Error in Vector::AddMultiple: Vectors must have same size";
        #endif

        for (int i=0; i<size; i++)
            data[i] += factor*vec2.data[i];
    }

    Vector& operator+=(const Vector &vec2)
    {
        #ifndef NDEBUG
        if(this->size != vec2.size)
            throw "Error in Vector::operator+=: Vectors must have same size";
        #endif

        for (int i=0; i<size; i++)
            data[i] += vec2.data[i];

        return *this;
    }

    Vector& operator-=(const Vector &vec2)
    {
        #ifndef NDEBUG
        if(this->size != vec2.size)
            throw "Error in Vector::operator+=: Vectors must have same size";
        #endif
        
        for (int i=0; i<size; i++)
            data[i] -= vec2.data[i];

        return *this;
    }

    Vector& operator*=(T factor)
    {
        for (int i=0; i<size; i++)
            data[i] *= factor;

        return *this;
    }



    double Norm() const
    {
        double sum = 0;
        for (int i=0; i<size; i++)
        {
            sum += std::norm(data[i]);
        }
        return sqrt(sum);
    }

    double NormSqr() const
    {
        double sum = 0;
        for (int i=0; i<size; i++)
        {
            sum += std::norm(data[i]);
        }
        return sum;
    }

    void Print(std::ostream& os) const
    {
        os << "[";
        for (int i=0; i<size; i++)
            os << data[i] << " ";
        os << "]" << std::flush;
    }
};


template <>
void Vector<double>::
SetRandom(double min, double max, int seed)
{
    std::uniform_real_distribution<double> unif(min, max);
    std::default_random_engine re;
    re.seed(seed);
    for (int i = 0; i < size; i++)
        data[i] = unif(re);
}

template <>
void Vector<std::complex<double>>::
SetRandom(double min, double max, int seed)
{
    std::uniform_real_distribution<double> unif(min, max);
    std::default_random_engine re;
    re.seed(seed);
    for (int i = 0; i < size; i++)
        data[i] = std::complex<double>(unif(re), unif(re));
}

template <>
void Vector<int>::
SetRandom(double min, double max, int seed)
{
    int imin = min; int imax = max;
    std::uniform_int_distribution<> unif(imin, imax);
    std::default_random_engine re;
    re.seed(seed);
    for (int i = 0; i < size; i++)
        data[i] = unif(re);
}

template <typename T>
void Vector<T>:: 
SetRandom(double min, double max, int seed)
{
    std::cout << "SetRandom not implemented for template type" << std::endl;
    std::cout << "Setting to zero" << std::endl;
    SetAll(0.);
}

// inner product
template <typename T>
T InnerProduct(const Vector<T>& a, const Vector<T>& b)
{
    #ifndef NDEBUG
    if(a.Size() !=b.Size())
        throw "Error in InnerProduct: Vectors must have same size";
    #endif
    
    T product = 0;
    for (int i=0; i<a.Size(); i++)
        product += b(i)*Conjugate(a(i));

    return product;
}


} // namespace SC


template <typename T>
inline std::ostream& operator<<(std::ostream& os, const SC::Vector<T>& vec)
{
    vec.Print(os);
    return os;
}

template <typename T>
inline SC::Vector<T> operator*(const T& factor, const SC::Vector<T>& vec)
{
    SC::Vector<T> multvec(vec);
    multvec *= factor;
    return multvec;
}

template <typename T>
inline SC::Vector<T> operator+(const SC::Vector<T>& vec1, const SC::Vector<T>& vec2)
{
    SC::Vector<T> sumvec(vec1);
    sumvec += vec2;
    return sumvec;
}
