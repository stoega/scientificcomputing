#include <iostream>


int main()
{

    double d = -6.5;// 64 bit double
    unsigned char *c = reinterpret_cast<unsigned char*>(&d); // reinterpret as 8 1byte char

    printf("internal representation of %f is:\n", d);
    for (int i=0; i<8; i++)
    {
        printf("%02x ", c[i]);
    }
    std::cout << std::endl;
    return 0;
}