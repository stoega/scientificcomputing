#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>

int main()
{
    // write to standard output with std::cout, compare printf
    // char static_string [] = {"HelloWorld static"};
    // printf("%s with printf\n", static_string);
    // std::cout << "Using cout: "; 
    // std::cout << static_string;
    // std::cout << std::endl;

    // write double, int, bool to std::cout
    double d = 4.5;
    int j = 3;
    bool b = true;
    std::cout << d;
    std::cout << std::endl;

    std::cout << j;
    std::cout << std::endl;
    std::cout << b;
    std::cout << std::endl;

    // use std::string class for dynamic-size strings ----------
    std::string std_string;
    std_string = "HelloWorld std::string";
    std_string += " ... ";

    std::cout << std_string;
    std::cout << std::endl;
    // get c-type char* pointer from string
    printf("%s with printf\n", std_string.c_str());

    // operators -----------------------------------------------
    // std::cout << "hello"; is equivalent to
    operator<<(std::cout, "hello\n");

    // return type of operator<< is std::ostream&
    std::ostream& myout = operator<<(std::cout, "hello\n");

    // therefore, we can "<<" several strings/values/etc, e.g.
    operator<<(operator<<(std::cout, "hello "), "world\n");
    // or, more conveniently
    std::cout << "hello " << "world\n";


    // formatting ---------------------------------------------
    std::cout << std::fixed << std::setprecision(5) ;
    std::cout << d << std::endl;


    // sstream -- streaming to strings
    std::stringstream mystringstream;
    mystringstream << "Hello World";
    mystringstream << "\n... again\n";

    std::string s2 = mystringstream.str();
    std::cout << "Stringstream contains:\n";
    std::cout << s2;
    std::cout << std::endl;

    // convert int to string
    std::stringstream intstream;
    int theanswer = 42;
    intstream << theanswer;
    std::string s3 = intstream.str();

    // ofstream -- write to file
    std::ofstream outfile("outputfile.txt");
    // above defaults to:
    // std::ofstream outfile("outputfile.txt",std::ios_base::out);
    // append to existing text in file
    // std::ofstream outfile("outputfile.txt",std::ios_base::out|std::ios_base::app);
    outfile << "Hello File! 1234";
    outfile.close();

    // ifstream -- read from file
    std::ifstream infile("outputfile.txt"); // has been generated above
    std::string filestring1, filestring2;
    int i;
    char c;
    // read first string (up to space)
    infile >> filestring1;
    // read second string (up to space)
    infile >> filestring2;
    // read next item as char
    infile >> c;
    // read next item as int
    infile >> i;
    std::cout << "Read from file via infile: ";
    std::cout << filestring1;
    std::cout << filestring2;
    std::cout << c;
    std::cout << i;
    std::cout << std::endl;




}