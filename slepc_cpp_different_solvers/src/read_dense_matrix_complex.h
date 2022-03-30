/*
Author: Praveen Kalarickel Ramakrishnan
This file contains the function "CompMatrixM" to read in the
complex matrix from data files contating the real and imaginary parts
*/

#include <cassert>
#include <fstream>
#include <vector>

int number_of_lines_of_file(std::string inputfilename);
void CompMatrixM(int& n, std::string inputfilename_real, std::string inputfilename_imag, std::vector<double> &M)
{
    n = number_of_lines_of_file(inputfilename_real);
    M.resize(2*n*n);

    std::ifstream inputfile_real;
    inputfile_real.open(inputfilename_real.c_str());
    assert(inputfile_real.is_open());

    std::ifstream inputfile_imag;
    inputfile_imag.open(inputfilename_imag.c_str());
    assert(inputfile_imag.is_open());

    for(int iloop=0; iloop < 2*n*n; iloop++)
    {
        inputfile_real >> M[iloop];
        inputfile_imag >> M[++iloop];
    }
    inputfile_real.close();
    inputfile_imag.close();
}

int number_of_lines_of_file(std::string inputfilename)
{
   int N_lines = 0;
   std::string dummy_string;

   std::ifstream inputfile;
   inputfile.open(inputfilename.c_str());
   assert(inputfile.is_open());
   while(std::getline(inputfile, dummy_string))
	   ++N_lines;
   inputfile.close();
   return N_lines;

}
