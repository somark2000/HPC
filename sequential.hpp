#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>
#include <cstdlib>   // for rand()

#include <mpi.h>

using namespace std;  

namespace program_options {

  struct Options {
    unsigned int mpi_mode;  
    std::string name;
    size_t N;
    size_t M;
    size_t iters;
    void print() const {
      std::printf("mpi_mode: %u\n", mpi_mode);    
      std::printf("name: %s\n", name.c_str());
      std::printf("N: %zu\n", N);
      std::printf("M: %zu\n", M);
      std::printf("iters: %zu\n", iters);
    }
  };

  auto parse(int argc, char *argv[]) {
    if (argc != 6)
      throw std::runtime_error("unexpected number of arguments");
    Options opts;
    if (std::string(argv[1]) == std::string("seq"))
      opts.mpi_mode = 1;
    else if( std::string(argv[1]) == std::string("par"))
      opts.mpi_mode = 2;
    else
    throw std::runtime_error("invalid parameter for mpi_mode (valid are 'seq' and 'par')");
    opts.name = argv[2];
    if (std::sscanf(argv[3], "%zu", &opts.N) != 1 && opts.N >= 2)
      throw std::runtime_error("invalid parameter for N");
    if (std::sscanf(argv[4], "%zu", &opts.M) != 1 && opts.M >= 2)
      throw std::runtime_error("invalid parameter for M");
    if (std::sscanf(argv[5], "%zu", &opts.iters) != 1 && opts.iters != 0)
      throw std::runtime_error("invalid parameter for iters");
    return opts;
  }

} // namespace program_options

//matrix for storing our cells
typedef vector < vector < char > >  Matrix;

char Stensil(int nx, int ny, Matrix &M, int i, int j){
    // count all the alive neighbours 
    int sum = M[(i-2)][(j-2)] + M[(i-1)][(j)] + M[(i-1)][(j+1)] + M[(i)][(j-2)] + M[(i)][(j+1)] + M[(i+2)][(j-2)] + M[(i+2)][(j)] + M[(i+2)][(j+2)%ny];
    // cout<<"sum="<<sum<<"\n";
    sum -= (8 * '0'); // char -> int
    // first term is if the cell is alive and second term is if it is dead -> end result will determine whether the cell will be alive or not
    int survive = (M[i][j]-'0')*(sum==2 || sum==3) + ('1'-M[i][j])*(sum==3);
    // cout<<"sunt in stensil\t survive="<<survive<<"\t sum="<<sum<<"\n";
    return '0' + survive;
}

void UpdateMatrix(int nx, int ny, Matrix &Matric_Old, Matrix &Matrix_New){
    for (int i = 2; i < nx-2; i++)
        for (int j = 2; j < ny-2; j++){
            Matrix_New[i][j] = Stensil(nx, ny, Matric_Old, i, j);
        }
}

void SequentialIterator(int nx, int ny, Matrix &Matrix_Old, Matrix &Matrix_New, int generations){
    for (int generation = 1; generation <= generations; generation++){
        UpdateMatrix (nx,  ny,  Matrix_Old, Matrix_New);
        Matrix_Old = Matrix_New;
    }
}

//an input generator for at least one type of initial configurations
void InputGenerator(int nx, int ny, Matrix &Matrix_Generated, Matrix &Matrix_New, float chance = 0.28){
    // we denote the dead cells with char '0' and the alive cells with char '1'
    // using char as it is a byte datatype for the communication part down the road
    for (int i = 0; i < nx; i++){
        for (int j = 0; j < ny; j++){
            Matrix_Generated[i][j] = '0' + ( float(rand()) / RAND_MAX <= chance );
            Matrix_New[i][j] = '0';
        }
    }
}

void PrintOutput(int nx, int ny, Matrix M){
    int contor = 0;
    for (int i = 0; i < nx; i++)
        for (int j = 0; j < nx; j++)
            contor += M[i][j] - '0';
    cout << contor << " x 1" << endl;
    cout << (nx * ny - contor) << " x 0" << endl;
}

int sequential(program_options::Options opts){
    Matrix Matrix_New, Matrix_Old;
    Matrix_Old.resize(opts.N+2);
    Matrix_New.resize(opts.N+2);
    for(size_t i=0;i<opts.N;++i){
        Matrix_New[i].resize(opts.M+2);
        Matrix_Old[i].resize(opts.M+2);
    }
    InputGenerator(opts.N+2, opts.M+2, Matrix_Old, Matrix_New);
    SequentialIterator(opts.N+2, opts.M+2, Matrix_Old, Matrix_New, opts.iters);
    // cout<<"sunt in sw\n";
    PrintOutput(opts.N+2, opts.M+2, Matrix_New);
    return 0;
}