#ifndef SEQUENTIAL_H
#define SEQUENTIAL_H

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
#include "options.hpp"

using namespace std;  

namespace sequential_solver{
    
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
        Matrix_Old.resize(opts.N+4);
        Matrix_New.resize(opts.N+4);
        for(size_t i=0;i<opts.N+4;++i){
            Matrix_New[i].resize(opts.M+4);
            Matrix_Old[i].resize(opts.M+4);
        }
        cout<<"Sequential solver started----------------------------------------------------\n";
        InputGenerator(opts.N+4, opts.M+4, Matrix_Old, Matrix_New);
        cout<<"Input generation finished----------------------------------------------------\n";
        SequentialIterator(opts.N+4, opts.M+4, Matrix_Old, Matrix_New, opts.iters);
        cout<<"Iteration finished-----------------------------------------------------------\n";
        PrintOutput(opts.N+4, opts.M+4, Matrix_New);
        return 0;
    }
} // namespace sequential_solver



#endif
