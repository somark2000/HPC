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
#include <sstream>

#include <mpi.h>
#include "options.hpp"

using namespace std;  

//matrix for storing our cells
typedef vector < vector < char > >  Matrix;

namespace sequential_solver{

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

    void UpdateMatrix(int nx, int ny, Matrix &Matrix_Old, Matrix &Matrix_New){
        Matrix_New = Matrix_Old;
        for (int i = 2; i < nx-2; i++)
            for (int j = 2; j < ny-2; j++){
                Matrix_New[i][j] = Stensil(nx, ny, Matrix_Old, i, j);
            }
    }

    void UpdateBorder(int nx, int ny, Matrix &Matrix_Old){
        for (int i = 2; i < nx - 2; i++){
            Matrix_Old[i][ny-2] = Matrix_Old[i][2];
            Matrix_Old[i][ny-1] = Matrix_Old[i][3];
            Matrix_Old[i][0] = Matrix_Old[i][ny-4];
            Matrix_Old[i][1] = Matrix_Old[i][ny-3];
        }
        for (int i = 2; i < ny - 2; i++){
            Matrix_Old[nx-2][i] = Matrix_Old[2][i];
            Matrix_Old[nx-1][i] = Matrix_Old[3][i];
            Matrix_Old[0][i] = Matrix_Old[nx-4][i];
            Matrix_Old[1][i] = Matrix_Old[nx-3][i];
        }
        Matrix_Old[0][0] = Matrix_Old[nx-4][ny-4];
        Matrix_Old[0][1] = Matrix_Old[nx-4][ny-3];
        Matrix_Old[1][0] = Matrix_Old[nx-3][ny-4];
        Matrix_Old[1][1] = Matrix_Old[nx-3][ny-3];
        
        Matrix_Old[0][ny-2] = Matrix_Old[nx-4][2];
        Matrix_Old[0][ny-1] = Matrix_Old[nx-4][3];
        Matrix_Old[1][ny-2] = Matrix_Old[nx-3][2];
        Matrix_Old[1][ny-1] = Matrix_Old[nx-3][3];

        Matrix_Old[nx-2][0] = Matrix_Old[2][ny-4];
        Matrix_Old[nx-2][1] = Matrix_Old[2][ny-3];
        Matrix_Old[nx-1][0] = Matrix_Old[3][ny-4];
        Matrix_Old[nx-1][1] = Matrix_Old[3][ny-3];

        Matrix_Old[nx-2][ny-2] = Matrix_Old[2][2];
        Matrix_Old[nx-2][ny-1] = Matrix_Old[2][3];
        Matrix_Old[nx-1][ny-2] = Matrix_Old[3][2];
        Matrix_Old[nx-1][ny-1] = Matrix_Old[3][3];
    }

    void SequentialIterator(int nx, int ny, Matrix &Matrix_Old, Matrix &Matrix_New, int generations){
        for (int generation = 1; generation <= generations; generation++){
            UpdateBorder (nx, ny, Matrix_Old);
            UpdateMatrix (nx, ny, Matrix_Old, Matrix_New);
            Matrix_Old = Matrix_New;
        }
    }

    //an input generator for at least one type of initial configurations
    void InputGenerator(int nx, int ny, Matrix &Matrix_Generated, Matrix &Matrix_New, float chance, int verification){
        // we denote the dead cells with char '0' and the alive cells with char '1'
        // using char as it is a byte datatype for the communication part down the road
        if(verification==1){
            //read specs to .txt file
            ifstream file;
            ostringstream oss;
            oss << "seq_" << nx-4 << "x" << ny-4 << "_gen_-1.txt";
            string name = oss.str();
            file.open(name);
            if (!file.is_open()) {
                std::cerr << "File could not be opened!" << std::endl;
            }
            string line;
            int i=0;
            char c;
            while (getline(file, line) && i<nx-4){
                for (int k = 0; k < ny-4; k++){ //do not consider the endl at the end
                    c = line[k];        
                    Matrix_Generated[i+2][k+2] = c;
                    Matrix_New[i+2][k+2] = '0';
                }
                i++;
            }
        }
        else{
            for (int i = 2; i < nx-2; i++){
                for (int j = 2; j < ny-2; j++){
                    Matrix_Generated[i][j] = '0' + ( float(rand()) / RAND_MAX <= chance );
                    Matrix_New[i][j] = '0';
                }
            }
        }
    }

    void PrintOutput(int nx, int ny, Matrix M){
        int contor = 0;
        for (int i = 0; i < nx; i++)
            for (int j = 0; j < ny; j++)
                contor += M[i+2][j+2] - '0';
        cout << contor << " x 1" << endl;
        cout << (nx * ny - contor) << " x 0" << endl;
        float survival = float(contor)/(nx * ny) ;
        cout << "survival rate is " << survival <<"\n";
    }

    void SaveState(int nx, int ny, Matrix M, int iters){
        //write specs to .txt file
        ofstream file;
        ostringstream oss;
        oss << "seq_" << nx-4 << "x" << ny-4 << "_gen_"<< iters << ".txt";
        string name = oss.str();
        file.open(name);
        if (!file.is_open()) {
            std::cerr << "File could not be opened!" << std::endl;
        }

        for (int i = 0; i < nx-4; i++){
            for (int j = 0; j < ny-4; j++)
                file << M[i+2][j+2];
            file << '\n';
        }
        file.close();
    }

    void GlobalInputGenerator(program_options::Options opts){
        //read specs to .txt file
        ifstream file;
        ostringstream oss;
        oss << "seq_" << opts.N << "x" << opts.M << "_gen_-1.txt";
        string name = oss.str();
        file.open(name);
        if (!file.is_open()) {
            Matrix Matrix_Old;
            Matrix_Old.resize(opts.N+4);
            for(size_t i=0;i<opts.N+4;++i){
                Matrix_Old[i].resize(opts.M+4);
            }
            for (size_t i = 0; i < opts.N+4; i++){
                for (size_t j = 0; j < opts.M+4; j++){
                    if( i > 1 and i < opts.N + 2 and j > 1 and j < opts.M + 2)
                    Matrix_Old[i][j] = '0' + ( float(rand()) / RAND_MAX <= opts.chance );
                    else
                    Matrix_Old[i][j] = '0';
                }
            }
            SaveState(opts.N+4,opts.M+4,Matrix_Old,-1);
        }
        
    }

    int sequential(program_options::Options opts){
        // time measurement points
        double starttime, endtime, runtime;
        Matrix Matrix_New, Matrix_Old;
        Matrix_Old.resize(opts.N+4);
        Matrix_New.resize(opts.N+4);
        for(size_t i=0;i<opts.N+4;++i){
            Matrix_New[i].resize(opts.M+4);
            Matrix_Old[i].resize(opts.M+4);
        }
        cout<<"Sequential solver started----------------------------------------------------\n";
        InputGenerator(opts.N+4, opts.M+4, Matrix_Old, Matrix_New,opts.chance,opts.verification);
        cout<<"Input generation finished----------------------------------------------------\n";

        //saving initial state
        if(opts.verification==true) SaveState(opts.N+4,opts.M+4,Matrix_Old,0);
        
        // start timer
        starttime = MPI_Wtime();
        
        SequentialIterator(opts.N+4, opts.M+4, Matrix_Old, Matrix_New, opts.iters);
        cout<<"Iteration finished-----------------------------------------------------------\n";
        
        // end timer
        endtime = MPI_Wtime();
        runtime = (endtime-starttime)*1000000; //runtime transformed in microseconds
        std::cout << "The sequential program executed in " << runtime << " microseconds.\n";
        //write specs to .txt file
        ofstream file;
        ostringstream oss;
        oss << "seq_" << opts.N << "x" << opts.M << "_times.dat";
        string name = oss.str();
        file.open(name,ios::app);
        file << "The sequential program executed in " << runtime << " microseconds.\n";

        //saving final state
        if(opts.verification==true) SaveState(opts.N+4,opts.M+4,Matrix_Old,opts.iters);
        PrintOutput(opts.N, opts.M, Matrix_New);
        return 0;
    }
} // namespace sequential_solver



#endif
