#ifndef PARALEL_H
#define PARALEL_H

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

namespace paralel_solver{

    //matrix for storing our cells
    typedef vector < vector < char > >  Matrix;

    char Stensil(int nx, int ny, Matrix &M, int i, int j){
        // count all the alive neighbours 
        int sum = M[(i-2)][(j-2)] + M[(i-1)][(j)] + M[(i-1)][(j+1)] + M[(i)][(j-2)] + M[(i)][(j+1)] + M[(i+2)][(j-2)] + M[(i+2)][(j)] + M[(i+2)][(j+2)];
        // cout<<"sum="<<sum<<"\n";
        sum -= (8 * '0'); // char -> int
        // first term is if the cell is alive and second term is if it is dead -> end result will determine whether the cell will be alive or not
        int survive = (M[i][j]-'0')*(sum==2 || sum==3) + ('1'-M[i][j])*(sum==3);
        // cout<<"sunt in stensil\t survive="<<survive<<"\t sum="<<sum<<"\n";
        return '0' + survive;
    }

    void PrintOutput(int nx, int ny, Matrix M, int rank){
        int contor = 0;
        for (int i = 0; i < nx; i++)
            for (int j = 0; j < nx; j++)
                contor += M[i][j] - '0';
        cout<<"for rank "<<rank<<" we have\n";
        cout << contor << " x 1" << endl;
        cout << (nx * ny - contor) << " x 0" << endl;
    }

    // return number of cols and rows in a subdomain characterized by (cartesian) coordinates
    std::tuple<int, int> GetSubdomainDims(std::array<int, 2> coord, int row_div, int row_mod, int col_div, int col_mod){
    int N_row = row_div;
    int N_col = col_div;

    // distribute resuming points evenly along domains
    if(coord[0] < col_mod){ N_col += 1;}
    if(coord[1] < row_mod){ N_row += 1;}

    // add ghost layers (4 ghost layers per axis (2-2))  since periodic topology
    N_row += 4;
    N_col += 4;
    
    return {N_col, N_row};
    }

    void UpdateMatrix(int nx, int ny, Matrix &Matrix_Old, Matrix &Matrix_New){
        for (int i = 2; i < nx-2; i++)
            for (int j = 2; j < ny-2; j++){
                Matrix_New[i][j] = Stensil(nx, ny, Matrix_Old, i, j);
            }
    }

    void CommunicateBorders(int nx, int ny, Matrix &Matrix_Old, Matrix &Matrix_New, int rang, MPI_Comm comm_2d, int r, int c){
        //start comm
        // MPI_Barrier(comm_2d);
        MPI_Status status1, status2;
        # define mnax(a, b) ((a > b) ? (a) : (b))
        int MAX_DIM = mnax (nx, ny);
        char *Buffer1_s = (char*)malloc(sizeof(char)*MAX_DIM);
        char *Buffer2_s = (char*)malloc(sizeof(char)*MAX_DIM);
        char *Buffer1_n = (char*)malloc(sizeof(char)*MAX_DIM);
        char *Buffer2_n = (char*)malloc(sizeof(char)*MAX_DIM);
        char *Buffer1_e = (char*)malloc(sizeof(char)*MAX_DIM);
        char *Buffer2_e = (char*)malloc(sizeof(char)*MAX_DIM);
        char *Buffer1_w = (char*)malloc(sizeof(char)*MAX_DIM);
        char *Buffer2_w = (char*)malloc(sizeof(char)*MAX_DIM);
        // find neighbours
        int dest, source;
        int rank_n, rank_s, rank_e, rank_w;
        MPI_Cart_shift(comm_2d, 1, 1, &rank_n, &rank_s); //dir 1 is Y
        MPI_Cart_shift(comm_2d, 0, 1, &rank_w, &rank_e); //dir 0 is X
        //find corner neighbours
        int rank_ne, rank_nw, rank_se, rank_sw;
        array<int, 2> coord_corner = {-1, -1};

        // cout<<"I am rang "<<rang<<" and my neighbeours N S E W are:"<< rank_n<<" "<<rank_s<<" "<< rank_e<<" "<<rank_w<<"\n";
        
        //communicate west

        // send west
        // Buffer1 = Col 2 because of ghost layers
        for (int i = 2; i < nx - 2; i++)
            Buffer1_w[i-2] = Matrix_Old[i][2];
        // Buffer2 = Col 3 because of ghost layers
        for (int i = 2; i < nx - 2; i++)
            Buffer2_w[i-2] = Matrix_Old[i][3];
        // send two cols to west
        dest = rank_w;
        MPI_Send(Buffer1_w, nx-4, MPI_CHAR, dest, 0, comm_2d);
        MPI_Send(Buffer2_w, nx-4, MPI_CHAR, dest, 0, comm_2d);
        // cout<<"sent west from "<< rang<< " to "<<dest<<"\n";

        // send south
        // Buffer1 = ROW nx-4 because of ghost layers
        for (int i = 2; i < ny - 2; i++)
            Buffer1_s[i-2] = Matrix_Old[nx-4][i];
        // Buffer2 = ROW nx-3 because of ghost layers
        for (int i = 2; i < ny - 2; i++)
            Buffer2_s[i-2] = Matrix_Old[nx-3][i];
        // send two cols to west
        dest = rank_s;
        MPI_Send(Buffer1_s, ny-4, MPI_CHAR, dest, 0, comm_2d);
        MPI_Send(Buffer2_s, ny-4, MPI_CHAR, dest, 0, comm_2d);
        // cout<<"sent south from "<< rang <<" to "<<dest<<"\n";

        // send east
        // Buffer1 = Col ny-4 because of ghost layers
        for (int i = 2; i < nx - 2; i++)
            Buffer1_e[i-2] = Matrix_Old[i][ny-4];
        // Buffer2 = Col ny-3 because of ghost layers
        for (int i = 2; i < nx - 2; i++)
            Buffer2_e[i-2] = Matrix_Old[i][ny-3];
        // send two cols to west
        dest = rank_e;
        MPI_Send(Buffer1_e, nx-4, MPI_CHAR, dest, 0, comm_2d);
        MPI_Send(Buffer2_e, nx-4, MPI_CHAR, dest, 0, comm_2d);
        // cout<<"sent west from "<< rang <<" to "<<dest<<"\n";

        // send north
        // Buffer1 = ROW 2 because of ghost layers
        for (int i = 2; i < ny - 2; i++)
            Buffer1_n[i-2] = Matrix_Old[2][i];
        // Buffer2 = ROW 3 because of ghost layers
        for (int i = 2; i < ny - 2; i++)
            Buffer2_n[i-2] = Matrix_Old[3][i];
        // send two cols to west
        dest = rank_n;
        MPI_Send(Buffer1_n, ny-4, MPI_CHAR, dest, 0, comm_2d);
        MPI_Send(Buffer2_n, ny-4, MPI_CHAR, dest, 0, comm_2d);
        // cout<<"sent north from "<< rang <<" to "<<dest<<"\n";

        // recieves start here

        // recieve west
        source = rank_w;
        MPI_Recv(Buffer1_w, nx-4, MPI_CHAR, source, 0, comm_2d, &status1);
        MPI_Recv(Buffer2_w, nx-4, MPI_CHAR, source, 0, comm_2d, &status2);
        // Buffer1 = Col 0 it is a ghost layers
        for (int i = 0; i < nx - 4; i++)
            Matrix_Old[i+2][0] = Buffer1_w[i];
        // Buffer2 = Col 1 it is a ghost layers
        for (int i = 0; i < nx - 4; i++)
            Matrix_Old[i+2][1] = Buffer2_w[i];
        // cout<<"recieved  west from "<< source<<" to "<<rang<<"\n";

        // recieve south
        source = rank_s;
        MPI_Recv(Buffer1_s, ny-4, MPI_CHAR, source, 0, comm_2d, &status1);
        MPI_Recv(Buffer2_s, ny-4, MPI_CHAR, source, 0, comm_2d, &status2);
        // Buffer1 = ROW nx-2 it is a ghost layers
        for (int i = 0; i < ny - 4; i++){
            // cout<<"\t recieve south\t"<<i<<" "<<nx-1<<"\n";
            Matrix_Old[nx-2][i+2] = Buffer1_s[i];
        }
        // Buffer2 = ROW nx-1 it is a ghost layers
        for (int i = 0; i < ny - 4; i++)
            Matrix_Old[nx-1][i+2] = Buffer2_s[i];
        // cout<<"recieved south from "<< source<<" to "<<rang<<"\n";
            
        // recieve east
        source = rank_e;
        MPI_Recv(Buffer1_e, nx-4, MPI_CHAR, source, 0, comm_2d, &status1);
        MPI_Recv(Buffer2_e, nx-4, MPI_CHAR, source, 0, comm_2d, &status2);
        // Buffer1 = Col ny-2 it is a ghost layers
        for (int i = 0; i < nx - 4; i++)
            Matrix_Old[i+2][ny-2] = Buffer1_e[i];
        // Buffer2 = Col ny-1 it is a ghost layers
        for (int i = 0; i < nx - 4; i++)
            Matrix_Old[i+2][ny-1] = Buffer2_e[i];
        // cout<<"recieved east from "<< source<<" to "<<rang<<"\n";

        // recieve north
        source = rank_n;
        MPI_Recv(Buffer1_n, ny-4, MPI_CHAR, source, 0, comm_2d, &status1);
        MPI_Recv(Buffer2_n, ny-4, MPI_CHAR, source, 0, comm_2d, &status2);
        // Buffer1 = ROW 0 it is a ghost layers
        for (int i = 0; i < ny - 4; i++)
            Matrix_Old[0][i+2] = Buffer1_n[i];
        // Buffer2 = ROW 1 because of ghost layers
        for (int i = 0; i < ny - 4; i++)
            Matrix_Old[1][i+2] = Buffer2_n[i];
        // cout<<"recieved north from "<< source<<" to "<<rang<<"\n";
        
        //send corners

        char *Buffer_sw = (char*)malloc(sizeof(char)*4);
        char *Buffer_se = (char*)malloc(sizeof(char)*4);
        char *Buffer_nw = (char*)malloc(sizeof(char)*4);
        char *Buffer_ne = (char*)malloc(sizeof(char)*4);
        //coordinates of rank in virtual topology
        array<int, 2> coord = {-1, -1};    
        MPI_Cart_coords(comm_2d, rang, 2, std::data(coord));
        
        //send SW corner
        coord_corner = {coord[0]-1, coord[1]+1};
        if(coord_corner[0] < 0) coord_corner[0]+=c;
        else if (coord_corner[0] > c) coord_corner[0]-=c;
        if(coord_corner[1] < 0) coord_corner[1]+=r;
        else if (coord_corner[1] > r) coord_corner[1]-=r;
        MPI_Cart_rank(comm_2d, std::data(coord_corner), &rank_sw);
        dest = rank_sw;
        // cout<<"preparing to sent south west from "<< rang<< " to "<<dest<<" with coord_corner col: "<<coord_corner[0]<<" and row: "<<coord_corner[1]<<"\n";
        Buffer_sw[0] = Matrix_Old[nx-4][2];
        Buffer_sw[1] = Matrix_Old[nx-4][3]; 
        Buffer_sw[2] = Matrix_Old[nx-3][2];
        Buffer_sw[3] = Matrix_Old[nx-3][3];
        MPI_Send(Buffer_sw, 4, MPI_CHAR, dest, 0, comm_2d);

        //send SE corner
        coord_corner = {coord[0]+1, coord[1]+1};
        if(coord_corner[0] < 0) coord_corner[0]+=c;
        else if (coord_corner[0] > c) coord_corner[0]-=c;
        if(coord_corner[1] < 0) coord_corner[1]+=r;
        else if (coord_corner[1] > r) coord_corner[1]-=r;
        MPI_Cart_rank(comm_2d, std::data(coord_corner), &rank_se);
        dest = rank_se;
        // cout<<"preparing to sent south east from "<< rang<< " to "<<dest<<" with coord_corner col: "<<coord_corner[0]<<" and row: "<<coord_corner[1]<<"\n";
        Buffer_se[0] = Matrix_Old[nx-4][ny-4];
        Buffer_se[1] = Matrix_Old[nx-4][ny-3]; 
        Buffer_se[2] = Matrix_Old[nx-3][ny-4];
        Buffer_se[3] = Matrix_Old[nx-3][ny-3];
        MPI_Send(Buffer_se, 4, MPI_CHAR, dest, 0, comm_2d);

        //send NE corner
        coord_corner = {coord[0]+1, coord[1]-1};
        if(coord_corner[0] < 0) coord_corner[0]+=c;
        else if (coord_corner[0] > c) coord_corner[0]-=c;
        if(coord_corner[1] < 0) coord_corner[1]+=r;
        else if (coord_corner[1] > r) coord_corner[1]-=r;
        MPI_Cart_rank(comm_2d, std::data(coord_corner), &rank_ne);
        dest = rank_ne;
        Buffer_ne[0] = Matrix_Old[2][ny-4];
        Buffer_ne[1] = Matrix_Old[2][ny-3]; 
        Buffer_ne[2] = Matrix_Old[3][ny-4];
        Buffer_ne[3] = Matrix_Old[3][ny-3];
        MPI_Send(Buffer_ne, 4, MPI_CHAR, dest, 0, comm_2d);
        // cout<<"preparing to sent north east from "<< rang<< " to "<<dest<<" with coord_corner col: "<<coord_corner[0]<<" and row: "<<coord_corner[1]<<"\n";

        //send NW corner
        coord_corner = {coord[0]-1, coord[1]-1};
        if(coord_corner[0] < 0) coord_corner[0]+=c;
        else if (coord_corner[0] > c) coord_corner[0]-=c;
        if(coord_corner[1] < 0) coord_corner[1]+=r;
        else if (coord_corner[1] > r) coord_corner[1]-=r;
        MPI_Cart_rank(comm_2d, std::data(coord_corner), &rank_nw);
        dest = rank_nw;
        Buffer_nw[0] = Matrix_Old[2][2];
        Buffer_nw[1] = Matrix_Old[2][3]; 
        Buffer_nw[2] = Matrix_Old[3][2];
        Buffer_nw[3] = Matrix_Old[3][3];
        MPI_Send(Buffer_nw, 4, MPI_CHAR, dest, 0, comm_2d);
        // cout<<"preparing to sent north west from "<< rang<< " to "<<dest<<" with coord_corner col: "<<coord_corner[0]<<" and row: "<<coord_corner[1]<<"\n";
        
        //recieve corners
        
        //recieve SW corner
        source = rank_sw;
        MPI_Recv(Buffer_sw, 4, MPI_CHAR, source, 0, comm_2d, &status1);
        Matrix_Old[nx-2][0] = Buffer_sw[0];
        Matrix_Old[nx-2][1] = Buffer_sw[1]; 
        Matrix_Old[nx-1][0] = Buffer_sw[2];
        Matrix_Old[nx-1][1] = Buffer_sw[3];
        // cout<<"recieved to corner SW from "<<source<<" to "<<rang<<"\n";

        //recieve SE corner
        source = rank_se;
        MPI_Recv(Buffer_se, 4, MPI_CHAR, source, 0, comm_2d, &status1);
        Matrix_Old[nx-2][ny-2] = Buffer_se[0];
        Matrix_Old[nx-2][ny-1] = Buffer_se[1];
        Matrix_Old[nx-1][ny-2] = Buffer_se[2]; 
        Matrix_Old[nx-1][ny-1] = Buffer_se[3];
        // cout<<"recieved to corner SE from "<<source<<" to "<<rang<<"\n";

        //recieve NE corner
        source = rank_ne;
        MPI_Recv(Buffer_ne, 4, MPI_CHAR, source, 0, comm_2d, &status1);
        Matrix_Old[0][ny-2] = Buffer_ne[0];
        Matrix_Old[0][ny-1] = Buffer_ne[1]; 
        Matrix_Old[1][ny-2] = Buffer_ne[2];
        Matrix_Old[1][ny-1] = Buffer_ne[3];
        // cout<<"recieved to corner NE from "<<source<<" to "<<rang<<"\n";

        //recieve NW corner
        source = rank_nw;
        MPI_Recv(Buffer_nw, 4, MPI_CHAR, source, 0, comm_2d, &status1);
        Matrix_Old[0][0] = Buffer_nw[0];
        Matrix_Old[0][1] = Buffer_nw[1]; 
        Matrix_Old[1][0] = Buffer_nw[2];
        Matrix_Old[1][1] = Buffer_nw[3];
        // cout<<"recieved to corner Nw from "<<source<<" to "<<rang<<"\n";


        // clear
        free(Buffer1_n);
        free(Buffer2_n);
        free(Buffer1_e);
        free(Buffer2_e);
        free(Buffer1_w);
        free(Buffer2_w);
        free(Buffer1_s);
        free(Buffer2_s);
        free(Buffer_sw);
        free(Buffer_nw);
        free(Buffer_ne);
        free(Buffer_se);

    }

    void Iterate(int nx, int ny, Matrix &Matrix_Old, Matrix &Matrix_New, int generations, int rang, MPI_Comm comm_2d, int r, int c){
        for (int generation = 1; generation <= generations; generation++){
            //we first communicate to update the ghost layers for each subdomain after the last iteration
            //also since the ghost layers originally haven't been set in the first iteration we set them 
            CommunicateBorders (nx,  ny,  Matrix_Old, Matrix_New, rang, comm_2d, r, c);
            // cout << "rank "<<rang<<" communicated\n";
            UpdateMatrix (nx,  ny,  Matrix_Old, Matrix_New);
            // cout << "rank "<<rang<<" updated\n";
            Matrix_Old = Matrix_New;
        }
    }

    //an input generator for at least one type of initial configurations
    void  InputGenerator(int nx, int ny, Matrix &Matrix_Generated, Matrix &Matrix_New, float chance = 0.28){
        // we denote the dead cells with char '0' and the alive cells with char '1'
        // using char as it is a byte datatype for the communication part down the road
        // we have to leave empty the ghost layers and later communicate them with their neighbours
        for (int j = 0; j < ny; j++){
                Matrix_Generated[0][j] = '0';
                Matrix_New[0][j] = '0';
                Matrix_Generated[1][j] = '0';
                Matrix_New[1][j] = '0';
                Matrix_Generated[nx-2][j] = '0';
                Matrix_New[nx-2][j] = '0';
                Matrix_Generated[nx-1][j] = '0';
                Matrix_New[nx-1][j] = '0';
        }
        for (int i = 2; i < nx-2; i++){
            Matrix_Generated[i][0] = '0';
            Matrix_Generated[i][1] = '0';
            Matrix_New[i][0] = '0';
            Matrix_New[i][1] = '0';
            for (int j = 2; j < ny-2; j++){
                Matrix_Generated[i][j] = '0' + ( float(rand()) / RAND_MAX <= chance );
                Matrix_New[i][j] = '0';
            }
            Matrix_New[i][ny-2] = '0';
            Matrix_New[i][ny-1] = '0';
            Matrix_Generated[i][ny-2] = '0';
            Matrix_Generated[i][ny-1] = '0';
        }
    }

    //function that connects all the slices of the processes together in the end (for 2D)
    void CollectData(Matrix Matrixx, Matrix &result, int rank, int numprocesses, int N, int M, int row_div, int row_mod, int col_div, int col_mod,
        int r, int c, int N_row, int N_col, MPI_Comm comm_2d){
        MPI_Status status;
        // char *Buffer = (char*)malloc(sizeof(char)*(row_div+1)*(col_div+1)); //maximum size that could appear

        //initialize 2d array to hold the received data
        std::vector<std::vector<Matrix>> data;
        data.resize(r);
        for(int i=0;i<r;++i)
            data[i].resize(c);
        //collect everything at rank 0
        if(rank==0){
            data[0][0] = Matrixx;

            for (int i = 1; i < numprocesses; i++) {
                std::array<int, 2> rec_coords{ -1, -1 };

                MPI_Cart_coords(comm_2d, i, 2, std::data(rec_coords));
                std::tuple <int, int> rw_col = GetSubdomainDims(rec_coords, row_div, row_mod, col_div, col_mod);
                int N_col, N_row;
                N_col = std::get<0>(rw_col);
                N_row = std::get<1>(rw_col);
                char *Buffer = (char*)malloc(sizeof(char)*(N_row)*(N_col)); //maximum size that could appear
                int size = (N_row-4)*(N_col-4);

                MPI_Recv(Buffer, size, MPI_CHAR, i, 0, comm_2d, &status);
                // cout<<"recieved from rank "<<i<<" "<<N_row<<" "<<N_col<<" "<<"\n";
                Matrix recieved;
                recieved.resize(N_row);
                for(int j=0;j<N_row;++j)
                    recieved[j].resize(N_col);
                //rebuild the matrix of the subdomain without the ghostlayers
                for (int j = 0; j < N_row-4; j++){
                    for (int k = 0; k < N_col-4; k++){
                        recieved[j][k] = Buffer[ (j)*(N_col-4) + (k) ];
                    }
                }
                data[rec_coords[1]][rec_coords[0]] = recieved;
                // cout<<"recieved into matrix from rank "<<i<<" with size "<<N_row<<" "<<N_col<<" "<<size<<" data\n"<<std::string(Buffer,size)<<".\n";
                free(Buffer);
            }
            // build the final solution
            int rowoffset = 0;
            for(int i=0;i<r;i++){ //row of data
                // cout<<"survived colletion\n";
                int maxrow = row_div; //number of rows in the current subdomain set
                if(i<row_mod) maxrow++;
                for(int j=0;j<maxrow;j++){ //row of recieved
                    int coloffset = 0;
                    for(int k=0;k<c;k++){ //column of data
                        int maxcol = col_div; //number of rows in the current subdomain set
                        if(k<col_mod) maxcol++;
                        for(int l=0;l<maxcol;l++){
                            // cout<<"doing fun ind subdomains of row "<<i<<"and cols "<<k<<" of inner rows "<<j<<" and inner cols"<<l<<"\n";
                            // cout<<"to be inserted at "<<rowoffset + j<<" "<<coloffset + l<<"\n";
                            // cout<<"the data is "<< data[i][k][j][l] <<"\n";
                            // cout<<"finished at "<<i<<" "<<k<<"\t and pos "<<j<<" "<<l<<"\n to be inserted at "<<rowoffset + j<<" "<<coloffset + l<<"\n" ;
                            // data[i][k][j][l] = '0';
                            // result[ rowoffset + j ][ coloffset + l ]='0';
                            result[ rowoffset + j ][ coloffset + l ] = data[i][k][j][l];
                        }
                        coloffset+=maxcol;
                    }
                }
                rowoffset+=maxrow;
            }
        }
        else{
            char *Buffer = (char*)malloc(sizeof(char)*(N_row)*(N_col)); //maximum size that could appear
            // create a buffer with all the elements of the subdomain without the ghost layers
            for (int i = 0; i < N_row-4; i++){
                for (int j = 0; j < N_col-4; j++){
                    Buffer[ i*(N_col-4) + j ] = Matrixx[i+2][j+2];
                }
            }
            int size = (N_row-4)*(N_col-4);
            // cout<<"sent from rank "<<rank<<" with size "<<N_row<<" "<<N_col<<" "<<size<<" data\n"<<std::string(Buffer,size)<<".\n";
            // printf("%.*s \n",size,Buffer);
            MPI_Send(Buffer,size, MPI_CHAR,0,0,comm_2d);
            free(Buffer);
        }
        // return result;
    }

    int paralel(program_options::Options opts, int rank, int numprocesses){
        // let openMPI take care of splitting into suitable subdomain grid
        int ndims = 2; //2D grid
        // std::array<int, 2> dims = {1, 12};
        std::array<int, 2> dims = {0, 0};
        MPI_Dims_create(numprocesses, ndims, std::data(dims));

        int r = dims[1]; //number or row-subdomain the domain was split
        int c = dims[0]; //number or col-subdomain the domain was split


        //create cartesian virtual geometry
        std::array<int, 2> periods = {true, true}; // periodic topology
        // int reorder = opts.reorder;
        int reorder = false;
        MPI_Comm comm_2d;
        MPI_Cart_create(MPI_COMM_WORLD, ndims, std::data(dims), std::data(periods), reorder, &comm_2d);

        if(rank == 0){
            cout << "2D decomposition done with " << r << " rows and " << c << " columns." << std::endl;
        }
        // cout << "rank " << rank << " is active\n";

        int row_div, row_mod, col_div, col_mod;
        row_div = opts.N / r; //minimum number of rows per slice
        row_mod = opts.N % r; //number of slices with max number of rows
        col_div = opts.M / c; //minimum number of cols per slice
        col_mod = opts.M % c; //number of slices with max number of cols
        /*
        n=17
        17/3 = 5
        17%3 = 2 
        -> 2 rand cu 6 si 3-2 cu 5
        n=18
        18/3=6
        18%3=0
        -> 0 rand cu 7 si 3-0 cu 6
        n=19
        19/3=6
        19%3=1
        -> 1 rand cu 7 si 3-1 cu 6
        */

        //coordinates of rank in virtual topology
        array<int, 2> coord = {-1, -1};    
        MPI_Cart_coords(comm_2d, rank, ndims, std::data(coord));
        // cout<<"I am rang "<<rank<<"and I am at row "<<coord[1]<<" and col "<<coord[0]<<"\n";

        // number of columns and rows in this subdomain (with added ghost layers)
        std::tuple <int, int> rw_col = GetSubdomainDims(coord, row_div, row_mod, col_div, col_mod);
        size_t N_col, N_row;
        N_col = std::get<0>(rw_col);
        N_row = std::get<1>(rw_col);

        //print the parsed options from the command line
        if(rank == 0){
            opts.print();
        }

        //generating the setup for current subdomain
        Matrix Matrix_New, Matrix_Old;
        Matrix_Old.resize(N_row);
        Matrix_New.resize(N_row);
        // cout<<"rank= "<<rank<<" nrow= "<<N_row<<" ncol= "<<N_col<<"\n";
        for(size_t i=0;i<N_row;++i){
            Matrix_New[i].resize(N_col);
            Matrix_Old[i].resize(N_col);
        }
        InputGenerator(N_row, N_col, Matrix_Old, Matrix_New);
        cout<<"Generation of the domain "<< rank <<" successfull------------------------------------------------" << std::endl;
        
        Iterate(N_row, N_col, Matrix_Old, Matrix_New, opts.iters, rank, comm_2d,r,c);
        cout<<"Iteration of the domain "<< rank <<" successfull------------------------------------------------" << std::endl;
        Matrix result;
        result.resize(opts.N);
        for(size_t i=0;i<opts.N;++i)
            result[i].resize(opts.M);
        CollectData(Matrix_New,result,rank,numprocesses,opts.N,opts.N,row_div,row_mod,col_div,col_mod,r,c,N_row,N_col,comm_2d);
        if (rank==0){
            PrintOutput(opts.N, opts.M, result,0);
        }

        return 0;
    }
}// namespace paralel_solver

#endif
