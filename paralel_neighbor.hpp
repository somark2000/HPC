#ifndef PARALELN_H
#define PARALELN_H

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

//matrix for storing our cells
typedef vector < vector < char > >  Matrix;

namespace paraleln_solver{


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

    void CreateGraphComm(int rank,int &rank_n, int &rank_s, int &rank_e, int &rank_w, int &rank_ne, int &rank_nw, int &rank_se, int &rank_sw,
                        MPI_Comm comm_2d, MPI_Comm &comm_graph, int r, int c){
        // fint main neighbors
        MPI_Cart_shift(comm_2d, 1, 1, &rank_n, &rank_s); //dir 1 is Y
        MPI_Cart_shift(comm_2d, 0, 1, &rank_w, &rank_e); //dir 0 is X
        
        //coordinates of rank in virtual topology
        array<int, 2> coord = {-1, -1};    
        MPI_Cart_coords(comm_2d, rank, 2, std::data(coord));
        
        //find corner neighbors
        array<int, 2> coord_corner = {-1, -1};
        coord_corner = {coord[0]-1, coord[1]+1};
        if(coord_corner[0] < 0) coord_corner[0]+=c;
        else if (coord_corner[0] > c) coord_corner[0]-=c;
        if(coord_corner[1] < 0) coord_corner[1]+=r;
        else if (coord_corner[1] > r) coord_corner[1]-=r;
        MPI_Cart_rank(comm_2d, std::data(coord_corner), &rank_sw);
        coord_corner = {coord[0]+1, coord[1]-1};
        if(coord_corner[0] < 0) coord_corner[0]+=c;
        else if (coord_corner[0] > c) coord_corner[0]-=c;
        if(coord_corner[1] < 0) coord_corner[1]+=r;
        else if (coord_corner[1] > r) coord_corner[1]-=r;
        MPI_Cart_rank(comm_2d, std::data(coord_corner), &rank_ne);
        coord_corner = {coord[0]+1, coord[1]+1};
        if(coord_corner[0] < 0) coord_corner[0]+=c;
        else if (coord_corner[0] > c) coord_corner[0]-=c;
        if(coord_corner[1] < 0) coord_corner[1]+=r;
        else if (coord_corner[1] > r) coord_corner[1]-=r;
        MPI_Cart_rank(comm_2d, std::data(coord_corner), &rank_se);
        coord_corner = {coord[0]-1, coord[1]+1};
        if(coord_corner[0] < 0) coord_corner[0]+=c;
        else if (coord_corner[0] > c) coord_corner[0]-=c;
        if(coord_corner[1] < 0) coord_corner[1]+=r;
        else if (coord_corner[1] > r) coord_corner[1]-=r;
        MPI_Cart_rank(comm_2d, std::data(coord_corner), &rank_sw);

        //create new communicator
        int sources[1],weights_s[1]; 
        int destinations[8], weights[8]; 

        /* get my main communication partners */ 
        destinations[0] = rank_w; weights[0] = 2; 
        destinations[1] = rank_s; weights[1] = 2; 
        destinations[2] = rank_e; weights[2] = 2; 
        destinations[3] = rank_n; weights[3] = 2; 
        
        /* get my communication partners along diagonals */ 
        destinations[4] = rank_se; weights[4] = 1; 
        destinations[5] = rank_sw; weights[5] = 1; 
        destinations[6] = rank_ne; weights[6] = 1; 
        destinations[7] = rank_nw; weights[7] = 1; 
        
        sources[0] = rank; 
        weights_s[0] = 2;
        MPI_Dist_graph_create_adjacent(comm_2d, 1, sources, weights_s, 8, destinations, weights, MPI_INFO_NULL, 0, &comm_graph); 
    }


    void CommunicateBorders(int nx, int ny, Matrix &Matrix_Old, Matrix &Matrix_New, int rang, MPI_Comm comm_2d, int r, int c){
        //start comm
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
        char *Buffer_sw = (char*)malloc(sizeof(char)*9);
        char *Buffer_se = (char*)malloc(sizeof(char)*9);
        char *Buffer_nw = (char*)malloc(sizeof(char)*9);
        char *Buffer_ne = (char*)malloc(sizeof(char)*9);
        
        // find neighbours
        int rank_n, rank_s, rank_e, rank_w;
        int rank_ne, rank_nw, rank_se, rank_sw;
        int displacement=0;

        //create new communicator for the hood architeture
        MPI_Comm comm_graph;
        CreateGraphComm(rang,rank_n,rank_s,rank_e,rank_w,rank_ne,rank_nw,rank_se,rank_sw,comm_2d,comm_graph,r,c);

        // Prepare data for MPI_Alltoallv
        int sendcounts[8];  // Number of elements to send to each process
        int recvcounts[8];  // Number of elements to receive from each process
        int sdispls[8];     // Displacements in send buffer
        int rdispls[8];     // Displacements in receive buffer
        char *Buffer_send = (char*)malloc(sizeof(char)*(4*(nx+ny)));
        char *Buffer_recv = (char*)malloc(sizeof(char)*(4*(nx+ny)));
        // cout<<"I am rang "<<rang<<" and my neighbeours N S E W are:"<< rank_n<<" "<<rank_s<<" "<< rank_e<<" "<<rank_w<<"\n";

        // send west
        // Buffer1 = Col 2 because of ghost layers
        // Buffer2 = Col 3 because of ghost layers
        for (int i = 2; i < nx - 2; i++){
            Buffer1_w[i-2] = Matrix_Old[i][2];
            Buffer2_w[i-2] = Matrix_Old[i][3];
        }
        sendcounts[0] = 2*(nx-4);
        recvcounts[0] = 2*(nx-4);
        sdispls[0] = displacement;
        rdispls[0] = displacement;
        displacement += nx-4;
        strcpy(Buffer_send + displacement,Buffer1_w);
        displacement += nx-4;
        strcpy(Buffer_send + displacement,Buffer2_w);
        
        // send south
        // Buffer1 = ROW nx-4 because of ghost layers
        // Buffer2 = ROW nx-3 because of ghost layers
        for (int i = 2; i < ny - 2; i++){
            Buffer1_s[i-2] = Matrix_Old[nx-4][i];
            Buffer2_s[i-2] = Matrix_Old[nx-3][i];
        }
        sendcounts[1] = 2*(ny-4);
        recvcounts[1] = 2*(ny-4);
        sdispls[1] = displacement;
        rdispls[1] = displacement;
        strcpy(Buffer_send + displacement,Buffer1_s);
        displacement += ny-4;
        strcpy(Buffer_send + displacement,Buffer2_s);
        displacement += ny-4;

        // send east
        // Buffer1 = Col ny-4 because of ghost layers
        // Buffer2 = Col ny-3 because of ghost layers
        for (int i = 2; i < nx - 2; i++){
            Buffer1_e[i-2] = Matrix_Old[i][ny-4];
            Buffer2_e[i-2] = Matrix_Old[i][ny-3];
        }
        sendcounts[2] = 2*(nx-4);
        recvcounts[2] = 2*(nx-4);
        sdispls[2] = displacement;
        rdispls[2] = displacement;
        strcpy(Buffer_send + displacement,Buffer1_e);
        displacement += nx-4;
        strcpy(Buffer_send + displacement,Buffer2_e);
        displacement += nx-4;
        
        // send north
        // Buffer1 = ROW 2 because of ghost layers
        // Buffer2 = ROW 3 because of ghost layers
        for (int i = 2; i < ny - 2; i++){
            Buffer1_n[i-2] = Matrix_Old[2][i];
            Buffer2_n[i-2] = Matrix_Old[3][i];
        }
        sendcounts[3] = 2*(ny-4);
        recvcounts[3] = 2*(ny-4);
        sdispls[3] = displacement;
        rdispls[3] = displacement;
        strcpy(Buffer_send + displacement,Buffer1_n);
        displacement += ny-4;
        strcpy(Buffer_send + displacement,Buffer2_n);
        displacement += ny-4;

        // send corners
        //send SW corner
        Buffer_sw[0] = Matrix_Old[nx-4][2];
        Buffer_sw[1] = Matrix_Old[nx-4][3]; 
        Buffer_sw[2] = Matrix_Old[nx-3][2];
        Buffer_sw[3] = Matrix_Old[nx-3][3];
        sendcounts[4] = 4;
        recvcounts[4] = 4;
        sdispls[4] = displacement;
        rdispls[4] = displacement;
        strcpy(Buffer_send + displacement,Buffer_sw);
        displacement += 4;
        // cout<<"preparing to sent south west from "<< rang<< " to "<<dest<<" with coord_corner col: "<<coord_corner[0]<<" and row: "<<coord_corner[1]<<"\n";

        //send SE corner
        Buffer_se[0] = Matrix_Old[nx-4][ny-4];
        Buffer_se[1] = Matrix_Old[nx-4][ny-3]; 
        Buffer_se[2] = Matrix_Old[nx-3][ny-4];
        Buffer_se[3] = Matrix_Old[nx-3][ny-3];
        sendcounts[5] = 4;
        recvcounts[5] = 4;
        sdispls[5] = displacement;
        rdispls[5] = displacement;
        strcpy(Buffer_send + displacement,Buffer_se);
        displacement += 4;
        // cout<<"preparing to sent south east from "<< rang<< " to "<<dest<<" with coord_corner col: "<<coord_corner[0]<<" and row: "<<coord_corner[1]<<"\n";

        //send NE corner
        Buffer_ne[0] = Matrix_Old[2][ny-4];
        Buffer_ne[1] = Matrix_Old[2][ny-3]; 
        Buffer_ne[2] = Matrix_Old[3][ny-4];
        Buffer_ne[3] = Matrix_Old[3][ny-3];
        sendcounts[6] = 4;
        recvcounts[6] = 4;
        sdispls[6] = displacement;
        rdispls[6] = displacement;
        strcpy(Buffer_send + displacement,Buffer_ne);
        displacement += 4;
        // cout<<"preparing to sent north east from "<< rang<< " to "<<dest<<" with coord_corner col: "<<coord_corner[0]<<" and row: "<<coord_corner[1]<<"\n";

        //send NW corner
        Buffer_nw[0] = Matrix_Old[2][2];
        Buffer_nw[1] = Matrix_Old[2][3]; 
        Buffer_nw[2] = Matrix_Old[3][2];
        Buffer_nw[3] = Matrix_Old[3][3];
        sendcounts[7] = 4;
        recvcounts[7] = 4;
        sdispls[7] = displacement;
        rdispls[7] = displacement;
        strcpy(Buffer_send + displacement,Buffer_nw);
        displacement += 4; // end
        // cout<<"preparing to sent north west from "<< rang<< " to "<<dest<<" with coord_corner col: "<<coord_corner[0]<<" and row: "<<coord_corner[1]<<"\n";
        
        cout<<"I am rank "<<rang<<" before comm\n";
        // Perform MPI_Alltoallv
        MPI_Neighbor_alltoallv(Buffer_send, sendcounts, sdispls, MPI_CHAR, Buffer_recv, recvcounts, rdispls, MPI_CHAR, comm_graph);
        cout<<"I am rank "<<rang<<" after comm\n";
        
        
        // recieves start here
        displacement = 0;
        // recieve west
        strncpy(Buffer1_w,Buffer_recv + displacement,nx-4);
        displacement += nx-4;
        strncpy(Buffer2_w,Buffer_recv + displacement,nx-4);
        displacement += nx-4;
        // Buffer1 = Col 0 it is a ghost layers
        // Buffer2 = Col 1 it is a ghost layers
        for (int i = 0; i < nx - 4; i++){
            Matrix_Old[i+2][0] = Buffer1_w[i];
            Matrix_Old[i+2][1] = Buffer2_w[i];
        }
        // cout<<"recieved  west from "<< source<<" to "<<rang<<"\n";

        // recieve south
        strncpy(Buffer1_s,Buffer_recv + displacement,ny-4);
        displacement += ny-4;
        strncpy(Buffer2_s,Buffer_recv + displacement,ny-4);
        displacement += ny-4;
        // Buffer1 = ROW nx-2 it is a ghost layers
        // Buffer2 = ROW nx-1 it is a ghost layers
        for (int i = 0; i < ny - 4; i++){
            Matrix_Old[nx-2][i+2] = Buffer1_s[i];
            Matrix_Old[nx-1][i+2] = Buffer2_s[i];
        }
        // cout<<"recieved south from "<< source<<" to "<<rang<<"\n";
            
        // recieve east
        strncpy(Buffer1_e,Buffer_recv + displacement,nx-4);
        displacement += nx-4;
        strncpy(Buffer2_e,Buffer_recv + displacement,nx-4);
        displacement += nx-4;
        // Buffer1 = Col ny-2 it is a ghost layers
        // Buffer2 = Col ny-1 it is a ghost layers
        for (int i = 0; i < nx - 4; i++){
            Matrix_Old[i+2][ny-2] = Buffer1_e[i];
            Matrix_Old[i+2][ny-1] = Buffer2_e[i];
        }
        // cout<<"recieved east from "<< source<<" to "<<rang<<"\n";

        // recieve north
        strncpy(Buffer1_n,Buffer_recv + displacement,ny-4);
        displacement += ny-4;
        strncpy(Buffer2_n,Buffer_recv + displacement,ny-4);
        displacement += ny-4;
        // Buffer1 = ROW 0 it is a ghost layers
        // Buffer2 = ROW 1 because of ghost layers
        for (int i = 0; i < ny - 4; i++){
            Matrix_Old[0][i+2] = Buffer1_n[i];
            Matrix_Old[1][i+2] = Buffer2_n[i];
        }
        // cout<<"recieved north from "<< source<<" to "<<rang<<"\n";
        
        //recieve corners
        //recieve SW corner
        strncpy(Buffer_sw,Buffer_recv + displacement,size_t(4));
        displacement += 4;
        Matrix_Old[nx-2][0] = Buffer_sw[0];
        Matrix_Old[nx-2][1] = Buffer_sw[1]; 
        Matrix_Old[nx-1][0] = Buffer_sw[2];
        Matrix_Old[nx-1][1] = Buffer_sw[3];
        // cout<<"recieved to corner SW from "<<source<<" to "<<rang<<"\n";

        //recieve SE corner
        strncpy(Buffer_se,Buffer_recv + displacement,size_t(4));
        displacement += 4;
        Matrix_Old[nx-2][ny-2] = Buffer_se[0];
        Matrix_Old[nx-2][ny-1] = Buffer_se[1];
        Matrix_Old[nx-1][ny-2] = Buffer_se[2]; 
        Matrix_Old[nx-1][ny-1] = Buffer_se[3];
        // cout<<"recieved to corner SE from "<<source<<" to "<<rang<<"\n";

        //recieve NE corner
        strncpy(Buffer_ne,Buffer_recv + displacement,size_t(4));
        displacement += 4;
        Matrix_Old[0][ny-2] = Buffer_ne[0];
        Matrix_Old[0][ny-1] = Buffer_ne[1]; 
        Matrix_Old[1][ny-2] = Buffer_ne[2];
        Matrix_Old[1][ny-1] = Buffer_ne[3];
        // cout<<"recieved to corner NE from "<<source<<" to "<<rang<<"\n";

        //recieve NW corner
        strncpy(Buffer_nw,Buffer_recv + displacement,size_t(4));
        displacement += 4;
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

    int paralel_neighbor(program_options::Options opts, int rank, int numprocesses){
        // let openMPI take care of splitting into suitable subdomain grid
        int ndims = 2; //2D grid
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
        cout<<"Generation of the domain "<< rank <<" successfull------------------------------------------------\n";
        
        Iterate(N_row, N_col, Matrix_Old, Matrix_New, opts.iters, rank, comm_2d,r,c);
        cout<<"Iteration of the domain "<< rank <<" successfull------------------------------------------------\n";
        Matrix result;
        result.resize(opts.N);
        for(size_t i=0;i<opts.N;++i)
            result[i].resize(opts.M);
        // CollectData(Matrix_New,result,rank,numprocesses,opts.N,opts.N,row_div,row_mod,col_div,col_mod,r,c,N_row,N_col,comm_2d);
        if (rank==0){
            PrintOutput(opts.N, opts.M, result,0);
        }

        return 0;
    }
}// namespace paraleln_solver

#endif
