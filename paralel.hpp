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
    // unsigned int reorder;
    //TODO: scaling options, chance of alive, 

    void print() const {
      std::printf("mpi_mode: %u\n", mpi_mode);    
      std::printf("name: %s\n", name.c_str());
      std::printf("N: %zu\n", N);
      std::printf("M: %zu\n", M);
      std::printf("iters: %zu\n", iters);
    //   std::printf("reoder: %u\n", reorder);
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
    // if (std::sscanf(argv[6], "%u", &opts.reorder) == 1 || opts.reorder == 0)
    //     throw std::runtime_error("invalid parameter for reorder flag (valid are 0 and 1) and recieved");
    return opts;
  }

} // namespace program_options

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

// return number of cols and rows in a subdomain characterized by (cartesian) coordinates
std::tuple<int, int> GetSubdomainDims(std::array<int, 2> coord, int row_ndom, int row_bdom, int col_ndom, int col_bdom, int rowparts, int colparts){
  int N_row = row_ndom;
  int N_col = col_ndom;

  // distribute resuming points evenly along domains
  if(coord[0] < col_bdom){ N_col += 1;}
  if(coord[1] < row_bdom){ N_row += 1;}

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

void CommunicateBorders(int nx, int ny, Matrix &Matrix_Old, Matrix &Matrix_New, int rang, array<int, 2> coord, MPI_Comm comm_2d, int r, int c){
    //start comm
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Status status1, status2;

    # define mnax(a, b) ((a > b) ? (a) : (b))
    int MAX_DIM = mnax (nx, ny);
    char *Buffer1 = (char*)malloc(sizeof(char)*MAX_DIM);
    char *Buffer2 = (char*)malloc(sizeof(char)*MAX_DIM);
    // find neighbours
    int dest, source;
    int rank_n, rank_s, rank_e, rank_w;
    MPI_Cart_shift(comm_2d, 0, 1, &rank_n, &rank_s);
    MPI_Cart_shift(comm_2d, 1, 1, &rank_w, &rank_e);
    //find corner neighbours
    int rank_ne, rank_nw, rank_se, rank_sw;
    array<int, 2> coord_corner = {-1, -1};


    // send borders

    // send west
    // Buffer1 = Col 2 because of ghost layers
    for (int i = 0; i < nx - 1; i++)
        Buffer1[i] = Matrix_Old[i][2];
    // Buffer2 = Col 3 because of ghost layers
    for (int i = 0; i < nx - 1; i++)
        Buffer2[i] = Matrix_Old[i][3];
    // send two cols to west
    dest = rank_w;
    MPI_Send(Buffer1, nx, MPI_CHAR, dest, 0, MPI_COMM_WORLD);
    MPI_Send(Buffer2, nx, MPI_CHAR, dest, 0, MPI_COMM_WORLD);
    cout<<"sent west from "<< rang<< "to"<<dest<<"\n";

    // send south
    // Buffer1 = ROW nx-3 because of ghost layers
    for (int i = 0; i < ny - 1; i++)
        Buffer1[i] = Matrix_Old[nx-3][i];
    // Buffer2 = ROW nx-4 because of ghost layers
    for (int i = 0; i < ny - 1; i++)
        Buffer2[i] = Matrix_Old[nx-4][i];
    // send two cols to west
    dest = rank_s;
    MPI_Send(Buffer1, ny, MPI_CHAR, dest, 0, MPI_COMM_WORLD);
    MPI_Send(Buffer2, ny, MPI_CHAR, dest, 0, MPI_COMM_WORLD);
    
    // send east
    // Buffer1 = Col ny-3 because of ghost layers
    for (int i = 0; i < nx - 1; i++)
        Buffer1[i] = Matrix_Old[i][ny-3];
    // Buffer2 = Col ny-4 because of ghost layers
    for (int i = 0; i < nx - 1; i++)
        Buffer2[i] = Matrix_Old[i][ny-4];
    // send two cols to west
    dest = rank_e;
    MPI_Send(Buffer1, nx, MPI_CHAR, dest, 0, MPI_COMM_WORLD);
    MPI_Send(Buffer2, nx, MPI_CHAR, dest, 0, MPI_COMM_WORLD);

    // send north
    // Buffer1 = ROW 2 because of ghost layers
    for (int i = 0; i < ny - 1; i++)
        Buffer1[i] = Matrix_Old[2][i];
    // Buffer2 = ROW 3 because of ghost layers
    for (int i = 0; i < ny - 1; i++)
        Buffer2[i] = Matrix_Old[3][i];
    // send two cols to west
    dest = rank_n;
    MPI_Send(Buffer1, ny, MPI_CHAR, dest, 0, MPI_COMM_WORLD);
    MPI_Send(Buffer2, ny, MPI_CHAR, dest, 0, MPI_COMM_WORLD);
    
    // cout << "rank "<<rang<<" has neighbors "<<rank_e<<" "<<rank_w<<" "<<rank_n<<" "<<rank_s<<"\n";
    // cout << "rank "<<rang<<" has diagonal neighbors "<<rank_se<<" "<<rank_sw<<" "<<rank_ne<<" "<<rank_nw<<"\n";
    
    // recieve the data from neighbours

    // recieve west
    source = rank_w;
    MPI_Recv(Buffer1, nx, MPI_CHAR, source, 0, MPI_COMM_WORLD, &status1);
    MPI_Recv(Buffer2, nx, MPI_CHAR, source, 0, MPI_COMM_WORLD, &status2);
    // Buffer1 = Col 0 it is a ghost layers
    for (int i = 0; i < nx - 1; i++)
        Matrix_Old[i][0] = Buffer1[i];
    // Buffer2 = Col 1 it is a ghost layers
    for (int i = 0; i < nx - 1; i++)
        Matrix_Old[i][1] = Buffer2[i];
    cout<<"recieved  west from "<< source<<" to "<<rang<<"\n";

    // recieve south
    source = rank_s;
    MPI_Recv(Buffer1, ny, MPI_CHAR, dest, 0, MPI_COMM_WORLD, &status1);
    MPI_Recv(Buffer2, ny, MPI_CHAR, dest, 0, MPI_COMM_WORLD, &status2);
    // Buffer1 = ROW nx-1 it is a ghost layers
    cout<<"inainte de for 192\n MAXDIM="<<MAX_DIM<<"\n";
    for (int i = 0; i < ny - 1; i++){
        cout<<"\t recieve south\t"<<i<<" "<<nx-1<<"\n";
        Matrix_Old[nx-1][i] = Buffer1[i];
    }
    cout<<"fupa de for 197\n";
    // Buffer2 = ROW nx-2 it is a ghost layers
    for (int i = 0; i < ny - 1; i++)
        Matrix_Old[nx-2][i] = Buffer2[i];
    
    // recieve east
    source = rank_e;
    MPI_Recv(Buffer1, nx, MPI_CHAR, source, 0, MPI_COMM_WORLD, &status1);
    MPI_Recv(Buffer2, nx, MPI_CHAR, source, 0, MPI_COMM_WORLD, &status2);
    // Buffer1 = Col ny-1 it is a ghost layers
    for (int i = 0; i < nx - 1; i++)
        Matrix_Old[i][ny-1] = Buffer1[i];
    // Buffer2 = Col ny-2 it is a ghost layers
    for (int i = 0; i < nx - 1; i++)
        Matrix_Old[i][ny-2] = Buffer2[i];

    // recieve north
    source = rank_n;
    MPI_Recv(Buffer1, ny, MPI_CHAR, source, 0, MPI_COMM_WORLD, &status1);
    MPI_Recv(Buffer2, ny, MPI_CHAR, source, 0, MPI_COMM_WORLD, &status2);
    // Buffer1 = ROW 0 it is a ghost layers
    for (int i = 0; i < ny - 1; i++)
        Matrix_Old[0][i] = Buffer1[i];
    // Buffer2 = ROW 1 because of ghost layers
    for (int i = 0; i < ny - 1; i++)
        Matrix_Old[1][i] = Buffer2[i];
    //send corners
    
    //send SW corner
    coord_corner = {coord[0]-1, coord[1]-1};
    if(coord_corner[0] < 0) coord_corner[0]+=r;
    else if (coord_corner[0] > r) coord_corner[0]-=r;
    if(coord_corner[1] < 0) coord_corner[1]+=c;
    else if (coord_corner[1] > c) coord_corner[1]-=c;
    cout<<"coord_corner "<<coord_corner[0]<<" "<<coord_corner[1]<<"\n";
    MPI_Cart_rank(comm_2d, std::data(coord_corner), &rank_sw);
    dest = rank_sw;
    cout<<"preparing to sent south west from "<< rang<< "to"<<dest<<"\n";
    Buffer1[0] = Matrix_Old[nx-4][2];
    Buffer1[1] = Matrix_Old[nx-4][3]; 
    Buffer1[2] = Matrix_Old[nx-3][2];
    Buffer1[3] = Matrix_Old[nx-3][3];
    MPI_Send(Buffer1, 4, MPI_CHAR, dest, 0, MPI_COMM_WORLD);

    //send SE corner
    coord_corner = {coord[0]-1, coord[1]+1};
    if(coord_corner[0] < 0) coord_corner[0]+=r;
    else if (coord_corner[0] > r) coord_corner[0]-=r;
    if(coord_corner[1] < 0) coord_corner[1]+=c;
    else if (coord_corner[1] > c) coord_corner[1]-=c;
    MPI_Cart_rank(comm_2d, std::data(coord_corner), &rank_se);
    dest = rank_se;
    Buffer1[0] = Matrix_Old[nx-4][ny-4];
    Buffer1[1] = Matrix_Old[nx-4][ny-3]; 
    Buffer1[2] = Matrix_Old[nx-3][ny-4];
    Buffer1[3] = Matrix_Old[nx-3][ny-3];
    MPI_Send(Buffer1, 4, MPI_CHAR, dest, 0, MPI_COMM_WORLD);

    //send NE corner
    coord_corner = {coord[0]+1, coord[1]+1};
    if(coord_corner[0] < 0) coord_corner[0]+=r;
    else if (coord_corner[0] > r) coord_corner[0]-=r;
    if(coord_corner[1] < 0) coord_corner[1]+=c;
    else if (coord_corner[1] > c) coord_corner[1]-=c;
    MPI_Cart_rank(comm_2d, std::data(coord_corner), &rank_ne);
    dest = rank_ne;
    Buffer1[0] = Matrix_Old[2][ny-4];
    Buffer1[1] = Matrix_Old[2][ny-3]; 
    Buffer1[2] = Matrix_Old[3][ny-4];
    Buffer1[3] = Matrix_Old[3][ny-3];
    MPI_Send(Buffer1, 4, MPI_CHAR, dest, 0, MPI_COMM_WORLD);

    //send NW corner
    coord_corner = {coord[0]+1, coord[1]-1};
    if(coord_corner[0] < 0) coord_corner[0]+=r;
    else if (coord_corner[0] > r) coord_corner[0]-=r;
    if(coord_corner[1] < 0) coord_corner[1]+=c;
    else if (coord_corner[1] > c) coord_corner[1]-=c;
    MPI_Cart_rank(comm_2d, std::data(coord_corner), &rank_nw);
    dest = rank_nw;
    Buffer1[0] = Matrix_Old[2][2];
    Buffer1[1] = Matrix_Old[2][3]; 
    Buffer1[2] = Matrix_Old[3][2];
    Buffer1[3] = Matrix_Old[3][3];
    MPI_Send(Buffer1, 4, MPI_CHAR, dest, 0, MPI_COMM_WORLD);
    
    //recieve corners
    
    //recieve SW corner
    source = rank_sw;
    MPI_Recv(Buffer1, 4, MPI_CHAR, source, 0, MPI_COMM_WORLD, &status1);
    Matrix_Old[nx-1][0] = Buffer1[0];
    Matrix_Old[nx-1][1] = Buffer1[1]; 
    Matrix_Old[nx-2][0] = Buffer1[2];
    Matrix_Old[nx-2][1] = Buffer1[3];

    //recieve SE corner
    source = rank_se;
    MPI_Recv(Buffer1, 4, MPI_CHAR, source, 0, MPI_COMM_WORLD, &status1);
    Matrix_Old[nx-1][ny-1] = Buffer1[0];
    Matrix_Old[nx-1][ny-2] = Buffer1[0]; 
    Matrix_Old[nx-2][ny-1] = Buffer1[0];
    Matrix_Old[nx-2][ny-2] = Buffer1[0];

    //recieve NE corner
    source = rank_ne;
    MPI_Recv(Buffer1, 4, MPI_CHAR, source, 0, MPI_COMM_WORLD, &status1);
    Matrix_Old[0][ny-2] = Buffer1[0];
    Matrix_Old[0][ny-1] = Buffer1[0]; 
    Matrix_Old[1][ny-2] = Buffer1[0];
    Matrix_Old[1][ny-1] = Buffer1[0];

    //recieve NW corner
    source = rank_nw;
    MPI_Recv(Buffer1, 4, MPI_CHAR, source, 0, MPI_COMM_WORLD, &status1);
    Matrix_Old[0][0] = Buffer1[0];
    Matrix_Old[0][1] = Buffer1[1]; 
    Matrix_Old[1][0] = Buffer1[2];
    Matrix_Old[1][1] = Buffer1[3];

    // clear
    free(Buffer1);
    free(Buffer2);

}

void Iterate(int nx, int ny, Matrix &Matrix_Old, Matrix &Matrix_New, int generations, int rang, array<int, 2> coord, MPI_Comm comm_2d, int r, int c){
    for (int generation = 1; generation <= generations; generation++){
        //we first communicate to update the ghost layers for each subdomain after the last iteration
        //also since the ghost layers originally haven't been set in the first iteration we set them 
        CommunicateBorders (nx,  ny,  Matrix_Old, Matrix_New, rang, coord, comm_2d, r, c);
        cout << "rank "<<rang<<" communicated\n";
        UpdateMatrix (nx,  ny,  Matrix_Old, Matrix_New);
        cout << "rank "<<rang<<" updated\n";
        Matrix_Old = Matrix_New;
    }
}

//an input generator for at least one type of initial configurations
void  InputGenerator(int nx, int ny, Matrix &Matrix_Generated, Matrix &Matrix_New, float chance = 0.28){
    // we denote the dead cells with char '0' and the alive cells with char '1'
    // using char as it is a byte datatype for the communication part down the road
    // we have to leave empty the ghost layers and later communicate them with their neighbours
    for (int i = 2; i < nx-2; i++){
        for (int j = 2; j < ny-2; j++){
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

int paralel(program_options::Options opts, int rank, int numprocesses){
    // let openMPI take care of splitting into suitable subdomain grid
    int ndims = 2; //2D grid
    std::array<int, 2> dims = {0, 0};
    MPI_Dims_create(numprocesses, ndims, std::data(dims));

    int r = dims[1]; //number or rows the domain was split
    int c = dims[0]; //number or cols the domain was split


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

    //row_ndom are the number of rows and columns per slice
    //row_bdom is the number of slices that have +1 rows and columns
    int row_ndom, row_bdom, col_ndom, col_bdom;
    row_ndom = opts.N/r;
    row_bdom = opts.N%r;
    col_ndom = opts.M/c;
    col_bdom = opts.M%c;

    //coordinates of rank in virtual topology
    array<int, 2> coord = {-1, -1};    
    MPI_Cart_coords(comm_2d, rank, ndims, std::data(coord));

    // number of columns and rows in this subdomain (with added ghost layers)
    std::tuple <int, int> rw_col = GetSubdomainDims(coord, row_ndom, row_bdom, col_ndom, col_bdom, r, c);
    int N_col, N_row;
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
    cout<<"nrow= "<<N_row<<" ncol= "<<N_col<<"\n";
    for(size_t i=0;i<opts.N;++i){
        Matrix_New[i].resize(N_col);
        Matrix_Old[i].resize(N_col);
    }
    InputGenerator(N_row, N_col, Matrix_Old, Matrix_New);
    cout<<"Generation of the domain "<< rank <<" successfull------------------------------------------------" << std::endl;
    
    Iterate(N_row, N_col, Matrix_Old, Matrix_New, opts.iters, rank, coord, comm_2d,r,c);
    cout<<"Iteration of the domain "<< rank <<" successfull------------------------------------------------" << std::endl;
    // PrintOutput(opts.N, opts.M, Matrix_New);

    return 0;
}