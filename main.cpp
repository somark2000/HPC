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
#include "paralel.hpp"
#include "sequential.hpp"
#include "paralel_neighbor.hpp"
using namespace std;  

// check if n is a prime number greater than 2
bool check_prime(int n){
  for(int i=2; i<n; i++){
    if( n%i == 0){
      return false;
    }
  }
  return true;
}

int main(int argc, char **argv) try{
  //start timer (before MPI Init)
  auto start1 = std::chrono::system_clock::now();

  //MPI start
  MPI_Init(&argc, &argv);

  // parse args
  auto opts = program_options::parse(argc, argv);
  
  //define each process their rank and the total number of processes
  int rank, num_processes;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &num_processes);
  if(opts.mpi_mode == 1){
    if(rank == 0)
        std::printf("Sequential solver------------------------------------------------------------\n");
    sequential_solver::sequential(opts);
  }
  else if (opts.mpi_mode == 2 && !check_prime(num_processes)){
    if(rank == 0){
        opts.print();
        std::printf("2D paralel decomposition with %i subdomains----------------------------------\n", num_processes);
    }
    paralel_solver::paralel(opts, rank, num_processes);
  }
  else if (opts.mpi_mode == 3 && !check_prime(num_processes)){
    if(rank == 0){
        opts.print();
        std::printf("2D paralel-hood decomposition with %i subdomains----------------------------------\n", num_processes);
    }
    paraleln_solver::paralel_neighbor(opts, rank, num_processes);
  } 
  else {
    std::cout<<"decomposition probem not possible!\n";
  }

  if (rank == 0){
    cout<<"Paralel terminated\n";
    auto end1 = std::chrono::system_clock::now();
    auto runtime1 = std::chrono::duration_cast<std::chrono::milliseconds>(end1 - start1); //runtime with init
    std::cout << "The program executed in " << runtime1.count() << " milliseconds." << std::endl;

    // //write specs to .txt file
    // std::ofstream file;
    // file.open("runtime_data.csv", std::ios::app);
    // if (!file.is_open()) {
    //   std::cerr << "File could not be opened!" << std::endl;
    //   return 1;
    // }
    // file << opts.mpi_mode << ", " << opts.N << ", " << num_processes << ", " << runtime1.count() << std::endl;
    // file.close();
  }
  MPI_Finalize();
}
catch (std::exception &e) {
  std::cout << e.what() << std::endl;
  return EXIT_FAILURE;
}