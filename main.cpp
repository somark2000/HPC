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
  //MPI start
  MPI_Init(&argc, &argv);

  // parse args
  auto opts = program_options::parse(argc, argv);
  
  //pre generate the input for all the runs for verification
  if(opts.verification == 1) sequential_solver::GlobalInputGenerator(opts);
  
  //define each process their rank and the total number of processes
  int rank, num_processes;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &num_processes);
  if(opts.mpi_mode == 1){
    if(rank == 0)
        std::printf("Sequential solver------------------------------------------------------------\n");
    sequential_solver::sequential(opts);
  }
  else if (opts.mpi_mode == 2/* && !check_prime(num_processes)*/){
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
    cout<<"Solver terminated\n";
  }
  MPI_Finalize();
}
catch (std::exception &e) {
  std::cout << e.what() << std::endl;
  return EXIT_FAILURE;
}