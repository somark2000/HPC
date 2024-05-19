#ifndef OPTIONS_H
#define OPTIONS_H

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


namespace program_options {

  struct Options {
    unsigned int mpi_mode;
    std::string name;
    size_t N;
    size_t M;
    size_t iters;
    // unsigned int reorder;
    // unsigned int scaling;
    // float chance;
    // unsigned int verification;

    void print() const {
      std::printf("mpi_mode: %u\n", mpi_mode);    
      std::printf("name: %s\n", name.c_str());
      std::printf("N: %zu\n", N);
      std::printf("M: %zu\n", M);
      std::printf("iters: %zu\n", iters);
    //   std::printf("reoder: %u\n", reorder);
    //   std::printf("scaling: %u\n", scaling);
    //   std::printf("chance of alive: %u\n", chance);
    //   std::printf("verification: %u\n", verification);
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
    else if( std::string(argv[1]) == std::string("parn"))
      opts.mpi_mode = 3;
    else
    throw std::runtime_error("invalid parameter for mpi_mode (valid are 'seq' and 'par')");
    opts.name = argv[2];
    if (std::sscanf(argv[3], "%zu", &opts.N) != 1 && opts.N >= 2)
      throw std::runtime_error("invalid parameter for N");
    if (std::sscanf(argv[4], "%zu", &opts.M) != 1 && opts.M >= 2)
      throw std::runtime_error("invalid parameter for M");
    if (std::sscanf(argv[5], "%zu", &opts.iters) != 1 && opts.iters != 0)
      throw std::runtime_error("invalid parameter for iters");
    // if (std::sscanf(argv[6], "%u", &opts.reorder) != 1 && opts.reorder != 0)
    //     throw std::runtime_error("invalid parameter for reorder flag (valid are 0 and 1) and recieved");
    // if (std::string(argv[7]) == std::string("off"))
    //   opts.scaling = 0;
    // else if( std::string(argv[7]) == std::string("strong"))
    //   opts.scaling = 1;
    // else if( std::string(argv[7]) == std::string("weak"))
    //   opts.scaling = 2;
    // else
    // throw std::runtime_error("invalid parameter for scaling (valid are 'off', 'strong' and 'weak')");
    // if (std::sscanf(argv[8], "%zu", &opts.chance) > 1 && opts.chance < 0)
    //   throw std::runtime_error("invalid parameter for chance of alive");
    // if (std::string(argv[9]) == std::string("off"))
    //   opts.verification = 0;
    // else if( std::string(argv[9]) == std::string("on"))
    //   opts.verification = 1;
    // else
    // throw std::runtime_error("invalid parameter for verification (valid are 'off' and 'on')");
    return opts;
  }

} // namespace program_options

#endif
