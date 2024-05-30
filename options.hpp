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
    int reorder;
    int scaling;
    unsigned int scaling_config;
    float chance;
    int verification;

    void print() const {
      std::printf("mpi_mode: %u\n", mpi_mode);    
      std::printf("name: %s\n", name.c_str());
      std::printf("N: %zu\n", N);
      std::printf("M: %zu\n", M);
      std::printf("iters: %zu\n", iters);
      std::printf("reoder: %u\n", reorder);
      std::printf("scaling: %u\n", scaling);
      std::printf("scaling_config: %u\n", scaling_config);
      std::printf("chance of alive: %f\n", chance);
      std::printf("verification: %u\n", verification);
    }
  };

  auto parse(int argc, char *argv[]) {
    if (argc != 11)
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
    if (std::string(argv[6]) == std::string("false"))
      opts.reorder = 0;
    else if( std::string(argv[6]) == std::string("true"))
      opts.reorder = 1;
    else throw std::runtime_error("invalid parameter for reorder flag (valid are true and false)");
    if (std::string(argv[7]) == std::string("off"))
      opts.scaling = 0;
    else if( std::string(argv[7]) == std::string("strong"))
      opts.scaling = 1;
    else if( std::string(argv[7]) == std::string("weak"))
      opts.scaling = 2;
    else
    throw std::runtime_error("invalid parameter for scaling (valid are 'off', 'strong' and 'weak')");
    if (std::sscanf(argv[8], "%u", &opts.scaling_config) != 1 && opts.scaling_config >= 0 && opts.scaling_config < 5)
      throw std::runtime_error("invalid parameter for scaling_config");
    if (std::sscanf(argv[9], "%f", &opts.chance) > 1 && opts.chance < 0)
      throw std::runtime_error("invalid parameter for chance of alive");
    if (std::string(argv[10]) == std::string("off"))
      opts.verification = 0;
    else if( std::string(argv[10]) == std::string("on"))
      opts.verification = 1;
    else
    throw std::runtime_error("invalid parameter for verification (valid are 'off' and 'on')");
    return opts;
  }

} // namespace program_options

#endif
