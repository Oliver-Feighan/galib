//
// Created by Oliver Feighan on 8/31/20.
//

#ifndef GAUS_GENETIC_ALGORITHM_H
#define GAUS_GENETIC_ALGORITHM_H

#include "genetic_algorithm/genes.h"

namespace gaus::genetic_algorithm {

struct Solution {

  arma::umat genes;
  int n_generations;

};

Solution
find_solution(const std::function<double(arma::umat)> & fitness_function,
               int colony_size,
               arma::SizeMat gene_length,
               double selection_rate,
               double mutation_rate);

}

#endif //GAUS_GENETIC_ALGORITHM_H
