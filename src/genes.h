//
// Created by Oliver Feighan on 8/31/20.
//

#ifndef GAUS_GENES_H
#define GAUS_GENES_H

#include <armadillo>

namespace gaus::genetic_algorithm {

class Cell {
private:
  arma::umat
  make_genes(arma::SizeMat gene_size);

public:
  int tag; //an identifer for the cell
  arma::umat genes; //parameters that define behaviour
  double fitness; //to be populated by the fitness function

  void
  mutate(double mutation_rate);

  Cell(int t, arma::SizeMat gene_size, int f) {
    tag = t;
    genes = make_genes(gene_size);
    fitness = f;
  }

};

class Colony {
private:

  std::vector<Cell>
  make_colony(int n_cells, arma::SizeMat gene_size);

public:
  std::vector<Cell> cells;
  int colony_size;
  double mutation_rate;
  double current_top_fitness = NAN;

  Colony(int size, arma::SizeMat gene_size, double m_rate) {
    cells = make_colony(size, gene_size);
    colony_size = size;
    assert(m_rate < 1);
    mutation_rate = m_rate;
  }

  void
  find_fitnesses(const std::function<double(arma::umat)> & fitness_function);

  void
  sort_fitnesses();

  void
  mutate_colony();

  void
  make_next_generation(double selection);

};


}

#endif //GAUS_GENES_H
