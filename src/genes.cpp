#include "genes.h"

namespace gaus::genetic_algorithm {

arma::umat
Cell::make_genes(const arma::SizeMat gene_size) {

  return
      arma::randi<arma::umat>(gene_size, arma::distr_param(0, 1));

}

void
Cell::mutate(const double mutation_rate) {

  const int n_mutations = std::round(mutation_rate * genes.n_elem);

  const auto random_indices =
      arma::randi<arma::uvec>(n_mutations,
                              arma::distr_param(0, genes.n_elem - 1));

  genes(random_indices) =
      arma::randi<arma::uvec>(n_mutations, arma::distr_param(0, 1));

}

std::vector<Cell>
Colony::make_colony(const int n_cells,
                    const arma::SizeMat gene_size) {

  std::vector<Cell> initial_cells;

  for (int i = 0; i < n_cells; i++) {
    initial_cells.push_back(Cell(i, gene_size, NAN));

  }

  return initial_cells;
}

void
Colony::find_fitnesses(
    const std::function<double(arma::umat)> & fitness_function) {
  for (auto & cell : cells) {
    cell.fitness = fitness_function(cell.genes);

  }
}

void
Colony::sort_fitnesses() {
  std::sort(cells.begin(), cells.end(), [](Cell & a, Cell & b) {
    return a.fitness < b.fitness;

  });

  current_top_fitness = cells[0].fitness;
}

void
Colony::mutate_colony() {
  for (auto & cell : cells) {
    cell.mutate(mutation_rate);

  }
}

void
Colony::make_next_generation(const double selection){
  const int n_survivors = std::round(selection * colony_size);

  //make survivors
  auto start = cells.begin();
  auto end = cells.begin() + n_survivors;

  const std::vector<Cell> survivors(start, end);

  std::vector<arma::umat> genetic_material = {};

  for(int i = 0; i < n_survivors; i++){

    const int gene_length = survivors[i].genes.n_rows;

    const auto top_half = survivors[i].genes.head_rows(gene_length / 2);
    const auto bottom_half = survivors[i].genes.tail_rows(gene_length / 2);

    genetic_material.push_back(top_half);
    genetic_material.push_back(bottom_half);

  }

  for(int i = 0; i < colony_size; i++){

    const auto a =
        arma::randi<int>(arma::distr_param(0, n_survivors * 2 - 1));

    const auto b =
        arma::randi<int>(arma::distr_param(0, n_survivors * 2 - 1));

    const auto gen_a = genetic_material[a];
    const auto gen_b = genetic_material[b];

    const arma::umat new_genes =
        arma::join_vert(genetic_material[a], genetic_material[b]);

    cells[i].genes = new_genes;
  }

}

}