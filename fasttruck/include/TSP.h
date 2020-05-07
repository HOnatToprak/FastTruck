//
// Created by onat on 1.05.2020.
//

#ifndef FASTTRUCK_TSP_H
#define FASTTRUCK_TSP_H

#include <iostream>
#include <string>
#include "armadillo"

namespace FastTruck
{
namespace TSP
{
    using namespace arma;

    struct TSPInit {
        fmat distances;
        int last_generation = 0;
        int population_size = 0;
        float elite_rate = 0.0f;
        float crossover_rate = 0.0f;
        float mutation_rate = 0.0f;
    };

    struct TSPResult {
        fvec best_route;
        double best_route_cost;
    };

    TSPResult solve(const TSPInit initializers);
    fmat generate_random_population(const int population_size,const  int number_of_genes);
    vec create_fitness_matrix(const fmat& population,const fmat& distances);
    inline double calculate_fitness(const subview_row<float>& route, const fmat& distances);
    double calculate_route_cost(const subview_row<float>&, const fmat& distances);
    fmat select_parents_roulette_wheel(const vec& fitnesses);
    fmat crossover_parents_ox(const fmat& population, const fmat& parent_indices, const float crossover_rate, const int elite_size);
    void apply_mutations(fmat& population, float mutation_rate);

} // namespace TSP
} // namespace FastTruck
#endif //FASTTRUCK_TSP_H
