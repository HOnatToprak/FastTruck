//
// Created by onat on 1.05.2020.
//

#include <iostream>
#include <TSP.h>
#include <unordered_map>
#include "Utils.h"

namespace FastTruck
{
namespace TSP
{

    using namespace FastTruck::Utils;
    TSPResult solve(const TSPInit init) {
        timeitinit;
        int elite_size = (int)(init.population_size * init.elite_rate);

        // Create initial population
        int number_of_genes = init.distances.n_cols;
        fmat population = generate_random_population(init.population_size, number_of_genes);
        vec fitnesses = create_fitness_matrix(population, init.distances);

        // Sort population by fitnesses
        uvec sorted_indices = sort_index(fitnesses, "descend");
        fitnesses = reorder_with_indices(fitnesses, sorted_indices);
        population = reorder_with_indices(population, sorted_indices);

        for(int gen = 0; gen < init.last_generation; ++gen){
            fmat parent_indices = select_parents_roulette_wheel(fitnesses);
            fmat children = crossover_parents_ox(population, parent_indices, init.crossover_rate, elite_size);
            apply_mutations(children, init.mutation_rate);
            vec children_fitnesses = create_fitness_matrix(children, init.distances);

            fitnesses.resize(init.population_size + elite_size);
            population.resize(init.population_size + elite_size, population.n_cols);
            population.tail_rows(init.population_size) = children;

            fitnesses.tail_rows(init.population_size) = children_fitnesses;

            sorted_indices = sort_index(fitnesses, "descend");
            fitnesses = reorder_with_indices(fitnesses, sorted_indices);
            population = reorder_with_indices(population, sorted_indices);

            fitnesses.resize(init.population_size);
            population.resize(init.population_size, population.n_cols);

            std::cout << "Generation " << gen + 1 << " => " <<
            calculate_route_cost(population.row(0), init.distances)
            << "\n";

        }

        return TSPResult{
            population.row(0),
            calculate_route_cost(population.row(0), init.distances),
        };
    }

    fmat generate_random_population(const int population_size,const int number_of_genes) {
        fmat m(number_of_genes, population_size);

        for(int i = 0; i < number_of_genes; ++i){
            m.row(i).fill(i);
        }

        for(int i = 0; i < population_size; ++i){
            //TODO: Refactor to std::random
            //TODO: Try armadillo shuffle on subviewcols
            std::random_shuffle(m.begin_col(i), m.end_col(i));
        }

        return m.t();
    }

    vec create_fitness_matrix(const fmat& population, const fmat& distances) {
        vec fitnesses(population.n_rows);
        for(int i = 0; i < population.n_rows; ++i)
            fitnesses(i) = calculate_fitness(population.row(i), distances);
        return fitnesses;
    }

    double calculate_route_cost(const subview_row<float>& route, const fmat& distances) {
        double sum = 0;
        for(int i = 0; i < route.n_cols - 1; ++i){
            sum += distances(route(i), route(i+1));
        }
        sum += distances(route.back(), route(0));
        return sum;
    }

    double calculate_fitness(const subview_row<float>& route, const fmat& distances) {
        return 1 / calculate_route_cost(route, distances);
    }

    fmat select_parents_roulette_wheel(const vec& fitnesses) {
        // Random generator with weights as fitnesses
        static std::mt19937_64 gen(rand());
        std::discrete_distribution<> d(fitnesses.begin(), fitnesses.end());

        // n * 2 parens created and coupled
        fmat parents(fitnesses.n_rows, 2);
        auto fill_parent = [&d](float& p){ p = d(gen);};
        parents.for_each(fill_parent);

        return parents;
    }

    fmat crossover_parents_ox(const fmat& population, const fmat& parent_indices, const float crossover_rate, const int elite_size) {
        fmat children(population.n_rows, population.n_cols);

        // For selecting cutting points
        std::mt19937_64 gen(rand());
        std::uniform_int_distribution<> di(0, population.n_cols - 1);

        //crossover_rate
        std::uniform_real_distribution<float> df;

        // For checking collisions;
        //std::unordered_map<float, bool> check(population.n_cols);
        bool* checks = (bool*) alloca(sizeof(bool) * population.n_cols);
        memset(checks, 0, sizeof(bool) * population.n_cols);

        const int gene_size = population.n_cols;
        for(int i = 0; i < population.n_rows; ++i){
            const auto parent1 = population.row(parent_indices(i, 0));
            const auto parent2 = population.row(parent_indices(i, 1));
            auto child = children.row(i);
            if(parent_indices(i,0) < elite_size  || df(gen) < crossover_rate){


                int cut1 = di(gen), cut2 = di(gen);
                if(cut1 > cut2) std::swap(cut1, cut2);

                for(int j = cut1; j <= cut2; ++j)
                    checks[(int)parent1(j)] = true;

                // Reduce parent2
                frowvec new_parent_2(gene_size - cut2 + cut1 - 1);
                for(int j = 0, k= 0; j < gene_size; ++j) {
                    if (!checks[(int)parent2(j)]) {
                        new_parent_2(k++) = parent2(j);
                    }
                }

                int iter = 0;
                for(int j = 0; j < cut1; ++j)
                    child(j) = new_parent_2(iter++);

                for(int j = cut1; j <= cut2; ++j)
                    child(j) = parent1(j);

                for(int j = cut2 + 1; j < gene_size; ++j)
                    child(j) = new_parent_2(iter++);

                //checks.clear();
                memset(checks, 0, sizeof(bool) * population.n_cols);
            }
            else
            {
                child = parent1;
            }
        }
        return children;
    }

    void apply_mutations(fmat& population, float mutation_rate) {
        static std::mt19937_64 gen(rand());
        std::uniform_int_distribution<> dint(1, population.n_cols - 2);
        std::uniform_real_distribution<> df;

        for(int i = 0; i < population.n_rows; ++i){
            if(df(gen) < mutation_rate){
                int cut1 = dint(gen), cut2 = dint(gen);
                if(cut1 > cut2) std::swap(cut1, cut2);
                population.row(i).cols(cut1, cut2) = reverse(population.row(i).cols(cut1, cut2));
            }
        }
    }

    frowvec crossover_ox(const subview_row<float>& parent1, const subview_row<float>& parent2){
        // For selecting cutting points
        std::mt19937_64 gen(rand());
        std::uniform_int_distribution<> di(1, parent1.n_cols - 1);

        // For checking collisions;
        std::unordered_map<float, bool> checks(parent1.n_cols);

        int cut1 = di(gen), cut2 = di(gen);
        if(cut1 > cut2) std::swap(cut1, cut2);

        for(int i = cut1; i < cut2; ++i)
            checks[parent1(i)] = true;

        // Reduce parent2
        frowvec new_parent_2(parent1.n_cols - cut2 + cut1 - 1);
        for(int i = 0, j = 0; i < parent1.n_cols; ++i) {
            if (!checks[parent2(i)])
                new_parent_2(j++) = parent2(i);
        }

        frowvec child(parent1.n_cols);
        int iter = 0;
        for(int i = 0; i < cut1; ++i)
            child(i) = new_parent_2(iter++);

        for(int i = cut1; i <= cut2; ++i)
            child(i) = parent1(i);

        for(int i = cut2 + 1; i < parent1.n_cols; ++i)
            child(i) = new_parent_2(iter++);

        return child;
    }




} // namespace TSP
} // namespace FastTruck