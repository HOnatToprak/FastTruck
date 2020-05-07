#include "TSP.h"
#include <chrono>
#include <fstream>
#include <vector>
#include "armadillo"


int main() {
    using namespace arma;
    using namespace FastTruck::TSP;

    srand(time(0));

    fmat distances;
    distances.load("resources/att48_d.txt");


    TSPInit init{
        distances,
        1000000,
        100,
        0.01f,
        0.75f,
        0.05f,
    };

    solve(init);

    return 0;
}

