//
// Created by onat on 5.05.2020.
//

#ifndef FASTTRUCK_ROOT_UTILS_H
#define FASTTRUCK_ROOT_UTILS_H

#define timeitinit auto t1 = std::chrono::high_resolution_clock::now();\
                   auto t2 = std::chrono::high_resolution_clock::now()
#define timeit(x) t1 = std::chrono::high_resolution_clock::now();\
                  x;\
                  t2 = std::chrono::high_resolution_clock::now();\
                  std::cout << "Function " #x " : "\
                  << std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count()\
                  << " microseconds\n"

#include "armadillo"
#include "ostream"

namespace FastTruck
{
namespace Utils
{

    using namespace arma;

    template<typename T>
    T reorder_with_indices(const T& m, const uvec& indices) {
        T x(m.n_rows, m.n_cols);
        for(int i = 0; i < indices.n_rows; ++i)
            x.row(i) = m.row(indices(i));
        return x;
    }



    inline void iosp(){std :: cout << "\n";}
} // namespace FastTruck
} // namespace Utils
#endif //FASTTRUCK_ROOT_UTILS_H
