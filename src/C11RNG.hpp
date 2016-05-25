// binomial_distribution
#include <random>

class C11RNG{

public:

    std::mt19937 gen;
    C11RNG(){

        std::random_device std_rand;
        gen = std::mt19937(std_rand());

    };

    int draw_binomial(int n, double p){
        std::binomial_distribution<int> distribution(n, p);
        return distribution(gen);
    };

    arma::vec draw_binomial(int n, double p, int size){
        std::binomial_distribution<int> distribution(n, p);

        arma::vec result(size);
        for (int i = 0; i < size; ++i)
        {
            result(i) = (double) distribution(gen);
        }
        return result;
    };
};