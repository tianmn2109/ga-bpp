#include "ga.hh"
#include <assert.h>
#include <iomanip>

int main(int argc, char **argv)
{
    srand((unsigned)time(NULL));
    // initialization
    struct bpp b;
    init(b, argv[1], 400);

    // baseline
    cout << "Baseline:" << endl;
    b.nb = firstFit(b) + 100; 
    cout << b.nb - 100 << endl;
    cout << "#Objs: " << b.n << endl;
    cout << "#Bins: " << b.nb << endl;
    
    // printWeight(b);
    // printCapacity(b);

    // genetic algorithm
    cout << "Genetic Algorithm:" << endl;

    // init
    struct phenotype zero(b.nb, b.n, b.c);
    for (int i = 0; i < 6; i++) {
        deepcopy(tmp[i], zero, b);
    }
    
    initialize(population, b);
    statistics(population);

    // iterate
    int i, j;
    int cnt = 1000;
    while (cnt-- > 0) {
        if (cnt % 10 == 0) {
            cout << setw(5) << flush << '\r' << "#" << cnt;
        }
        i = select();
        j = select();
        // cout << i << " " << j << endl;
        deepcopy(tmp[0], population[i], b);
        while (j == i)
            j = select();
        deepcopy(tmp[1], population[j], b);
        crossOver(tmp[0], tmp[1], tmp[2], b);
        crossOver(tmp[1], tmp[0], tmp[3], b);
        mutate(tmp[0], tmp[4], b);
        mutate(tmp[1], tmp[5], b);

        // check the validity of tmp
        assert(isValidSolution(tmp[0], b) == 1);
        assert(isValidSolution(tmp[1], b) == 1);
        assert(isValidSolution(tmp[2], b) == 1);
        assert(isValidSolution(tmp[3], b) == 1);
        assert(isValidSolution(tmp[4], b) == 1);
        assert(isValidSolution(tmp[5], b) == 1);
        
        // select 2 best phenotype and add it to the population
        int best_so_far = b.nb;
        int best_so_far_id;
        for (int k = 0; k < 6; k++) {
            if (tmp[k].num < best_so_far) {
                best_so_far = tmp[k].num;
                best_so_far_id = k;
            }
        }
        deepcopy(population[i], tmp[best_so_far_id], b);
        best_so_far = b.nb;
        tmp[best_so_far_id].num = b.nb;

        for (int k = 0; k < 6; k++) {
            if (tmp[k].num < best_so_far) {
                best_so_far = tmp[k].num;
                best_so_far_id = k;
            }
        }
        deepcopy(population[j], tmp[best_so_far_id], b);
        // statistics(population);
    }

    cout << endl << "After iterations: " << endl;
    statistics(population);
    // printCapacity(population[0], b);
    // calcCap(population[0], b);

    delete[] b.w;
    delete[] b.cap;
    for (int i = 0; i < b.nb; i++) {
        delete[] b.alloc[i];
    }
    return 0;
}
