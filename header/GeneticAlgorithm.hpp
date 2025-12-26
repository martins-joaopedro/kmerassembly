#ifndef GA_HPP
#define GA_HPP

#include "Chromosome.hpp"
#include <bits/stdc++.h>
using namespace std;

class GeneticAlgoritm {

    public:
        // clusters
        int islandsNumber;
        int genesNumber;
        int chromosomesNumber;
        int MIN_OVERLAPING;
        int eras = 0;
        vector<uint8_t> usedReads;

        GeneticAlgoritm(int chromosomesNumber, int genesNumber, int islandsNumber, int MIN_OVERLAPING);
        vector<Chromosome> createPopulation(int islandNumber, int era);
        void endMiscenegationPhase(vector<Chromosome> &pop, int islandNumber, int era);
        vector<Chromosome> crossover(Chromosome c1, Chromosome c2);
        void start();

};


#endif //GA_HPP