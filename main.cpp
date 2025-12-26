#include <bits/stdc++.h>
#include <random>
#include <omp.h>
#include "header/GeneticAlgorithm.hpp"
#include "header/Cluster.hpp"

using namespace std;

int main() {
    
    omp_set_num_threads(16);

    int K = 6;
    int EXPLORATION = 150;
    const int MAX_CLUSTERS = 5;
    const double THRESHOLD = 0.2;  
    const int MIN_OVERLAPING = 50;

    Cluster cluster(MAX_CLUSTERS, EXPLORATION, THRESHOLD, K);
    //cluster.clusterize();
    
    int chromosomesNumber = 20;
    int genesNumber = 200;
    GeneticAlgoritm GA(chromosomesNumber, genesNumber, MAX_CLUSTERS, MIN_OVERLAPING);
    GA.start();

    return 0;
}