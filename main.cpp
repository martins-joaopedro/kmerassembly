#include <bits/stdc++.h>
#include <random>
#include "header/GeneticAlgorithm.hpp"
#include "header/Cluster.hpp"

using namespace std;

int main() {
    
    int K = 6;
    int EXPLORATION = 150;
    const int MAX_CLUSTERS = 5;
    const double THRESHOLD = 0.2;  

    Cluster cluster(MAX_CLUSTERS, EXPLORATION, THRESHOLD, K);
    //cluster.clusterize();
    
    int chromosomesNumber = 2;
    int genesNumber = 2;
    GeneticAlgoritm GA(chromosomesNumber, genesNumber, MAX_CLUSTERS);
    GA.start();

    return 0;
}