#ifndef CHROMOSOME_HPP
#define CHROMOSOME_HPP

#include <bits/stdc++.h>
using namespace std;

class Chromosome {
    public:
        int id;
        int genesNumber;
        int generation;
        int fitness;
        int minOverlap;
        string status;
        bool overlapped;
        
        vector<string> genes;
        vector<int> order;
        vector<vector<int>> overlaps;
        
        Chromosome(int id, int genesNumber, int generation, int minOverlap, string status);
        int computeFitness(vector<int>& order);
        void updateFitness(const vector<int>& ord);
        void mutate();
        void computeOverlapping();
        void insertionMove(vector<int>& order, int i, int j);
        int deltaInsertion(int i, int j);
        vector<string> getFormedContigs();
        void printOverlapMatrix();

};

#endif // CHROMOSOME_HPP