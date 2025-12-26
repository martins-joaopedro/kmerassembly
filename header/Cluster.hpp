#ifndef CLUSTER_HPP
#define CLUSTER_HPP

#include <bits/stdc++.h>
using namespace std;

class Cluster {
    
        int MAX_CLUSTERS;
        int EXPLORATION;
        int THRESHOLD;
        int K;
        bool initialized = false;
        map<string, double> centroid;
        int size = 0;
        double norm = 0.0;

    public:
        Cluster(int MAX_CLUSTERS, int EXPLORATION, int THRESHOLD, int K);
        void update(const map<string, double>& readVec); 
        void initialize(const map<string, double>& readVec);
        double cosineSimilarity(const map<string, double>& readVec, const Cluster& c);
        map<string, double> normalizeKmerFreq(const map<string, int>& kmerCounts);
        void clusterize();

};

#endif //CLUSTER_HPP