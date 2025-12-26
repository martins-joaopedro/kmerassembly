#ifndef READER_HPP
#define READER_HPP

#include <bits/stdc++.h>
using namespace std;

struct ReadIndex {
    int id;
    // local de início da sequencia
    uint64_t offset; 
    uint32_t length; 
};

class Reader {
    
    public:
        Reader();
        static vector<ReadIndex> indexFasta(const string& filename);
        static string loadReadAt(ifstream& fasta, uint64_t offset, uint32_t len);

};

#endif //READER_HPP