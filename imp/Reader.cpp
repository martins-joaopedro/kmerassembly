#include "../header/Reader.hpp"

Reader::Reader() {}

string Reader::loadReadAt(ifstream& fasta, uint64_t offset, uint32_t len) {
    fasta.clear();
    fasta.seekg(offset, ios::beg);

    string seq, line;
    seq.reserve(len);

    while (getline(fasta, line)) {
        if (!line.empty() && line.back() == '\r')
            line.pop_back();

        // para quando encontrar o proximo read
        if (!line.empty() && line[0] == '>')
            break;  

        // TODO: necessário??
        for (char c : line) {
            if (seq.size() == len)
                return seq;
            seq.push_back(c);
        }
    }

    return seq;
}

vector<ReadIndex> Reader::indexFasta(const string& filename) {
    
    ifstream file(filename, ios::binary);
    if (!file.is_open()) 
        throw runtime_error("Erro ao abrir FASTA");

    vector<ReadIndex> index;
    string line;
    int id = 0;

    uint64_t seqOffset = 0;
    uint32_t seqLen = 0;
    bool readingSeq = false;

    while (true) {
        uint64_t lineStart = file.tellg();
        if (!getline(file, line))
            break;

        // tira o CRLF
        if (!line.empty() && line.back() == '\r')
            line.pop_back();

        if (!line.empty() && line[0] == '>') {
            // processa read anterior
            if (readingSeq) {
                index.push_back({id++, seqOffset, seqLen});
                seqLen = 0;
            }

            // offset da começa dps do header, depois do '\n'
            seqOffset = file.tellg(); 
            readingSeq = true;
        } else if (readingSeq) {
            seqLen += line.size();
        }
    }

    // ultima
    if (readingSeq) {
        index.push_back({id++, seqOffset, seqLen});
    }

    file.close();
    return index;
}