#include "../header/Cluster.hpp"
#include "../header/Reader.hpp"

Cluster::Cluster(int MAX_CLUSTERS, int EXPLORATION, int THRESHOLD, int K) {
    this->MAX_CLUSTERS = MAX_CLUSTERS;
    this->EXPLORATION = EXPLORATION;    
    this->THRESHOLD = THRESHOLD;    
    this->K = K;    
}

void Cluster::update(const map<string, double>& readVec) {
    for (auto& [kmer, f] : readVec) {
        centroid[kmer] = (centroid[kmer] * size + f) / (size + 1);
    }
    
    // Atualizar norm
    norm = 0.0;
    for (auto& [_, v] : centroid) {
        norm += v * v;
    }
    norm = sqrt(norm);
    size++;
}

void Cluster::initialize(const map<string, double>& readVec) {
    centroid = readVec;
    size = 1;
    norm = 0.0;
    for (auto& [_, v] : centroid) {
        norm += v * v;
    }
    norm = sqrt(norm);
    initialized = true;
}
    
double Cluster::cosineSimilarity(const map<string, double>& readVec, const Cluster& c) {
    double dot = 0.0;
    double readNorm = 0.0;
    
    // Calcular produto escalar
    for (auto& [kmer, f] : readVec) {
        readNorm += f * f;
        auto it = c.centroid.find(kmer);
        if (it != c.centroid.end()) {
            dot += f * it->second;
        }
    }
    
    readNorm = sqrt(readNorm);
    if (readNorm == 0.0 || c.norm == 0.0) return 0.0;
    
    return dot / (readNorm * c.norm);
}

map<string, double> Cluster::normalizeKmerFreq(const map<string, int>& kmerCounts) {
    map<string, double> normalized;
    if (kmerCounts.empty()) return normalized;
    
    // Calcular total
    int total = 0;
    for (auto& [_, count] : kmerCounts) {
        total += count;
    }
    
    // Normalizar
    for (auto& [kmer, count] : kmerCounts) {
        normalized[kmer] = (double)count / total;
    }
    
    return normalized;
}

void Cluster::clusterize() {

    Reader reader;
    vector<Cluster> clusters;
    vector<ReadIndex> index;
    vector<int> seeds = {0, 200, 500, 300, 400};

    for (int i = 0; i < MAX_CLUSTERS; i++) 
        clusters.emplace_back(MAX_CLUSTERS, EXPLORATION, THRESHOLD, K);

    string fastaFile = "instance/mycoa.fasta";
    namespace fs = std::filesystem;
    fs::remove_all("clusters");
    fs::create_directory("clusters");

    index = reader.indexFasta(fastaFile);

    ifstream fasta(fastaFile, ios::binary);
    if (!fasta.is_open()) {
        cerr << "Erro ao abrir FASTA\n";
        return;
    }

    int totalReads = 0;
    int clusteredReads = 0;
    int seedReads = 0;
    map<string, int> kmerCounts;
    
    for (const auto& r : index) {

        string read = reader.loadReadAt(fasta, r.offset, r.length);

        if (read.empty()) continue;
        if (read.find('N') != string::npos) continue;

        totalReads++;
        kmerCounts.clear();

        int readLength = read.size();
        for (int i = 0; i < EXPLORATION; i++) {
            if (i + K > readLength) break;

            kmerCounts[read.substr(i, K)]++;

            int pos = readLength - K - i;
            if (pos < 0) break;
            kmerCounts[read.substr(pos, K)]++;
        }

        if (kmerCounts.empty()) continue;

        auto normalizedVec = normalizeKmerFreq(kmerCounts);

        bool seeded = false;
        for (int c = 0; c < MAX_CLUSTERS; c++) {
            if (r.id == seeds[c] && !clusters[c].initialized) {

                clusters[c].initialize(normalizedVec);
                seeded = true;
                seedReads++;

                cout << "Seed " << c + 1
                     << " inicializado com read " << r.id << endl;

                ofstream out("clusters/" + to_string(c + 1) + ".txt", ios::app);
                out << r.id << " " << r.offset << " " << r.length << "\n";
                out.close();

                break;
            }
        }

        if (seeded) continue;

        // competição pelo melhor cluster
        int best = -1;
        double bestSim = -1.0;
        for (int c = 0; c < MAX_CLUSTERS; c++) {
            if (!clusters[c].initialized) continue;

            double sim = cosineSimilarity(normalizedVec, clusters[c]);
            if (sim > bestSim) {
                bestSim = sim;
                best = c;
            }
        }

        if (best != -1 && bestSim >= THRESHOLD) {
            clusters[best].update(normalizedVec);
            clusteredReads++;

            ofstream out("clusters/" + to_string(best + 1) + ".txt", ios::app);
            out << r.id << " " << r.offset << " " << r.length << "\n";
            out.close();

        } else {
            // vai pro arquivo 0
            ofstream out("clusters/0.txt", ios::app);
            out << r.id << " " << r.offset << " " << r.length << "\n";
            out.close();
        }
    }

    fasta.close();
    
    cout << "\n=== ESTATÍSTICAS ===\n";
    cout << "Total de reads: " << totalReads << endl;
    cout << "Seeds: " << seedReads << endl;
    cout << "Clusterizadas: " << clusteredReads << endl;
    cout << "Não clusterizadas: "
         << totalReads - clusteredReads - seedReads << endl;

    for (int i = 0; i < MAX_CLUSTERS; i++) {
        cout << "Cluster " << i << ": "
             << clusters[i].size << " reads\n";
    }
}