#include <bits/stdc++.h>
using namespace std; // minimal

int K = 6;
int EXPLORATION = 150;
const int MAX_CLUSTERS = 5;
const double THRESHOLD = 0.2;  // Aumentado para clusterizar mais reads

struct Cluster {
    bool initialized = false;
    map<string, double> centroid;
    int size = 0;
    double norm = 0.0;

    void update(const map<string, double>& readVec) {
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
    
    void initialize(const map<string, double>& readVec) {
        centroid = readVec;
        size = 1;
        norm = 0.0;
        for (auto& [_, v] : centroid) {
            norm += v * v;
        }
        norm = sqrt(norm);
        initialized = true;
    }
};

double cosineSimilarity(const map<string, double>& readVec, const Cluster& c) {
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

map<string, double> normalizeKmerFreq(const map<string, int>& kmerCounts) {
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

void clusterize() {
    vector<Cluster> clusters(MAX_CLUSTERS);
    vector<int> seeds = {0, 200, 500, 300, 400};
    
    bool RECREATE = true;
    #include <filesystem>
    namespace fs = std::filesystem;

    if(RECREATE) {
        fs::remove_all("clusters");
        fs::create_directory("clusters");
    }

    ifstream file("instance/mycoa.fasta");
    if (!file.is_open()) {
        cerr << "Erro ao abrir arquivo\n";
    }
    
    string line, read;
    int id = 0;
    map<string, int> kmerCounts;
    
    // Contador para debug
    int totalReads = 0;
    int clusteredReads = 0;
    int seedReads = 0;

    while (getline(file, line)) {
        if (line[0] == '>' and !read.empty()) {
            // se não tem N processa read
            if (read.find('N') == string::npos) {
                totalReads++;
                kmerCounts.clear();
                int len = read.size();
                
                // Extrair k-mers das extremidades
                for (int i = 0; i < EXPLORATION; i++) {
                    if (i + K > len) break;
                    kmerCounts[read.substr(i, K)]++;

                    int pos = len - K - i;
                    if (pos < 0) break;
                    kmerCounts[read.substr(pos, K)]++;
                }

                if (kmerCounts.empty()) {
                    read.clear();
                    continue;
                }
                
                // Normalizar frequências
                auto normalizedVec = normalizeKmerFreq(kmerCounts);
                
                // Inicialização dos clusters a partir dos seeds
                bool seeded = false;
                for (int c = 0; c < MAX_CLUSTERS; c++) {
                    if (id == seeds[c] && !clusters[c].initialized) {
                        clusters[c].initialize(normalizedVec);
                        seeded = true;
                        seedReads++;
                        cout << "Seed " << c+1 << " inicializado com read " << id << endl;
                        break;
                    }
                }
                
                // Competição normal
                if (!seeded) {
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
                        cout << "Read " << id << " -> Cluster " << best << " (similaridade: " << bestSim << ")" << endl;
                    
                        string name = "clusters/" + to_string(best+1) + ".txt";
                        ofstream outCluster(name, ios::app);

                        if (outCluster.is_open()) {
                            outCluster << ">" << endl;
                            outCluster << read << endl;
                            outCluster.close();
                        }
                    } else if (best != -1) {
                        cout << "Read " << id << " rejeitado (similaridade: "  << bestSim << " < " << THRESHOLD << ")" << endl;
                        string name = "clusters/0.txt";
                        ofstream outCluster(name, ios::app);

                        if (outCluster.is_open()) {
                            outCluster << ">" << endl;
                            outCluster << read << endl;
                            outCluster.close();
                        }
                    }
                }
                // aumenta toda vez que eu processo uma read
                id++;
            }
            // limpa toda vez ao processar
            read.clear();
        } else read += line;
    }
    
    file.close();
    
    // Estatísticas
    cout << "\n=== ESTATÍSTICAS ===" << endl;
    cout << "Total de reads: " << totalReads << endl;
    cout << "Reads usados como seed: " << seedReads << endl;
    cout << "Reads clusterizados: " << clusteredReads << endl;
    cout << "Reads não clusterizados: " << totalReads - clusteredReads - seedReads << endl;
    
    for (int i = 0; i < MAX_CLUSTERS; i++) {
        cout << "Cluster " << i << ": " << clusters[i].size << " reads" << endl;
    }
    
    // Salvar resultados
    ofstream out("instance/clusters.txt");
    if (!out.is_open()) {
        cerr << "Erro ao criar arquivo de saída\n";
    }
    
    for (int i = 0; i < MAX_CLUSTERS; i++) {
        out << "Cluster " << i << "\n";
        out << "Tamanho: " << clusters[i].size << "\n";
        out << "Top 10 k-mers no centroide:\n";
        
        // Ordenar k-mers por importância
        vector<pair<string, double>> sortedCentroid(
            clusters[i].centroid.begin(), 
            clusters[i].centroid.end()
        );
        sort(sortedCentroid.begin(), sortedCentroid.end(),
             [](const auto& a, const auto& b) { return a.second > b.second; });
        
        for (int j = 0; j < min(10, (int)sortedCentroid.size()); j++) {
            out << sortedCentroid[j].first << ":" 
                << fixed << setprecision(4) << sortedCentroid[j].second << " ";
        }
        out << "\n---------------------------\n";
    }
    
    out.close();
    cout << "\nFinalizado. Clusters salvos em 'instance/clusters.txt'\n";
}

struct Chromossome {
    // vetor de reads
    vector<string> genes;

    // avalia as sobreposições
    int getFitness() {
        for(string read : genes) {

        }
    }
}

int main() {
    
    //clusterize();



    return 0;
}