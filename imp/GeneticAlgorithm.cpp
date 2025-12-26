#include "../header/GeneticAlgorithm.hpp"
using namespace std;

GeneticAlgoritm::GeneticAlgoritm(int chromosomesNumber, int genesNumber, int islandsNumber, int MIN_OVERLAPING) {
    this->chromosomesNumber = chromosomesNumber;
    this->genesNumber = genesNumber;
    this->islandsNumber = islandsNumber;
    this->MIN_OVERLAPING = MIN_OVERLAPING;
} 

vector<Chromosome> GeneticAlgoritm::crossover(Chromosome c1, Chromosome c2) {

    // pega um numero entre 25% - 75% do cromossomo 
    int c1Genes = c1.genes.size();
    int c2Genes = c2.genes.size();
    int minGenesNumber = min(c1Genes, c2Genes);

    int random = rand() % (int)(0.5 * minGenesNumber);
    int point = (int)(0.25 * minGenesNumber) + random;

    int generation = 0;

    vector<Chromosome> offspring;
    offspring.emplace_back(c1.id, genesNumber, generation, MIN_OVERLAPING, "");
    offspring.emplace_back(c2.id, genesNumber, generation, MIN_OVERLAPING, "");
    
    // sempre padronizo para o tamanho inicial mesmo que ja tenha formado contig
    //c1.genes.resize(genesNumber);
    //c2.genes.resize(genesNumber);
    
    offspring[0].genes.resize(genesNumber);
    offspring[1].genes.resize(genesNumber);
    
    // TODO: a ordem aq importa? acabei de criar eles e to misturando os genes
    //offspring[0].order.resize(genesNumber);
    //offspring[1].order.resize(genesNumber);
    
    // copia a primeira parte de cada cromossomo 
    for (int i = 0; i < point; i++) {
        offspring[0].genes[i] = c1.genes[i];
        offspring[1].genes[i] = c2.genes[i];
        //offspring[0].order[i] = c1.order[i];
        //offspring[1].order[i] = c2.order[i];
    }

    // copia a segunda parte
    for (int i = point; i < genesNumber; i++) {
        offspring[1].genes[i] = c1.genes[i];
        offspring[0].genes[i] = c2.genes[i];
        //offspring[1].order[i] = c1.order[i];
        //offspring[0].order[i] = c2.order[i];
    }
        
    return offspring;
}

string loadReadAt(ifstream& fasta, uint64_t offset, uint32_t len) {
    // Limpa flags e posiciona no início da sequência
    fasta.clear();
    fasta.seekg(offset, ios::beg);
    
    string seq;
    string line;
    
    // Lê até completar o tamanho ou encontrar próxima header
    while (seq.size() < len && getline(fasta, line)) {
        // Se encontrar um novo header, para
        if (!line.empty() && line[0] == '>') {
            break;
        }
        
        // Remove Windows CRLF se existir
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }
        
        // Adiciona à sequência
        seq += line;
        
        // Se já leu o suficiente, para
        if (seq.size() >= len) {
            break;
        }
    }
    
    // Garante o tamanho exato (se leu mais)
    if (seq.size() > len) {
        seq.resize(len);
    }
    
    return seq;
}

vector<Chromosome> GeneticAlgoritm::createPopulation(int islandNumber, int era) {

    int clusterId = (islandNumber % 1) + 1;
    string filename = "clusters/" + to_string(clusterId) + ".txt";
    string fastaFile = "instance/mycoa.fasta";

    cout << "Lendo população do cluster " << clusterId 
         << " (ilha " << islandNumber << ", era " << era << ")\n";

    ifstream cluster(filename);
    ifstream fasta(fastaFile, ios::binary);

    if (!cluster.is_open()) {
        cerr << "Erro ao abrir cluster: " << filename << "\n";
        return {};
    }
    if (!fasta.is_open()) {
        cerr << "Erro ao abrir FASTA: " << fastaFile << "\n";
        return {};
    }

    vector<Chromosome> population;

    int id;
    uint64_t offset;
    uint32_t len;

    int processedReads = 0;
 
    Chromosome currentChromo(0, genesNumber, era, MIN_OVERLAPING, "");

    while (cluster >> id >> offset >> len) {

        // marco como read usada
        if (id >= usedReads.size())
            usedReads.resize(id + 1, 0);

        if (usedReads[id])
            continue;
        usedReads[id] = 1;

        string read = loadReadAt(fasta, offset, len);

        if (read.size() != len) {
            cerr << "ERRO: Read " << id
                 << " tamanho incorreto (" << read.size()
                 << " != " << len << "), descartando\n";
            continue;
        }

        currentChromo.genes.push_back(read);
        
        processedReads++;

        if ((int)currentChromo.genes.size() == genesNumber) {

            currentChromo.id = population.size();
            population.push_back(currentChromo);

            cout << "Cromossomo " << currentChromo.id
                 << " criado com " << genesNumber << " genes\n";

            if ((int)population.size() >= chromosomesNumber)
                break;

            currentChromo = Chromosome(
                population.size(), genesNumber, era, MIN_OVERLAPING, ""
            );
        }
    }

    if (!currentChromo.genes.empty()) {
        cerr << "AVISO: Cromossomo incompleto descartado ("
             << currentChromo.genes.size() << "/" << genesNumber << " genes)\n";
    }

    cout << "População criada: " << population.size()
         << " cromossomos, " << processedReads << " reads processadas\n";

    if (population.empty()) {
        cerr << "ERRO: população vazia\n";
        return {};
    }

    return population;
}


// TODO: cuidar para nao ultrapassar limites do arquivo 
void GeneticAlgoritm::endMiscenegationPhase(vector<Chromosome> &pop, int islandNumber, int era) {
    
    string filename = "clusters/" + to_string(1) + ".txt";
    string fastaFile = "instance/mycoa.fasta";
    ifstream cluster(filename);
    ifstream fasta(fastaFile, ios::binary);

    int id;
    uint64_t offset;
    uint32_t len;

    for(auto &cromo : pop) {
        cromo.groupContigs();

        cout << cromo.genes.size() << endl;
        
        // completo numero de genes
        while (cluster >> id >> offset >> len) {

            if(cromo.genes.size() >= genesNumber)
                break;

            // marco como read usada
            if (id >= usedReads.size())
                usedReads.resize(id + 1, 0);

            if (usedReads[id])
                continue;
            //usedReads[id] = 1;

            string read = loadReadAt(fasta, offset, len);

            if (read.size() != len) {
                cerr << "ERRO: Read " << id
                    << " tamanho incorreto (" << read.size()
                    << " != " << len << "), descartando\n";
                continue;
            }

            // adiciono mais genes
            cromo.genes.push_back(read);
        }
    }
}


void GeneticAlgoritm::start() {

    cout << "INICIANDO ALGORITMO GENETICO" << endl;

    ofstream outFile("out.txt");

    // iterações do GA procurando melhorias
    int GA_IT = 10;

    // marca um período para união dos contigs e inserção de novos genes;
    int MISCEGENATION = 5;

    // marca o período para a clusterização 
    int MAX_ERAS = 3;

    float ELITE_PERC = 0.3;
    bool DEBUG = true;

    for(int era=0; era<MAX_ERAS; era++) {

        // processo o GA-memetico pra todas as ilhas 
        for(int clusterId = 0; clusterId < 1; clusterId++) {
            // inicial população para aquela ilha
            vector<Chromosome> pop = createPopulation(clusterId, era);
    
            // evolui a população inicial e ao final 
            for(int it=0; it<GA_IT; it++) {
                
                // marca periodo para união em contigs e inserção de novos genes 
                for(int it=0; it<MISCEGENATION; it++) {

                    // ordeno os cromossomos
                    sort(pop.begin(), pop.end(), [](const Chromosome& a, const Chromosome& b) {
                        return a.fitness > b.fitness;
                    });
        
                    outFile << "POPULAÇÂO ORDENADA" << endl;
                    for (int i = 0; i < pop.size(); i++) {
                        outFile << pop[i].fitness << endl;
                    }
        
                    // CROSSOVER: um da elite e um geral
                    if (pop.size() >= 2) {
                        int idx1 = rand() % max(1, (int)(pop.size() * ELITE_PERC));
                        int idx2 = rand() % pop.size();
                        
                        while (idx2 == idx1) 
                            idx2 = rand() % pop.size();
                        
                        vector<Chromosome> offspring = crossover(pop[idx1], pop[idx2]);
                        
                        offspring[0].status = "GENERATION";
                        offspring[1].status = "GENERATION";
        
                        // substitui 2 piores candidatos da população -> dois ultimos indices
                        for(int i=0; i<2; i++) {
                            int id = chromosomesNumber -i -1;
                            pop[id] = offspring[i];
                        }
                    
                        outFile << "\n=== CROSSOVER ===" << endl;
                        outFile << "Pai 1 (índice " << idx1 << ") x Pai 2 (índice " << idx2 << ")" << endl;
                        outFile << "Ponto de corte: " << (int)(0.25 * genesNumber) << " a " << (int)(0.75 * genesNumber) << endl;
                        outFile << "Descendente 1 fitness: " << offspring[0].computeFitness(offspring[0].order) << endl;
                        outFile << "Descendente 2 fitness: " << offspring[1].computeFitness(offspring[1].order) << endl;
                    }
        
                    // MUTATION: todos cromossomos sofrem mutação individual -> busca local
                    #pragma omp parallel for
                    for (int i = 0; i < pop.size(); i++) {
                        outFile << "\n=== MUTANDO CROMOSSOMO " << i << " ===\n";
                        pop[i].mutate();
                    
                        if(DEBUG) {
                            outFile << ">>>>>> INICIAL " << endl;
                            for(auto gene : pop[i].genes) outFile << gene << endl;     
                            
                            //outFile << "\nOrdem final: ";
                            //for (int x : pop[i].order) outFile << x << " ";
        
                            outFile << "\nFitness final: " << pop[i].fitness << endl;
                            
                            vector<string> contigs = pop[i].getFormedContigs();
                            outFile << ">>>>>> " << contigs.size() << "contigs formados" << endl;
                            for(auto contig : contigs) outFile << contig << endl;
                        }
                    }
                }
            
                // final da etapa de miscigenação
                // formo contigs e completo todos os cromossomo com novos genes 
                outFile << "MISCIGENAÇÃO" << endl;

                endMiscenegationPhase(pop, clusterId, era);

                for(auto cromo : pop) 
                    for(auto gene : cromo.genes) 
                        outFile << gene << endl;
            }

            outFile.close();
            
            // salva todos os resultados obtidos nas eras 
            string file = "eras/" + to_string(era) + ".txt"; 
            ofstream eraFile(file);
            
            int totalGenes = 0;
            for (int i = 0; i < pop.size(); i++) {
                for(auto gene : pop[i].genes) {
                    eraFile << ">" << endl << gene << endl;
                    totalGenes++;
                }
            }

            eraFile.seekp(0);
            eraFile << totalGenes << endl << ">" << endl;
            eraFile.close();
        }
    }
}