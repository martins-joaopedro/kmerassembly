#include "header/GeneticAlgorithm.hpp"
using namespace std;

GeneticAlgoritm::GeneticAlgoritm(int chromosomesNumber, int genesNumber, int islandsNumber) {
    this->chromosomesNumber = chromosomesNumber;
    this->genesNumber = genesNumber;
    this->islandsNumber = islandsNumber;
} 

vector<Chromosome> GeneticAlgoritm::crossover(Chromosome c1, Chromosome c2) {

    // pega um numero entre 25% - 75% do cromossomo 

    int c1Genes = c1.genes.size();
    int c2Genes = c2.genes.size();
    int minGenesNumber = min(c1Genes, c2Genes);

    int random = rand() % (int)(0.5 * minGenesNumber);
    int point = (int)(0.25 * minGenesNumber) + random;

    int generation = 0;
    int minOverlap = 5;
    vector<Chromosome> offspring;
    offspring.emplace_back(c1.id, genesNumber, generation, 5, "");
    offspring.emplace_back(c2.id, genesNumber, generation, 5, "");
    
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

// REFACTOR: nao ta fazendo sentido a remontagem com base na era, to montando os mesmos cromossomos, 
vector<Chromosome> GeneticAlgoritm::createPopulation(int islandNumber, int era) {

    string filename;

    // nao vou clusterizar a cada era
    /* if(era == 0) {
        filename = "eras/" + to_string(islandNumber) + ".txt";
    }   */
    filename = "clusters/" + to_string(islandNumber) + ".txt";

    ifstream file(filename);
    
    if (!file.is_open()) {
        cerr << "Erro ao abrir arquivo\n";
        return vector<Chromosome>();
    }

    // TODO: resolver esse mockup de numero do começo do cluster
    string line, read;
    int totalReads = 5000;
    int interval = chromosomesNumber * genesNumber;
    
   
    
    cout << totalReads << endl;

    vector<Chromosome> population;
    
    int id = 0;
    int generation = 0;
    int MINOVERLAP = 5;
    
    Chromosome chromo(0, genesNumber, generation, MINOVERLAP, "");

    // marca a quantidade de reads adicionadas
    int processedReads = 0;

    // marca o valor atual da read sendo lida
    int currentRead = 0;

    // se 100 / 3 = 33 interval
    int K = (totalReads / interval);
    K = 1;
    
    while (getline(file, line)) {
        //cout << processedReads << endl;
        // toda vez que encontrar uma read processa
        if(line[0] == '>') {
            if(!read.empty()) {
                if(currentRead % K == 0) {
                    cout << currentRead << endl;
                    chromo.genes.push_back(read);
                    processedReads++;
                }
                 
                // sempre atualizo
                currentRead++;
                read.clear();
            }
        } else read += line;

        // quando eu terminar de ler as quantidades eu add na pop
        if(processedReads % genesNumber == 0 && !chromo.genes.empty()) {
            population.push_back(chromo);
            chromo = Chromosome(id, genesNumber, generation, MINOVERLAP, "");
            id++;
        }

        if(processedReads > 0 and processedReads == chromosomesNumber*genesNumber) {
            cout << "população calculada" << endl;
            break;
        }
    }

    file.close();
    return population;
}

void GeneticAlgoritm::start() {

    cout << "INICIANDO ALGORITMO GENETICO" << endl;

    ofstream outFile("out.txt");

    int MAX_GA_IT = 5;
    int MAX_ERAS = 3;
    float ELITE_PERC = 0.3;
    bool DEBUG = true;

    for(int era=0; era<MAX_ERAS; era++) {
        for(int clusterId = 0; clusterId < 1; clusterId++) {

            // inicial população para aquela ilha
            vector<Chromosome> pop = createPopulation(clusterId, era);
    
            // evolui a população inicial   
            for(int it=0; it<MAX_GA_IT; it++) {
                
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
                    int idx1 = rand() % min(1, (int)(pop.size() * ELITE_PERC));
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
                int chromoId = 1;
                for(Chromosome& chromo : pop) {           
                    outFile << "\n=== MUTANDO CROMOSSOMO " << chromoId << " ===\n";
                    chromo.mutate();
                    chromoId++;
    
                    if(DEBUG) {
                        outFile << ">>>>>> INICIAL " << endl;
                        for(auto gene : chromo.genes) outFile << gene << endl;     
                        
                        //outFile << "\nOrdem final: ";
                        //for (int x : chromo.order) outFile << x << " ";
    
                        outFile << "\nFitness final: " << chromo.fitness << endl;
                        
                        vector<string> contigs = chromo.getFormedContigs();
                        outFile << ">>>>>> " << contigs.size() << "contigs formados" << endl;
                        for(auto contig : contigs) outFile << contig << endl;
                    }
                }
    
            }

            outFile.close();
            
            //salva a elite dos cromossomos de cada cluster nas eras
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