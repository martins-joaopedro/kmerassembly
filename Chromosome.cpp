#include "header/Chromosome.hpp"
using namespace std;
Chromosome::Chromosome(int id, int genesNumber, int generation, int minOverlap, string status)
{
    this->id;
    this->genesNumber = genesNumber;
    this->generation = generation;
    this->minOverlap = minOverlap;
    this->status = status;
    this->overlapped = false;
    this->fitness = 0;

    // define ordem inicial;
    order.resize(genesNumber);
    for (int i = 0; i < genesNumber; i++)
        order[i] = i;
}

int Chromosome::computeFitness(vector<int> &ord)
{
    // so avalia uma vez, depois te ter preenchido os genes
    if(!overlapped)
        computeOverlapping();

    int sobreposition = 0;
    for (int i = 0; i + 1 < ord.size(); i++)
        sobreposition += overlaps[ord[i]][ord[i + 1]];
    
    return sobreposition;
}

void Chromosome::updateFitness(const vector<int>& ord) {
    int sobreposition = 0;
    for (int i = 0; i + 1 < ord.size(); i++)
        sobreposition += overlaps[ord[i]][ord[i + 1]];
    this->fitness = sobreposition;
}

// TODO: verificar se com ponteiro é melhor e se uma outra estrutura é melhor
// TODO: computar overlapping por kmer
void Chromosome::computeOverlapping()
{
    int MIN_OVERLAP = 5;
    int n = genes.size();
    overlaps.assign(n, vector<int>(n, 0));

    for (int A = 0; A < n; A++)
    {
        for (int B = 0; B < n; B++)
        {
            if (A == B)
                continue;

            int maxOverlap = 0;
            int minLen = min(genes[A].size(), genes[B].size());
            if (minLen < MIN_OVERLAP) continue;

            // testa se esta contido
            if (genes[A].find(genes[B]) != string::npos) {
                overlaps[A][B] = genes[B].size();
            } else {
                // procura o maior overlap entre sufixo de A e prefixo de B
                for (int len = MIN_OVERLAP; len <= minLen; len++)
                {
                    string prefixB = genes[B].substr(0, len);
                    string suffixA = genes[A].substr(genes[A].size() - len);
    
                    if (suffixA == prefixB)
                        maxOverlap = len;
                }
                overlaps[A][B] = maxOverlap;
            }
        }
    }
    overlapped = true;
}

void Chromosome::insertionMove(vector<int> &order, int i, int j)
{
    if (i == j)
        return;

    int gene = order[i];
    order.erase(order.begin() + i);

    if (j > i)
        j--;
    order.insert(order.begin() + j, gene);
}

// TODO: revisar
int Chromosome::deltaInsertion(int i, int j)
{
    if (i == j)
        return 0;

    vector<int> newOrder = order;
    insertionMove(newOrder, i, j);

    return computeFitness(newOrder) - fitness;
}

void Chromosome::mutate()
{

    bool improved = true;
    int iterations = 0;
    const int MAX_ITERATIONS = 1000;

    while (improved && iterations < MAX_ITERATIONS)
    {
        improved = false;
        iterations++;

        for (int i = 0; i < genesNumber; i++)
        {
            for (int j = 0; j < genesNumber; j++)
            {
                if (i == j)
                    continue;

                int delta = deltaInsertion(i, j);
                // se proporcionar melhoria realmente aplica
                if (delta > 0)
                {
                    insertionMove(order, i, j);
                    fitness += delta;
                    improved = true;

                    // DEBUG
                    cout << "Iteração " << iterations << ": Melhoria encontrada! Delta = " << delta << ", Fitness atual = " << fitness << endl;
                    break;
                }
            }
            if (improved)
                break;
        }
    }

    if (iterations >= MAX_ITERATIONS)
        cout << "Atingiu o limite máximo de iterações (" << MAX_ITERATIONS << ")" << endl;
    
    updateFitness(order);
}

// retorna o contig resolvido com as overlays aplicadas
vector<string> Chromosome::getFormedContigs()
{
    vector<string> contigs;

    if (order.empty())
        return contigs;

    string current = genes[order[0]];

    for (int i = 1; i < genesNumber; i++)
    {
        int a = order[i - 1];
        int b = order[i];

        int ov = overlaps[a][b];

        if (ov >= minOverlap)
        {
            // continua o contig
            current += genes[b].substr(ov);
        }
        else
        {
            // quebra de contig
            contigs.push_back(current);
            current = genes[b];
        }
    }

    contigs.push_back(current);

   /*  //TODO: verificar
    // atualizo meus genes pros contigs formados
    vector<string> newGenes;
    for(auto contig : contigs)
        newGenes.push_back(contig);
    genes = newGenes;

    genesNumber = genes.size();
    order.resize(genesNumber);
    iota(order.begin(), order.end(), 0);
    overlapped = false; */

    return contigs;
}

void Chromosome::printOverlapMatrix()
{
    cout << "\nMatriz de Overlaps (A->B):\n";
    cout << "   ";
    for (int i = 0; i < genesNumber; i++)
        cout << i << " ";
    cout << "\n";

    for (int i = 0; i < genesNumber; i++)
    {
        cout << i << ": ";
        for (int j = 0; j < genesNumber; j++)
        {
            cout << overlaps[i][j] << " ";
        }
        cout << endl;
    }
}