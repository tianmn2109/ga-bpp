/*
 * =====================================================================================
 *
 *       Filename:  ga.hh
 *
 *    Description:  genetic algorithm to solve bin packing problem
 *
 *        Version:  1.0
 *        Created:  05/01/13 22:08:07
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Wei Wei, 
 *   Organization:  SJTU
 *
 * =====================================================================================
 */

#include <iostream>  // std::cout
#include <fstream>
#include <vector>    // std::vector
#include <algorithm> // std::find
#include <utility>
using namespace std;

#define POP_SIZE 20
#define GENERATION_NUM 1000

bool comp(const pair<int, int> l, const pair<int, int> r) {
    return l.second < r.second;
}

struct bpp {
    int n;
    int nb;
    int c;
    int *w;
    int *cap;
    int **alloc;
};

int init(struct bpp &b, char *filename, int nb)
{
    b.nb = nb;
    ifstream fin(filename);
    fin >> b.n;
    fin >> b.c;
    b.w = new int[b.n];
    for (int i = 0; i < b.n; i++) {
        fin >> b.w[i];
    }
    fin.close();
    b.cap = new int[b.nb];
    b.alloc = new int*[b.nb];
    for (int i = 0; i < b.nb; i++) {
        b.cap[i] = b.c;
        b.alloc[i] = new int[b.n];
        for (int j = 0; j < b.n; j++) {
            b.alloc[i][j] = 0;
        }
    }
    return 0;
}

struct phenotype {
    vector<vector<int> > chrom;
    vector<int> cap;
    int num;
    phenotype() {
        num = 0;
        chrom.reserve(50);
        cap.reserve(50);
    }
    phenotype(int nb, int n, int c) {
        chrom.resize(nb);
        for (int i = 0; i < nb; i++) {
            chrom[i].resize(nb);
            fill(chrom[i].begin(), chrom[i].end(), 0);
        }
        cap.resize(n);
        fill(cap.begin(), cap.end(), c);
        num = 0;
    }
    ~phenotype() {
    }
};

struct phenotype population[POP_SIZE];
struct phenotype tmp[6]; 

void swap(int &a, int &b)
{
    int tmp = a;
    a = b;
    b = tmp;
}

// fisher-yates shuffle
void shuffle(int *array, int *index, int first, int last)
{
    int i, j;
    for (i = last - first; i >= 1; --i) {
        j = rand() % (i + 1) + first; // random integer with 0<=j<=i
        swap(array[i], array[j]);
        swap(index[i], index[j]);
    }
}

void initialize(struct phenotype *p, const struct bpp b)
{
    int i, j, k;
    // int *cap = new int[b.nb];
    int *tmp = new int[b.n];
    int *index = new int[b.n];
    for (i = 0; i < POP_SIZE; i++) {
        p[i].chrom.resize(b.nb);
        for (j = 0; j < b.nb; j++) {
            p[i].chrom[j].resize(b.n);
            for (k = 0; k < b.n; k++) {
                p[i].chrom[j][k] = 0;
            }
        }
    }
    for (i = 0; i < POP_SIZE; i++) {
        for (j = 0; j < b.n; j++) {
            index[j] = j;
        }
        // fill(cap, cap + b.nb, b.c);
        p[i].cap.resize(b.nb);
        fill(p[i].cap.begin(), p[i].cap.end(), b.c);
        memcpy(tmp, b.w, sizeof(int) * b.n);
        shuffle(tmp, index, 0, b.n - 1);
        for (j = 0; j < b.n; j++) {
            for (k = 0; k < b.nb; k++) {
                if (tmp[j] <= p[i].cap[k]) {
                    p[i].chrom[k][index[j]] = 1;
                    p[i].cap[k] -= tmp[j];
                    if (p[i].num < k)
                        p[i].num = k;
                    break;
                }
            }
        }
        p[i].num++;
    }
    delete[] tmp;
    delete[] index;
}

int select()
{
    return rand() % POP_SIZE;
}

int crossOver(const struct phenotype p1, const struct phenotype p2, struct phenotype &p, struct bpp b)
{
    int i, j, k;
    p.num = 0;
    // randomly pick num / 2 bins
    // Fisher-Yates shuffle
    vector<int> a;
    for (i = 0; i < p1.num; i++)
        a.push_back(i);
    for (i = p1.num - 1; i>= 1; i--) {
        j = rand() % (i + 1);
        swap(a[j], a[i]);
    }
    // fill
    fill(p.cap.begin(), p.cap.end(), b.c);
    for (i = 0; i < b.nb; i++) {
        fill(p.chrom[i].begin(), p.chrom[i].end(), 0);
    }
    for (i = 0; i < p1.num / 2; i++) {
        p.cap[i] = p1.cap[a[i]];
        p.chrom[i] = p1.chrom[a[i]];
    }
    // get the rest of objects
    int flag;
    vector<int> rest;
    for (i = 0; i < b.n; i++) {
        flag = 0;
        for (j = p1.num / 2; j < p1.num; j++) {
            if (p1.chrom[a[j]][i] == 1)
                flag = 1;
        }
        if (flag == 1) { 
            rest.push_back(i);
        }
    }

    /* 
    // print rest
    cout << endl << "REST:" << endl;
    sort(rest.begin(), rest.end());
    for (int i = 0; i < rest.size(); i++) {
        cout << rest[i] << " ";
    }
    cout << endl;
    */

    // put the rest objects into the bins
    // first approach
    // let's first try simple first fit
    // it can not have better performance than ffd in a single run!
    /* p.num = 0;
    // cout << rest.size() << endl;
    for (int i = 0; i < rest.size(); i++) {
        for (int j = 0; j < b.nb; j++) {
            if (b.w[rest[i]] <= p.cap[j]) {
                p.cap[j] -= b.w[rest[i]];
                p.chrom[j][rest[i]] = 1;
                if (p.num < j)
                    p.num = j;
                break;
            }
        }
    }
    p.num++;
    // cout << p.num << endl;
    */

    // second approach
    // first try to preserve p2's allocation to the maximum
    // get a list of p2's objs
    vector<vector<int> > objs;
    vector<int>::iterator it;
    for (i = 0; i < b.nb; i++) {
        vector<int> v;
        for (j = 0; j < b.n; j++) {
            if (p2.chrom[i][j] == 1) {
                if (!rest.empty()) {
                    it = find(rest.begin(), rest.end(), j);
                    if (*it == j) {
                        v.push_back(j);
                    }
                }
            }
        }
        if (v.size() != 0)
            objs.push_back(v);
    }

    /* 
    // print objs
    cout << endl << "SHOW VEC:" << endl;
    for (int i = 0; i < objs.size(); i++) {
        for (int j = 0; j < objs[i].size(); j++) {
            cout << objs[i][j] << " ";
        }
        cout << endl;
    }*/
    

    // sort the items in non-decreasing order 
    vector<pair<int, int> > pairs;
    for (i = 0; i < objs.size(); i++) {
        int weight = 0;
        for (int j = 0; j < objs[i].size(); j++) {
            weight += b.w[objs[i][j]];
        }
        pair<int, int> tmp (i, weight);
        pairs.push_back(tmp);
    }
    sort(pairs.begin(), pairs.end(), comp);

    // reconstruct items
    // vector<vector<int> > items;
    // for (int i = 0; i < pairs.size(); i++) {
    //    items.push_front(objs[pairs[i]]);
    //}

    //rest.clear();
    //pairs.clear();
    //objs.clear();

    // put objs into bins using first fit decreasing
    for (i = pairs.size() - 1; i >= 0; i--) {
        // cout << "Pair: " << i << " " << pairs[i].first << " " << pairs[i].second << endl;
        for (j = 0; j < b.nb; j++) {
            if (pairs[i].second <= p.cap[j]) {
                // cout << p.cap[j] << endl;
                p.cap[j] -= pairs[i].second;
                for (k = 0; k < objs[pairs[i].first].size(); k++) {
                    p.chrom[j][objs[pairs[i].first][k]] = 1;
                    if (p.num < j)
                        p.num = j;
                }
                break;
            }
        }
    }
    p.num++;
    
    return 0;
}

void statistics(const struct phenotype *p)
{
    for (int i = 0; i < POP_SIZE; i++)
        cout << p[i].num << " ";
    cout << endl;
}

int firstFit(struct bpp &b)
{
    int tmp = 0;
    for (int i = 0; i < b.n; i++) {
        for (int j = 0; j < b.nb; j++) {
            if (b.w[i] <= b.cap[j]) {
                b.alloc[j][i] = 1;
                b.cap[j] -= b.w[i];
                if (tmp < j) {
                    tmp = j;
                }
                break;
            }
        }
    }
    return tmp + 1;
}

void printWeight(const struct phenotype p, const struct bpp b)
{
    cout << "WEIGHT:" << endl;
    for (int i = 0; i < b.n; i++)
        cout << b.w[i] << " ";
    cout << endl;
}

void printCapacity(const struct phenotype p, const struct bpp b)
{
    cout << "CAPACITY:" << endl;
    for (int i = 0; i < b.nb; i++)
        cout << p.cap[i] << " ";
    cout << endl;
}

void printAlloc(const struct phenotype p, const struct bpp b)
{
    cout << "ALLOCATION:" << endl;
    for (int i = 0; i < b.nb; i++) {
        for (int j = 0; j < b.n; j++) {
            if (j % 5 == 0)
                cout << " ";
            cout << p.chrom[i][j];
        }
        cout << endl;
    }
}

// calculate capacity of each bin based on obj allocation
void calcCap(const struct phenotype p, const struct bpp b)
{
    int *cap = new int[b.nb];
    fill(cap, cap + b.nb, b.c);
    for (int i = 0; i < b.n; i++) {
        for (int j = 0; j < b.nb; j++) {
            if (p.chrom[j][i] == 1) {
                cap[j] -= b.w[i];
            }
        }
    }
    for (int i = 0; i < b.nb; i++) {
        cout << cap[i] << " ";
    }
    cout << endl;
    delete[] cap;
}

// test whether a solution is valid or not
int isValidSolution(const struct phenotype p, const struct bpp b)
{
    int flag = 1;
    // test if each obj is put in exactly one bin
    for (int i = 0; i < b.n; i++) {
        int tmp = 0;
        for (int j = 0; j < b.nb; j++) {
            if (p.chrom[j][i] == 1)
                tmp++;
        }
        if (tmp != 1) {
            cout << "Each obj should be put into exactly one bin!" << endl;
            flag = 0;
        }
    }
    // test if each obj is put without violation
    int *cap = new int[b.nb];
    fill(cap, cap + b.nb, b.c);
    for (int i = 0; i < b.n; i++) {
        for (int j = 0; j < b.nb; j++) {
            if (p.chrom[j][i] == 1) {
                cap[j] -= b.w[i];
            }
        }
    }
    for (int i = 0; i < b.nb; i++) {
        if (cap[i] < 0) {
            cout << "Each obj must be put without violation!" << endl;
            flag = 0;
        }
    }
    delete[] cap;
    return flag;
}

void deepcopy(struct phenotype &p1, const struct phenotype p2, const struct bpp b)
{
    /*
    // initialize p1 if not defined
    if (p1.chrom == NULL) {
        // cout << "NULL Pointer" << endl;
        p1.chrom = new int*[b.nb];
        for (int i = 0; i < b.nb; i++)  {
            p1.chrom[i] = new int[b.n];
        }
    }
    if (p1.cap == NULL) {
        // cout << "NULL Pointer" << endl;
        p1.cap = new int[b.n];
    }

    for (int i = 0; i < b.nb; i++) {
        memcpy(p1.chrom[i], p2.chrom[i], sizeof(int) * b.n);
    }
    memcpy(p1.cap, p2.cap, sizeof(int) * b.nb);
    */
    p1.num = p2.num;
    p1.chrom = p2.chrom;
    p1.cap = p2.cap;
}

void mutate(struct phenotype &p1, struct phenotype &p, struct bpp b)
{
    p.num = 0;
    // select num / 10 bins randomly and then do first fit decreasing
    vector<int> a;
    for (int i = 0; i < p1.num; i++)
        a.push_back(i);
    for (int i = p1.num - 1; i >= 1; i--) {
        int j = rand() % (i + 1);
        swap(a[j], a[i]);
    }
    fill(p.cap.begin(), p.cap.end(), b.c);
    for (int i = 0; i < b.nb; i++) {
        fill(p.chrom[i].begin(), p.chrom[i].end(), 0);
    }
    for (int i = 0; i < 0.3 * p1.num; i++) {
        p.cap[i] = p1.cap[a[i]];
        p.chrom[i] = p1.chrom[a[i]];
    }
    // get the rest
    vector<int> rest;
    for (int i = 0.3 * p1.num; i < p1.num; i++) {
        for (int j = 0; j < b.n; j++) {
            if (p1.chrom[a[i]][j] == 1) {
                rest.push_back(j);
            }
        }
    }
    // sort rest first
    vector<pair<int, int> > pairs;
    for (int i = 0; i < rest.size(); i++) {
        pair<int, int> tmp (rest[i], b.w[rest[i]]);
        pairs.push_back(tmp);
    }
    sort(pairs.begin(), pairs.end(), comp);
    // do first fit decreasing to rest
    for (int i = pairs.size() - 1; i >= 0; i--) {
        for (int j = 0; j < b.nb; j++) {
            if (pairs[i].second <= p.cap[j]) {
                p.cap[j] -= pairs[i].second;
                p.chrom[j][pairs[i].first] = 1;
                if (p.num < j)
                    p.num = j;
                break;
            }
        }
    }
    p.num++;
}
