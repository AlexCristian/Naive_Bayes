//
//  NB.hpp
//  NB_C++
//
//  Created by Alexandru Cristian on 06/04/2017.
//
//

#ifndef NB_hpp
#define NB_hpp

#include <cstdio>
#include <fstream>
#include <iostream>
#include <queue>
#include <vector>
#include <thread>
#include <mutex>
#include "Genome.hpp"
#include "Class.hpp"
#include "Diskutil.hpp"
using namespace std;

template <class T>
class Class;

class NB {
private:
public:
    NB(int _kmer_size, path _save_dir, int _nthreads = 1);
    ~NB();

    void load();
    void save();

    path getSavedir();
    int getKmerSize();

    void addClass(Class<int>* cl);
    Class<int>* getClass(string cl_id);
    void addClassToUpdateQueue(Class<int>*);
    void processClassUpdates();

    void classify(vector<Genome*> reads);

    static int levenshteinDistance(string src, string tgt);

    static enum Debug{
      NO_LOG,
      LOG_SOME,
      LOG_ALL
    } debug_flag;

protected:
    unordered_map<string, Class<int>* > classes;
    int kmer_size, progress, total, nthreads;
    path save_dir;
    queue<Class<int>* > classesToProcess;
    mutex classQueueAccess;

    static int min(int a, int b, int c);
    void classifyThreadController(vector<Genome*>& reads);
    void trainThreadController();
};

#endif /* NB_hpp */
