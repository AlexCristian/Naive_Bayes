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
#include <unordered_set>
#include <functional>
#include <algorithm>
#include "Genome.hpp"
#include "Class.hpp"
#include "Diskutil.hpp"
using namespace std;

template <class T>
class Class;

class NB {
private:
public:
    /**
     * Instantiates a NB instance.
     * @param _kmer_size The kmer size to use.
     * @param _save_dir  The root directory containing savefiles for all
     *                   classes.
     * @param _nthreads  The number of threads with which to simultaneously
     *                   process classes.
     */
    NB(int _kmer_size, path _save_dir, int _nthreads = 1);

    ~NB();

    /**
     * Populates NB's list of classes by crawling through the root directory and
     * adding all savefiles available - without loading the class data itself.
     */
    void load();

    /**
     * Directs all classes to save their data to disk.
     * Classes that haven't been loaded will be skipped.
     */
    void save();

    /**
     * Getter for the path of the save directory.
     * @return A boost::path object pointing to the save directory.
     */
    path getSavedir();

    /**
     * Getter for the kmer size to be used within NB.
     * @return An integer - the current kmer size.
     */
    int getKmerSize();

    /**
     * Adds a class to the classifier.
     * @param cl The class to be added.
     */
    void addClass(Class<kmer_t>* cl);

    /**
     * @return The number of classes known to the classifier.
     */
    size_t numClasses();

    /**
     * Looks up a class in the classifier by its ID.
     * @param cl_id Unique class ID.
     */
    Class<kmer_t>* getClass(string cl_id);

    /**
     * Adds a class to the update queue, instead of processing updates right
     * away. A class already waiting in the queue is not added again.
     * @param cl Class to add to the queue.
     * @return   True if the class was added, false if it was already queued.
     */
    bool addClassToUpdateQueue(Class<kmer_t>* cl);

    /**
     * @param cl Class to look up.
     * @return   True if the class is waiting in the update queue.
     */
    bool isClassQueuedForUpdate(Class<kmer_t>* cl);

    /**
     * Queues every class for an update (with no new genomes, an update just
     * reloads and re-saves the class).
     */
    void queueAllClassesForUpdate();

    /**
     * Triggers all queued classes to update.
     */
    void processClassUpdates();

    /**
     * Given an array of reads, classifies them into available classes and
     * computes membership likelihood for each.
     * @param reads Array of reads to be classified.
     */
    void classify(vector<Genome*> reads);

    /**
     * EXPERIMENTAL
     * Calculates the Levenshtein edit distance between two strings.
     * @param  src Source string.
     * @param  tgt Target string to change source into.
     * @return     Minimum number of edits to change the source into the target.
     */
    static int levenshteinDistance(string src, string tgt);

    static enum Debug{
      NO_LOG,
      LOG_SOME,
      LOG_ALL
    } debug_flag;

protected:
    unordered_map<string, Class<kmer_t>* > classes;
    int kmer_size, progress, total, nthreads;
    path save_dir;
    queue<Class<kmer_t>* > classesToProcess;
    unordered_set<Class<kmer_t>* > classesQueuedForUpdate;
    mutex classQueueAccess;

    static int min(int a, int b, int c);
    void classifyThreadController(vector<Genome*>& reads,
                                  const vector<kmer_t>& batch_kmers);
    void forEachRead(vector<Genome*>& reads, function<void(Genome*)> fn);
    void trainThreadController();
};

#endif /* NB_hpp */
