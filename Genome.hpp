//
//  Genome.hpp
//  NB_C++
//
//  Created by Alexandru Cristian on 01/05/2017.
//
//

#ifndef Genome_hpp
#define Genome_hpp

#include <vector>
#include <utility>
#include <unordered_map>
#include <queue>
#include <mutex>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/sum_kahan.hpp>
#include <boost/accumulators/statistics.hpp>
#include "Diskutil.hpp"
using namespace std;
using namespace boost::accumulators;


template <class T>
class Class;

class Genome {
private:
public:
  typedef pair<double, Class<int>* > score;
  typedef priority_queue<score, vector<score>, less<score> > pqueue;

  /**
   * Instantiates a Genome instance.
   * @param  kmr_path Path to file containing serialized k-mers.
   * @param  sequence_path Path to file containing DNA sequence.
   */
  Genome(path kmr_path, path sequence_path);

  ~Genome();

  /**
   * Loads kmer counts from the file specified in the "kmr_path" instantiation
   * parameter. When triggered by multiple threads (on the same object), will
   * wait for the initial thread to finish loading data and then simply return
   * to all threads waiting.
   */
  void loadKmerCounts();

  /**
   * Reads FASTA sequence from disk using the file specified in the
   * "sequence_path" instantiation parameter.
   * NOT THREAD SAFE, at the moment.
   */
  void loadSequence();

  /**
   * Deletes the sequence and/or the kmer counts from memory.
   * Leaves the numerator queue intact.
   */
  void unload();

  /**
   * Computes optimal sequence alignment between the current instance and the
   * Genome supplied in the parameter. Using the Smith-Waterman algorithm.
   * @param  seq Genome instance to compare.
   * @return     The pair's alignment score (greater value is better).
   */
  long long int computeAlignmentScore(Genome* seq);

  /**
   * Takes in a Class instance, and returns the prior numerator; a score
   * correlated to the probability of the genome instance pertaining to the
   * supplied class.
   * @param  class Instance of Class object to use for the computation.
   * @return       Score correlated to probability of class membership.
   */
  void computeClassificationNumerator(Class<int>* cls);

  /**
   * Gets the character located at position "pos" in the DNA sequence.
   * @param  pos Integer, representing the position we want to query.
   * @return     The base associated to that position.
   */
  string charAt(int pos);

  /**
   * Getter for the FASTA sequence path.
   * @return     The path to the FASTA file.
   */
  string& getSequence();

  /**
   * Obtains a priority queue populated with class-probability pairs.
   * Classes most likely to hold the genome will be at the top of the queue.
   * @return     For each class, the probability of this genome pertaining
   *             to it.
   */
  pqueue getConfidences();

  /**
   * @return     The length of the FASTA sequence associated to this genome.
   */
  int size();

  /**
   * @return     A hashmap containing the counts (values) associated to each
   *             kmer (key).
   */
  unordered_map<int, int>& getKmerCounts();
protected:
  bool kmersLoaded = false, sequenceLoaded = false;
  pqueue numerator;
  mutex numeratorAccess, loadKmersLock, loadSequenceLock;
  path sequence_path, kmr_path;
  string *sequence;
  unordered_map<int, int> *kmer_counts;
  long long int max(long long int a,
                    long long int b,
                    long long int c,
                    long long int d);
};

#include "Class.hpp"
#include "NB.hpp"

#endif /* Genome_hpp */
