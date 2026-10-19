//
//  Genome.hpp
//  NB_C++
//
//  Created by Alexandru Cristian on 01/05/2017.
//
//

#ifndef Genome_hpp
#define Genome_hpp

#include <atomic>
#include <cstdint>
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

// 2 bits per base, so k-mers up to k = 32 fit without collisions.
typedef uint64_t kmer_t;

class Genome {
private:
public:
  typedef pair<double, Class<kmer_t>* > score;

  // Orders scores by value, breaking ties by class ID so results don't depend
  // on hash map order or thread timing.
  struct ScoreLess {
    bool operator()(const score& a, const score& b) const;
  };
  typedef priority_queue<score, vector<score>, ScoreLess> pqueue;

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
   * Maps this read's k-mers to positions in a batch-wide sorted k-mer list,
   * keeping the order in which scoring visits them, and frees the k-mer map.
   * Must be called (once, single-threaded) before scoreClass().
   * @param batch_kmers Sorted, unique k-mers of every read in the batch.
   */
  void indexKmers(const vector<kmer_t>& batch_kmers);

  /**
   * Scores this read against a class, given the class's log counts for the
   * batch k-mer list passed to indexKmers(). Thread-safe.
   * @param cls      The class being scored.
   * @param batch_lg batch_lg[i] is the class's log count of batch k-mer i.
   */
  void scoreClass(Class<kmer_t>* cls, const vector<double>& batch_lg);

  /**
   * Reads this genome's k-mer counts into a vector sorted by k-mer, without
   * keeping them in the Genome.
   * @param out Receives (k-mer, count) pairs, sorted and unique.
   */
  void loadSortedKmerCounts(vector<pair<kmer_t, int> >& out);

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
   * Returns the path on disk of the kmr count file registered for this read.
   * @return The file's path.
   **/
  path getKmrPath();

  /**
   * @return     The length of the FASTA sequence associated to this genome.
   */
  int size();

  /**
   * @return     A hashmap containing the counts (values) associated to each
   *             kmer (key).
   */
  unordered_map<kmer_t, int>& getKmerCounts();

  /**
   * When the STORE_ALL_NUMERATORS flag is false, returns the highest ranked class.
   */
  score getMaximum();

  static bool STORE_ALL_NUMERATORS;

  /**
   * Expected k-mer length in count files. When nonzero, k-mers of any other
   * length abort the program instead of being silently mis-encoded.
   */
  static unsigned int KMER_SIZE;

  /**
   * Whether scores include the class prior, log(1 + number of training
   * genomes in the class). On by default.
   */
  static bool USE_PRIOR;
protected:
  atomic<bool> kmersLoaded{false};
  bool sequenceLoaded = false;
  pqueue numerator;
  mutex numeratorAccess, loadKmersLock, loadSequenceLock;
  path sequence_path, kmr_path;
  string *sequence;
  unordered_map<kmer_t, int> *kmer_counts;

  // Set by indexKmers(): each k-mer's position in the batch list and its
  // count, in the order the k-mer map iterates (so sums match exactly).
  vector<uint32_t> batch_index;
  vector<int> batch_count;

  // Parses a jellyfish "dump -c" file, calling add(kmer, count) for each
  // valid k-mer.
  template <class F>
  void parseKmrFile(F add);
  long long int max(long long int a,
                    long long int b,
                    long long int c,
                    long long int d);

  double maximumNumerator;
  Class<kmer_t>* maximumNumeratorClass = NULL;
};

#include "Class.hpp"
#include "NB.hpp"

#endif /* Genome_hpp */
