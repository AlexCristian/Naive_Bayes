//
//  Class.hpp
//  NB_C++
//
//  Created by Alexandru Cristian on 06/04/2017.
//
//

#ifndef Class_hpp
#define Class_hpp

#include <fstream>
#include <sstream>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <queue>
#include <utility>
#include <cmath>
#include <cassert>
#include <cstring>
#include <cstdint>
#include <string>
#include <atomic>
#include <algorithm>
#include <boost/filesystem.hpp>
using namespace boost::filesystem;
using namespace std;

class Genome;

template <class T>
class Class {
private:
public:

    /**
     * Instantiates a Class instance.
     * @param id_       The unique identifier for this class.
     * @param kmer_size The kmer size used to create the savefile.
     * @param savefile  The path to the savefile on disk.
     */
    Class(string id_, int kmer_size, path savefile);

    ~Class();

    /**
     * Getter for the class ID.
     * @return A string containing the class unique ID.
     */
    string getId();

    /**
     * Checks that a savefile was written for this k-mer size and counting
     * mode, reading only its header. Exits with an error on mismatch.
     * Legacy savefiles (no header) pass with a warning.
     * @param source_path The path to the savefile on disk.
     */
    void checkSavefileHeader(path source_path);

    /**
     * Readies the instance for use and loads data from specified savefile.
     * @param source_path The path to the savefile on disk.
     */
    void load(path source_path);

    /**
     * Readies the instance for use and loads data from the set savefile.
     */
    void load();

    /**
     * Releases the class's k-mer data. It must be reloaded before further
     * use.
     */
    void unload();

    /**
     * Writes the class to a savefile, folding in any genomes trained since it
     * was loaded.
     * @param destination_path The path to the savefile on disk.
     */
    void save(path destination_path);

    /**
     * Writes the class to the set savefile.
     */
    void save();

    /**
     * Looks up the class's smoothed log count for each k-mer of a sorted,
     * duplicate-free list.
     * @param kmers Sorted, unique k-mers.
     * @param out   Resized to kmers.size(); out[i] is the log count of
     *              kmers[i], or 0 (log 1, the smoothing pseudo-count) if the
     *              class has never seen it.
     */
    void lookupFreqCounts_lg(const vector<T>& kmers, vector<double>& out);

    /**
     * Adds a genome's k-mer counts to the class.
     * @param genome The genome to be added to the class.
     */
    void addGenome(Genome* genome);

    /**
     * Queues a genome to be added to this class.
     * @param genome The genome to be trained on.
     */
    void queueGenome(Genome* genome);

    /**
     * Train this class on all of its queued genomes, waiting to be added.
     * Genomes already trained into the class (same file name) are skipped.
     */
    void addGenomesInQueue();

    /**
     * @return The log of (1 + the number of genomes in this class).
     */
    double getNGenomes_lg();

    /**
     * @return The log of (total k-mer count + vocabulary size), the smoothed
     *         denominator.
     */
    double getSumFreq_lg();

    /**
     * Gets a path to the savefile of this class.
     * @return Path to savefile.
     */
    path getSavefilePath();

    /**
     * @return True if the class is ready for use (has been initialized /
     *         loaded from disk).
     */
    bool loaded();

    /**
     * Whether k-mer counts are canonical (a k-mer and its reverse complement
     * counted together). Sets the vocabulary size used for smoothing.
     */
    static bool CANONICAL_KMERS;

    /**
     * Number of distinct k-mers of the given size, the "K" of Laplace
     * smoothing.
     */
    static double vocabularySize(int kmer_size);

    /** Savefile format written by this build. */
    static const int SAVEFILE_VERSION = 3;

protected:
  bool isLoaded = false;
  string id;
  int kmerSize;
  path savefile;
  vector<Genome*> genomes;
  queue<Genome*> queuedGenomes;

  // File names of genomes already trained into this class (kept in the
  // savefile), so re-running training on the same files doesn't count them
  // twice.
  unordered_set<string> trainedGenomes;

  typedef pair<double, bool> double_wflag;

  // Log data cached for classifying. The flag is false while plain counts
  // added since loading still need folding in.
  double_wflag ngenomes_lg, sumfreq_lg;

  // Plain counts added since loading.
  int ngenomes;
  long long int sumfreq;

  // Trained k-mers (sorted) and their smoothed log counts, log(1 + count).
  vector<T> kmers;
  vector<double> freqcnt_lg;

  // K-mer counts added since loading, sorted by k-mer, not yet folded into
  // freqcnt_lg.
  vector<pair<T, int> > pending;

  void addNGenomes(int ngenomes_);
  void addSumfreq(long long int sumfreq_);
  void foldInPending();

  void readTextBody(std::istream& in, bool legacy);
  void readBinaryBody(std::istream& in);

  [[noreturn]] void savefileError(path file, string what);
  static uint64_t serializeDouble(double value);
  static double deserializeDouble(uint64_t value);
  static double logAdd(double a, double b);
};

#include "Class.cpp"

#endif /* Class_hpp */
