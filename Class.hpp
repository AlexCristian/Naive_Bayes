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
#include <memory>
#include <mutex>
#include <condition_variable>
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
     * Like lookupFreqCounts_lg(), for the slice [qlo, qhi) of the query list
     * only, writing out[qlo..qhi). Threads may run disjoint slices at once.
     */
    void lookupFreqCounts_lg(const vector<T>& kmers, size_t qlo, size_t qhi,
                             vector<double>& out);

    /*
     * Parallel training of a loaded class:
     *   takeGenomesToTrain()      (one thread)
     *   beginUpdate(R, sample)    (one thread)
     *   addGenomeCounts(run)      (any threads, once per genome)
     *   commitRange(r)            (any threads, once per r in [0, R))
     *   save()                    (one thread)
     * Counts are split by k-mer range, so threads never touch the same
     * k-mer, and the saved values are bit-identical to adding the genomes
     * one at a time.
     */

    /**
     * Removes the queued genomes and returns those not yet trained into the
     * class, recording their names. Already-trained ones are skipped.
     */
    vector<Genome*> takeGenomesToTrain();

    /**
     * Splits the k-mer space into R ranges of about equal size, using the
     * class's own k-mers or, for a new class, a sample genome's.
     */
    void beginUpdate(size_t R, const vector<pair<T, int> >* sample);

    /**
     * Adds one genome's sorted k-mer counts. Thread-safe.
     */
    void addGenomeCounts(const vector<pair<T, int> >& counts);

    /**
     * Folds range r's new counts into the class. Thread-safe for distinct r.
     */
    void commitRange(size_t r);

    /** @return Number of k-mer ranges set by beginUpdate(). */
    size_t numRanges(){ return pendingRanges.size(); }

    /** @return Number of k-mers in the class. */
    size_t numKmers();

    /** @return Total size in bytes of the queued genomes' k-mer files. */
    uint64_t queuedBytes();

    /** @return Rough number of k-mers in the queued genomes' files. */
    uint64_t queuedKmerEstimate();

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

    /**
     * How many savefiles may be read or written at once (0: no limit). On a
     * spinning disk, several threads streaming different files at once make
     * it seek between them; one or two streams keep it reading sequentially
     * while other threads compute.
     */
    static void setIoSlots(int slots);

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

  // Parallel update state (see beginUpdate()): range boundaries, the new
  // counts per range, a lock per range, and the committed result per range.
  vector<T> splitters;
  vector<vector<pair<T, int> > > pendingRanges;
  unique_ptr<mutex[]> rangeLocks;
  vector<vector<T> > newKmers;
  vector<vector<double> > newFreqcnt_lg;
  atomic<int> addedGenomes{0};
  atomic<long long> addedSum{0};
  bool updating = false;

  void addNGenomes(int ngenomes_);
  void addSumfreq(long long int sumfreq_);
  size_t rangeStart(const vector<T>& sorted, size_t r);

  void readTextBody(std::istream& in, bool legacy);
  void readBinaryBody(std::istream& in);

  [[noreturn]] void savefileError(path file, string what);

  // Counting semaphore limiting concurrent savefile reads/writes.
  struct IoSlot {
    IoSlot();
    ~IoSlot();
  };
  static int io_slots, io_busy;
  static mutex io_lock;
  static condition_variable io_free;
  static uint64_t serializeDouble(double value);
  static double deserializeDouble(uint64_t value);
  static double logAdd(double a, double b);
};

#include "Class.cpp"

#endif /* Class_hpp */
