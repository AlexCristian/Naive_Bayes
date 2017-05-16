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
#include <vector>
#include <queue>
#include <utility>
#include <cmath>
#include <cassert>
#include <cstring>
#include <boost/filesystem.hpp>
using namespace boost::filesystem;
using namespace std;

class Genome;

template <class T>
class Class {
private:
public:
    Class(string id_, int kmer_size, path savefile);
    ~Class();

    string getId();

    string serialize();
    void deserialize(ifstream& in);

    void load(path source_path);
    void load();
    void unload();
    void save(path destination_path);
    void save();

    void computeBatchNumerators(vector<Genome*> genomes);
    void addGenome(Genome* genome, double confidence_lg = 0.0);
    void queueGenome(Genome* genome);
    void addGenomesInQueue();

    double getNGenomes_lg();
    double getFreqCount_lg(T feature_);
    double getSumFreq_lg();

protected:
  bool isLoaded = false;
  string id;
  path savefile;
  vector<Genome*> genomes;
  queue<Genome*> queuedGenomes;

  typedef pair<double, bool> double_wflag;

  // Log data cached for classifying
  double_wflag ngenomes_lg, sumfreq_lg;
  unordered_map<T, double_wflag> *freqcnt_lg;

  // Plain data, for updates
  int ngenomes;
  long long int sumfreq;
  unordered_map<T, int> *freqcnt;

  void addNGenomes_lg(double ngenomes_);
  void addSumfreq_lg(double sumfreq_);
  void addFreqCount_lg(T feature_, double count_);
  void addNGenomes(int ngenomes_);
  void addSumfreq(long long int sumfreq_);
  void addFreqCount(T feature_, int count_);
  bool existsInFreqCount(T feature_);

  unordered_map<T, double_wflag>& getFreqcnt_lg();
  unordered_map<T, int>& getFreqcnt();

  static uint64_t serializeDouble(double value);
  static double deserializeDouble(uint64_t value);
  static double logAdd(vector<double> exponents);
};

template <class T>
ostream& operator<<(ostream& out, Class<T>& cls);

#include "Class.cpp"

#endif /* Class_hpp */
