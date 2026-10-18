//
//  Class.cpp
//  NB_C++
//
//  Created by Alexandru Cristian on 06/04/2017.
//
//

#include "Genome.hpp"
#include "NB.hpp"

template <class T>
bool Class<T>::CANONICAL_KMERS = true;

template <class T>
double Class<T>::vocabularySize(int kmer_size){
  double all = pow(4.0, kmer_size);
  if(!CANONICAL_KMERS){
    return all;
  }
  // A k-mer and its reverse complement share one canonical form; only even-k
  // palindromes are their own reverse complement.
  double palindromes = (kmer_size % 2 == 0) ? pow(4.0, kmer_size / 2) : 0;
  return (all + palindromes) / 2;
}

template <class T>
Class<T>::Class(string id_, int kmer_size, path _savefile)  {
  ngenomes_lg = make_pair(0, false);

  // Laplace smoothing adds 1 per possible k-mer, so the denominator starts at
  // the vocabulary size and grows with the class's total k-mer count.
  sumfreq_lg = make_pair(log(vocabularySize(kmer_size)), false);

  ngenomes=0;
  sumfreq=0;

  id = id_;
  kmerSize = kmer_size;

  savefile = _savefile;
}

template <class T>
void Class<T>::savefileError(path file, string what){
  cerr<<"[ERROR] Savefile "<<file.native()<<" (class "<<id<<"): "<<what<<"\n";
  cerr.flush();
  // Other threads may be mid-update; their savefiles are written to a temp
  // file and renamed, so stopping here leaves every savefile intact.
  _Exit(1);
}

template <class T>
Class<T>::~Class() {
  for(vector<Genome*>::iterator iter = genomes.begin(); iter != genomes.end();
      iter++){
        delete *iter;
      }
}

template <class T>
void Class<T>::queueGenome(Genome* genome){
  queuedGenomes.push(genome);
}

template <class T>
void Class<T>::addGenomesInQueue(){
  while(!queuedGenomes.empty()){
    Genome* current = queuedGenomes.front();
    string name = current->getKmrPath().filename().native();
    if(!trainedGenomes.insert(name).second){
      cout<<"Skipping "<<name<<": already trained into class "<<id<<"\n";
      genomes.push_back(current);
    }else{
      addGenome(current);
    }
    current->unload();
    queuedGenomes.pop();
  }
}

template <class T>
double Class<T>::logAdd(vector<double> exponents){
  double result = 0, max;
  int maxId = -1;
  for(vector<double>::iterator num=exponents.begin();
    num != exponents.end(); num++){
    if(maxId == -1 || max < *num){
      maxId = 0;
      max = *num;
    }
  }

  for(vector<double>::iterator num=exponents.begin();
    num != exponents.end(); num++){
    result += exp(*num - max);
  }
  result = max + log(result);

  return result;
}

template <class T>
void Class<T>::addNGenomes(int ngenomes_){
  if(ngenomes_lg.second){
    ngenomes = 0;
    ngenomes_lg.second = false;
  }
  ngenomes += ngenomes_;
}

template <class T>
void Class<T>::addSumfreq(long long int sumfreq_){
  if(sumfreq_lg.second){
    sumfreq = 0;
    sumfreq_lg.second = false;
  }
  sumfreq += sumfreq_;
}

template <class T>
void Class<T>::addFreqCount(T feature_, int count_){
  if(getFreqcnt().find(feature_) == getFreqcnt().end()){
    getFreqcnt()[feature_] = count_;
  }else{
    getFreqcnt()[feature_] += count_;
  }
  if(getFreqcnt_lg().find(feature_) == getFreqcnt_lg().end()){
      getFreqcnt_lg()[feature_] = make_pair(0 ,false);
  }else{
      getFreqcnt_lg()[feature_].second = false;
  }
}

template <class T>
void Class<T>::addSumfreq_lg(double sumfreq_){
  vector<double> items;
  items.push_back(getSumFreq_lg());
  items.push_back(sumfreq_);
  sumfreq_lg.first = logAdd(items);
}

template <class T>
void Class<T>::addNGenomes_lg(double ngenomes_){
  vector<double> items;
  items.push_back(getNGenomes_lg());
  items.push_back(ngenomes_);
  ngenomes_lg.first = logAdd(items);
}

template <class T>
void Class<T>::addFreqCount_lg(T feature_, double count_){
  vector<double> items;
  items.push_back(getFreqCount_lg(feature_));
  items.push_back(count_);
  getFreqcnt_lg()[feature_].first = logAdd(items);
}

template <class T>
path Class<T>::getSavefilePath(){
  return savefile;
}

template <class T>
double Class<T>::getNGenomes_lg(){
  if(!ngenomes_lg.second){
    vector<double> terms;
    terms.push_back(ngenomes_lg.first);
    terms.push_back(log(ngenomes));
    ngenomes = 0;
    ngenomes_lg.first = logAdd(terms);
    ngenomes_lg.second = true;
  }
  return ngenomes_lg.first;
}

template <class T>
double Class<T>::getFreqCount_lg(T feature_){
  if(!existsInFreqCount(feature_)){
    return 0;
  }else{
    if(!getFreqcnt_lg()[feature_].second){
      vector<double> terms;
      terms.push_back(getFreqcnt_lg()[feature_].first);
      terms.push_back(log(getFreqcnt()[feature_]));
      getFreqcnt()[feature_] = 0;
      getFreqcnt_lg()[feature_].first = logAdd(terms);
      getFreqcnt_lg()[feature_].second = true;
    }
    return getFreqcnt_lg()[feature_].first;
  }
}

template <class T>
double Class<T>::getSumFreq_lg(){
  if(!sumfreq_lg.second){
    vector<double> terms;
    terms.push_back(sumfreq_lg.first);
    terms.push_back(log(sumfreq));
    sumfreq = 0;
    sumfreq_lg.first = logAdd(terms);
    sumfreq_lg.second = true;
  }
  return sumfreq_lg.first;
}

template <class T>
bool Class<T>::existsInFreqCount(T feature_){
  return (getFreqcnt().find(feature_) != getFreqcnt().end())
    || (getFreqcnt_lg().find(feature_) != getFreqcnt_lg().end());
}

template <class T>
string Class<T>::getId(){
  return id;
}

template <class T>
uint64_t Class<T>::serializeDouble(double value){
  uint64_t result;
  assert(sizeof(result) == sizeof(value));
  memcpy(&result, &value, sizeof(result));
  return result;
}

template <class T>
double Class<T>::deserializeDouble(uint64_t value){
  double result;
  assert(sizeof(result) == sizeof(value));
  memcpy(&result, &value, sizeof(result));
  return result;
}

template <class T>
void Class<T>::addGenome(Genome* genome, double confidence_lg){
  if(confidence_lg == 0.0)
    addNGenomes(1);
  else
    addNGenomes_lg(confidence_lg + log(1));

  long long int sum=0, new_features=0;
  unordered_map<T, int>& freqs = genome->getKmerCounts();
  for(typename unordered_map<T, int>::iterator iter = freqs.begin();
    iter != freqs.end(); iter++){
      sum += iter->second;
      if(existsInFreqCount(iter->first)){
        new_features++;
      }
      if(confidence_lg == 0.0)
        addFreqCount(iter->first, iter->second);
      else
        addFreqCount_lg(iter->first, confidence_lg + log(iter->second));
    }

  if(confidence_lg == 0.0)
    addSumfreq(sum);
  else
    addSumfreq_lg(confidence_lg + log(sum));
  genomes.push_back(genome);
}

template <class T>
void Class<T>::computeBatchNumerators(vector<Genome*> genomes){
  for(vector<Genome*>::iterator iter = genomes.begin(); iter != genomes.end();
    iter++){
      (*iter)->computeClassificationNumerator(this);
  }
}

template <class T>
void Class<T>::load(path source_path){
  if(isLoaded){
    return;
  }else{
    isLoaded = true;
  }

  freqcnt_lg = new unordered_map<T, double_wflag>;
  freqcnt = new unordered_map<T, int>;
  trainedGenomes.clear();

  // No savefile yet means a new, empty class. An existing savefile that
  // can't be read is an error, not an empty class.
  if(!exists(source_path)){
    return;
  }
  std::ifstream in(source_path.native());
  if(!in.good()){
    savefileError(source_path, "cannot open for reading");
  }
  deserialize(in);
  in.close();
}

template <class T>
void Class<T>::save(path destination_path){
  // Write to a temp file and rename it over the old savefile, so a crash or
  // full disk mid-write never leaves a truncated savefile behind.
  path tmp_path = path(destination_path.native() + ".tmp");
  std::ofstream out(tmp_path.native());
  out<<serialize();
  out.close();
  if(out.fail()){
    remove(tmp_path);
    savefileError(destination_path, "write failed (disk full?)");
  }
  boost::system::error_code ec;
  rename(tmp_path, destination_path, ec);
  if(ec){
    savefileError(destination_path, "could not replace savefile: " + ec.message());
  }
}

template <class T>
void Class<T>::save(){
  save(savefile);
}

template <class T>
void Class<T>::load(){
  load(savefile);
}

template <class T>
void Class<T>::unload(){
  delete freqcnt;
  delete freqcnt_lg;
  trainedGenomes.clear();

  isLoaded = false;
}

template <class T>
bool Class<T>::loaded(){
  return isLoaded;
}

template <class T>
unordered_map<T, pair<double, bool> >& Class<T>::getFreqcnt_lg(){
  return *freqcnt_lg;
}

template <class T>
unordered_map<T, int>& Class<T>::getFreqcnt(){
  return *freqcnt;
}

/*
 *  WARNING: Serialization may be platform-dependent, as C++11 does not
 *  guarantee all 'double' types to follow the IEEE-754 64-bit standard.
 */
template <class T>
string Class<T>::serialize(){
  ostringstream strs;
  /*
   * Serialization format spec (version 2):
   * NBSAVE <version> <k-mer size> <canonical 0/1>
   * ngenomes_lg sumfreq_lg freqcnt_lg.size() [kmer freqcnt_lg]...
   * GENOMES <count>
   * <one trained genome file name per line>
   * Doubles are stored as their IEEE-754 bit patterns (uint64).
   */

  strs<<"NBSAVE "<<SAVEFILE_VERSION<<" "<<kmerSize<<" "<<CANONICAL_KMERS<<"\n";
  strs<<serializeDouble(getNGenomes_lg())<<" ";
  strs<<serializeDouble(getSumFreq_lg())<<" "<<getFreqcnt_lg().size();
  for(typename unordered_map<T, double_wflag>::iterator iter = getFreqcnt_lg().begin();
    iter != getFreqcnt_lg().end(); iter++){
      strs<<" "<<iter->first<<" "<<serializeDouble(getFreqCount_lg(iter->first));
    }
  strs<<"\nGENOMES "<<trainedGenomes.size()<<"\n";
  for(typename unordered_set<string>::iterator iter = trainedGenomes.begin();
    iter != trainedGenomes.end(); iter++){
      strs<<*iter<<"\n";
    }

  return strs.str();
}

template <class T>
void Class<T>::checkSavefileHeader(path source_path){
  std::ifstream in(source_path.native());
  if(!in.good()){
    savefileError(source_path, "cannot open for reading");
  }
  string magic;
  in>>magic;
  if(magic != "NBSAVE"){
    return; // legacy savefile; deserialize() warns
  }
  int version, k, canonical;
  if(!(in>>version>>k>>canonical)){
    savefileError(source_path, "corrupt header");
  }
  if(version != SAVEFILE_VERSION){
    savefileError(source_path, "unsupported format version "
                  + to_string(version));
  }
  if(k != kmerSize){
    savefileError(source_path, "trained with k=" + to_string(k)
                  + " but running with k=" + to_string(kmerSize));
  }
  if((canonical != 0) != CANONICAL_KMERS){
    savefileError(source_path, string("trained with ")
                  + (canonical ? "canonical" : "non-canonical")
                  + " k-mers, but --canonical is "
                  + (CANONICAL_KMERS ? "1" : "0"));
  }
}

template <class T>
void Class<T>::deserialize(std::ifstream& in){
  uint64_t tmp;
  string first;
  bool legacy = false;

  in>>first;
  if(first == "NBSAVE"){
    int version, k, canonical;
    in>>version>>k>>canonical;
    // checkSavefileHeader() reports mismatches with a clear message.
    if(!in || version != SAVEFILE_VERSION || k != kmerSize
       || (canonical != 0) != CANONICAL_KMERS){
      checkSavefileHeader(savefile);
      savefileError(savefile, "corrupt header");
    }
    in>>tmp;
  }else{
    // Format 1: no header, no genome list, nothing to check k against.
    legacy = true;
    static atomic<bool> warned(false);
    if(!warned.exchange(true)){
      cerr<<"[WARNING] Legacy savefile format (e.g. "<<savefile.native()
          <<"): k-mer size and counting mode can't be verified, and "
          "already-trained genomes can't be recognized. Savefiles are "
          "upgraded when next trained.\n";
    }
    try{
      tmp = stoull(first);
    }catch(...){
      savefileError(savefile, "not a savefile");
    }
  }
  ngenomes_lg = make_pair(deserializeDouble(tmp), true);

  in>>tmp;
  sumfreq_lg = make_pair(deserializeDouble(tmp), true);

  size_t n;
  T kmer;
  in>>n;
  if(!in){
    savefileError(savefile, "truncated or corrupt (header values)");
  }

  for(size_t i=0; i<n; i++){
    if(!(in>>kmer>>tmp)){
      savefileError(savefile, "truncated or corrupt (read " + to_string(i)
                    + " of " + to_string(n) + " k-mers)");
    }
    getFreqcnt_lg()[kmer] = make_pair(deserializeDouble(tmp), true);
  }

  if(legacy){
    return;
  }
  string tag;
  size_t m;
  if(!(in>>tag>>m) || tag != "GENOMES"){
    savefileError(savefile, "truncated or corrupt (genome list missing)");
  }
  string name;
  getline(in, name); // rest of the GENOMES line
  for(size_t i=0; i<m; i++){
    if(!getline(in, name)){
      savefileError(savefile, "truncated or corrupt (read " + to_string(i)
                    + " of " + to_string(m) + " genome names)");
    }
    trainedGenomes.insert(name);
  }
}

template <class T>
ostream& operator<<(ostream& out, Class<T>& cls){
  out<<cls.serialize();
  return out;
}
