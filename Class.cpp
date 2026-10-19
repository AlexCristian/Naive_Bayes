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
  for(typename vector<Genome*>::iterator iter = genomes.begin();
      iter != genomes.end(); iter++){
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
double Class<T>::logAdd(double a, double b){
  // log(e^a + e^b), computed exactly as the original vector version did
  // (same max choice and summation order) so savefiles stay bit-identical.
  double max = a;
  if(max < b){
    max = b;
  }
  double result = 0;
  result += exp(a - max);
  result += exp(b - max);
  return max + log(result);
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
path Class<T>::getSavefilePath(){
  return savefile;
}

template <class T>
double Class<T>::getNGenomes_lg(){
  if(!ngenomes_lg.second){
    ngenomes_lg.first = logAdd(ngenomes_lg.first, log(ngenomes));
    ngenomes = 0;
    ngenomes_lg.second = true;
  }
  return ngenomes_lg.first;
}

template <class T>
double Class<T>::getSumFreq_lg(){
  if(!sumfreq_lg.second){
    sumfreq_lg.first = logAdd(sumfreq_lg.first, log(sumfreq));
    sumfreq = 0;
    sumfreq_lg.second = true;
  }
  return sumfreq_lg.first;
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
void Class<T>::addGenome(Genome* genome){
  addNGenomes(1);

  vector<pair<T, int> > counts;
  genome->loadSortedKmerCounts(counts);

  // Merge into the pending counts (both sorted by k-mer).
  long long int sum = 0;
  vector<pair<T, int> > merged;
  merged.reserve(pending.size() + counts.size());
  size_t i = 0, j = 0;
  while(i < pending.size() || j < counts.size()){
    if(j == counts.size()
       || (i < pending.size() && pending[i].first < counts[j].first)){
      merged.push_back(pending[i++]);
    }else{
      sum += counts[j].second;
      if(i < pending.size() && pending[i].first == counts[j].first){
        merged.push_back(make_pair(pending[i].first,
                                   pending[i].second + counts[j].second));
        i++;
      }else{
        merged.push_back(counts[j]);
      }
      j++;
    }
  }
  pending.swap(merged);

  addSumfreq(sum);
  genomes.push_back(genome);
}

template <class T>
void Class<T>::foldInPending(){
  if(pending.empty()){
    return;
  }
  // A k-mer's log count is log(1 + total count); the 1 is the smoothing
  // pseudo-count. New counts are folded in with one logAdd per save, exactly
  // like the original lazy update.
  vector<T> new_kmers;
  vector<double> new_lg;
  new_kmers.reserve(kmers.size() + pending.size());
  new_lg.reserve(kmers.size() + pending.size());
  size_t i = 0, j = 0;
  while(i < kmers.size() || j < pending.size()){
    if(j == pending.size()
       || (i < kmers.size() && kmers[i] < pending[j].first)){
      new_kmers.push_back(kmers[i]);
      new_lg.push_back(freqcnt_lg[i]);
      i++;
    }else{
      double previous = 0; // log(1): the pseudo-count of an unseen k-mer
      if(i < kmers.size() && kmers[i] == pending[j].first){
        previous = freqcnt_lg[i];
        i++;
      }
      new_kmers.push_back(pending[j].first);
      new_lg.push_back(logAdd(previous, log(pending[j].second)));
      j++;
    }
  }
  kmers.swap(new_kmers);
  freqcnt_lg.swap(new_lg);
  vector<pair<T, int> >().swap(pending);
}

template <class T>
void Class<T>::lookupFreqCounts_lg(const vector<T>& query, vector<double>& out){
  out.assign(query.size(), 0.0);
  size_t n = kmers.size(), q = query.size();
  if(n == 0 || q == 0){
    return;
  }
  // Binary-search each query k-mer when the class is much larger than the
  // query; otherwise walk both sorted lists once.
  if(q * (size_t)log2((double)n + 1) < n){
    for(size_t j = 0; j < q; j++){
      typename vector<T>::const_iterator it =
        lower_bound(kmers.begin(), kmers.end(), query[j]);
      if(it != kmers.end() && *it == query[j]){
        out[j] = freqcnt_lg[it - kmers.begin()];
      }
    }
  }else{
    size_t i = 0;
    for(size_t j = 0; j < q && i < n; j++){
      while(i < n && kmers[i] < query[j]){
        i++;
      }
      if(i < n && kmers[i] == query[j]){
        out[j] = freqcnt_lg[i];
      }
    }
  }
}

template <class T>
void Class<T>::checkSavefileHeader(path source_path){
  std::ifstream in(source_path.native(), ios::binary);
  if(!in.good()){
    savefileError(source_path, "cannot open for reading");
  }
  string magic;
  in>>magic;
  if(magic != "NBSAVE"){
    return; // legacy savefile; load() warns
  }
  int version, k, canonical;
  if(!(in>>version>>k>>canonical)){
    savefileError(source_path, "corrupt header");
  }
  if(version != 2 && version != SAVEFILE_VERSION){
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
void Class<T>::load(path source_path){
  if(isLoaded){
    return;
  }else{
    isLoaded = true;
  }

  trainedGenomes.clear();

  // No savefile yet means a new, empty class. An existing savefile that
  // can't be read is an error, not an empty class.
  if(!exists(source_path)){
    return;
  }
  std::ifstream in(source_path.native(), ios::binary);
  if(!in.good()){
    savefileError(source_path, "cannot open for reading");
  }

  string magic;
  in>>magic;
  if(magic != "NBSAVE"){
    // Format 1: no header, no genome list, nothing to check k against.
    static atomic<bool> warned(false);
    if(!warned.exchange(true)){
      cerr<<"[WARNING] Legacy savefile format (e.g. "<<source_path.native()
          <<"): k-mer size and counting mode can't be verified, and "
          "already-trained genomes can't be recognized. Savefiles are "
          "upgraded when next trained.\n";
    }
    in.clear();
    in.seekg(0);
    readTextBody(in, true);
    return;
  }

  int version, k, canonical;
  in>>version>>k>>canonical;
  if(!in || (version != 2 && version != SAVEFILE_VERSION) || k != kmerSize
     || (canonical != 0) != CANONICAL_KMERS){
    // checkSavefileHeader() reports mismatches with a clear message.
    checkSavefileHeader(source_path);
    savefileError(source_path, "corrupt header");
  }
  string rest;
  getline(in, rest); // end of the header line
  if(version == SAVEFILE_VERSION){
    readBinaryBody(in);
  }else{
    readTextBody(in, false);
  }
}

template <class T>
void Class<T>::readTextBody(std::istream& in, bool legacy){
  // Formats 1 and 2: decimal text, k-mers in hash order.
  uint64_t tmp;
  if(!(in>>tmp)){
    savefileError(savefile, "truncated or corrupt (header values)");
  }
  ngenomes_lg = make_pair(deserializeDouble(tmp), true);
  in>>tmp;
  sumfreq_lg = make_pair(deserializeDouble(tmp), true);

  size_t n;
  in>>n;
  if(!in){
    savefileError(savefile, "truncated or corrupt (header values)");
  }

  vector<pair<T, double> > entries;
  entries.reserve(n);
  T kmer;
  for(size_t i=0; i<n; i++){
    if(!(in>>kmer>>tmp)){
      savefileError(savefile, "truncated or corrupt (read " + to_string(i)
                    + " of " + to_string(n) + " k-mers)");
    }
    entries.push_back(make_pair(kmer, deserializeDouble(tmp)));
  }
  sort(entries.begin(), entries.end());
  kmers.resize(n);
  freqcnt_lg.resize(n);
  for(size_t i=0; i<n; i++){
    kmers[i] = entries[i].first;
    freqcnt_lg[i] = entries[i].second;
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
void Class<T>::readBinaryBody(std::istream& in){
  double ngenomes_v, sumfreq_v;
  uint64_t n;
  in.read((char*)&ngenomes_v, sizeof(ngenomes_v));
  in.read((char*)&sumfreq_v, sizeof(sumfreq_v));
  in.read((char*)&n, sizeof(n));
  if(!in){
    savefileError(savefile, "truncated or corrupt (header values)");
  }
  ngenomes_lg = make_pair(ngenomes_v, true);
  sumfreq_lg = make_pair(sumfreq_v, true);

  kmers.resize(n);
  freqcnt_lg.resize(n);
  in.read((char*)kmers.data(), n * sizeof(T));
  in.read((char*)freqcnt_lg.data(), n * sizeof(double));
  if(!in){
    savefileError(savefile, "truncated or corrupt (k-mer table)");
  }

  uint64_t m;
  in.read((char*)&m, sizeof(m));
  if(!in){
    savefileError(savefile, "truncated or corrupt (genome list missing)");
  }
  string name;
  for(uint64_t i=0; i<m; i++){
    if(!getline(in, name)){
      savefileError(savefile, "truncated or corrupt (read " + to_string(i)
                    + " of " + to_string(m) + " genome names)");
    }
    trainedGenomes.insert(name);
  }
}

template <class T>
void Class<T>::save(path destination_path){
  foldInPending();
  double ngenomes_v = getNGenomes_lg(), sumfreq_v = getSumFreq_lg();
  uint64_t n = kmers.size(), m = trainedGenomes.size();

  /*
   * Savefile format (version 3):
   * NBSAVE <version> <k-mer size> <canonical 0/1>\n      (text)
   * ngenomes_lg sumfreq_lg n                             (double, double, uint64)
   * kmers[n] (sorted)  freqcnt_lg[n]                     (k-mer type, double)
   * m, then m genome file names, one per line
   * Binary fields use the machine's byte order (little-endian on x86/ARM).
   *
   * Written to a temp file and renamed over the old savefile, so a crash or
   * full disk mid-write never leaves a truncated savefile behind.
   */
  path tmp_path = path(destination_path.native() + ".tmp");
  {
    std::ofstream out(tmp_path.native(), ios::binary);
    out<<"NBSAVE "<<SAVEFILE_VERSION<<" "<<kmerSize<<" "<<CANONICAL_KMERS<<"\n";
    out.write((const char*)&ngenomes_v, sizeof(ngenomes_v));
    out.write((const char*)&sumfreq_v, sizeof(sumfreq_v));
    out.write((const char*)&n, sizeof(n));
    out.write((const char*)kmers.data(), n * sizeof(T));
    out.write((const char*)freqcnt_lg.data(), n * sizeof(double));
    out.write((const char*)&m, sizeof(m));
    for(typename unordered_set<string>::iterator iter = trainedGenomes.begin();
      iter != trainedGenomes.end(); iter++){
        out<<*iter<<"\n";
      }
    out.close();
    if(out.fail()){
      remove(tmp_path);
      savefileError(destination_path, "write failed (disk full?)");
    }
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
  vector<T>().swap(kmers);
  vector<double>().swap(freqcnt_lg);
  vector<pair<T, int> >().swap(pending);
  trainedGenomes.clear();

  isLoaded = false;
}

template <class T>
bool Class<T>::loaded(){
  return isLoaded;
}
