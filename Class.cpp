//
//  Class.cpp
//  NB_C++
//
//  Created by Alexandru Cristian on 06/04/2017.
//
//

#include "Genome.hpp"
#include "NB.hpp"
#include <climits>

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
vector<Genome*> Class<T>::takeGenomesToTrain(){
  vector<Genome*> to_train;
  while(!queuedGenomes.empty()){
    Genome* current = queuedGenomes.front();
    queuedGenomes.pop();
    genomes.push_back(current); // owned by the class from here on
    string name = current->getKmrPath().filename().native();
    if(!trainedGenomes.insert(name).second){
      cout<<"Skipping "<<name<<": already trained into class "<<id<<"\n";
    }else{
      to_train.push_back(current);
    }
  }
  return to_train;
}

template <class T>
void Class<T>::addGenomesInQueue(){
  vector<Genome*> to_train = takeGenomesToTrain();
  vector<pair<T, int> > counts;
  if(!to_train.empty()){
    to_train[0]->loadSortedKmerCounts(counts);
  }
  beginUpdate(1, to_train.empty() ? NULL : &counts);
  for(size_t i = 0; i < to_train.size(); i++){
    if(i > 0){
      to_train[i]->loadSortedKmerCounts(counts);
    }
    addGenomeCounts(counts);
  }
  commitRange(0);
}

template <class T>
size_t Class<T>::numKmers(){
  return kmers.size();
}

template <class T>
uint64_t Class<T>::queuedBytes(){
  uint64_t total = 0;
  queue<Genome*> copy = queuedGenomes;
  while(!copy.empty()){
    boost::system::error_code ec;
    uintmax_t size = file_size(copy.front()->getKmrPath(), ec);
    if(!ec){
      total += size;
    }
    copy.pop();
  }
  return total;
}

template <class T>
size_t Class<T>::rangeStart(const vector<T>& sorted, size_t r){
  // Range r holds k-mers x with splitters[r-1] <= x < splitters[r].
  if(r == 0){
    return 0;
  }
  if(r > splitters.size()){
    return sorted.size();
  }
  return lower_bound(sorted.begin(), sorted.end(), splitters[r-1]) - sorted.begin();
}

template <class T>
void Class<T>::beginUpdate(size_t R, const vector<pair<T, int> >* sample){
  if(R < 1){
    R = 1;
  }
  splitters.clear();
  for(size_t i = 1; i < R; i++){
    if(kmers.size() >= R){
      splitters.push_back(kmers[i * kmers.size() / R]);
    }else if(sample != NULL && sample->size() >= R){
      splitters.push_back((*sample)[i * sample->size() / R].first);
    }
  }
  R = splitters.size() + 1;
  pendingRanges.assign(R, vector<pair<T, int> >());
  rangeLocks.reset(new mutex[R]);
  newKmers.assign(R, vector<T>());
  newFreqcnt_lg.assign(R, vector<double>());
  addedGenomes = 0;
  addedSum = 0;
  updating = true;
}

template <class T>
void Class<T>::addGenomeCounts(const vector<pair<T, int> >& counts){
  long long int sum = 0;
  for(size_t i = 0; i < counts.size(); i++){
    sum += counts[i].second;
  }
  size_t R = pendingRanges.size();
  size_t lo = 0;
  for(size_t r = 0; r < R; r++){
    size_t hi = (r + 1 < R)
      ? lower_bound(counts.begin(), counts.end(), make_pair(splitters[r], INT_MIN))
        - counts.begin()
      : counts.size();
    if(hi > lo){
      // Integer counts add up to the same totals in any order, so genomes
      // can be merged in whatever order threads get to them.
      lock_guard<mutex> guard(rangeLocks[r]);
      vector<pair<T, int> >& pending = pendingRanges[r];
      vector<pair<T, int> > merged;
      merged.reserve(pending.size() + (hi - lo));
      size_t i = 0, j = lo;
      while(i < pending.size() || j < hi){
        if(j == hi || (i < pending.size() && pending[i].first < counts[j].first)){
          merged.push_back(pending[i++]);
        }else if(i < pending.size() && pending[i].first == counts[j].first){
          merged.push_back(make_pair(pending[i].first,
                                     pending[i].second + counts[j].second));
          i++; j++;
        }else{
          merged.push_back(counts[j++]);
        }
      }
      pending.swap(merged);
    }
    lo = hi;
  }
  addedGenomes++;
  addedSum += sum;
}

template <class T>
void Class<T>::commitRange(size_t r){
  // A k-mer's log count is log(1 + total count); the 1 is the smoothing
  // pseudo-count. New counts are folded in with one logAdd per save, exactly
  // like the original lazy update.
  size_t i = rangeStart(kmers, r), iend = rangeStart(kmers, r + 1);
  vector<pair<T, int> >& pending = pendingRanges[r];
  vector<T>& out_k = newKmers[r];
  vector<double>& out_v = newFreqcnt_lg[r];
  out_k.reserve((iend - i) + pending.size());
  out_v.reserve((iend - i) + pending.size());
  size_t j = 0;
  while(i < iend || j < pending.size()){
    if(j == pending.size() || (i < iend && kmers[i] < pending[j].first)){
      out_k.push_back(kmers[i]);
      out_v.push_back(freqcnt_lg[i]);
      i++;
    }else{
      double previous = 0; // log(1): the pseudo-count of an unseen k-mer
      if(i < iend && kmers[i] == pending[j].first){
        previous = freqcnt_lg[i];
        i++;
      }
      out_k.push_back(pending[j].first);
      out_v.push_back(logAdd(previous, log(pending[j].second)));
      j++;
    }
  }
  vector<pair<T, int> >().swap(pending);
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
void Class<T>::lookupFreqCounts_lg(const vector<T>& query, vector<double>& out){
  out.assign(query.size(), 0.0);
  lookupFreqCounts_lg(query, 0, query.size(), out);
}

template <class T>
void Class<T>::lookupFreqCounts_lg(const vector<T>& query, size_t qlo,
                                   size_t qhi, vector<double>& out){
  for(size_t j = qlo; j < qhi; j++){
    out[j] = 0.0;
  }
  if(kmers.empty() || qhi <= qlo){
    return;
  }
  // Only the part of the class between this slice's first and last k-mer
  // can match.
  size_t ilo = lower_bound(kmers.begin(), kmers.end(), query[qlo]) - kmers.begin();
  size_t ihi = upper_bound(kmers.begin() + ilo, kmers.end(), query[qhi-1])
               - kmers.begin();
  size_t n = ihi - ilo, q = qhi - qlo;
  if(n == 0){
    return;
  }
  // Binary-search each query k-mer when that part of the class is much
  // larger than the slice; otherwise walk both sorted lists once.
  if(q * (size_t)log2((double)n + 1) < n){
    for(size_t j = qlo; j < qhi; j++){
      typename vector<T>::const_iterator it =
        lower_bound(kmers.begin() + ilo, kmers.begin() + ihi, query[j]);
      if(it != kmers.begin() + ihi && *it == query[j]){
        out[j] = freqcnt_lg[it - kmers.begin()];
      }
    }
  }else{
    size_t i = ilo;
    for(size_t j = qlo; j < qhi && i < ihi; j++){
      while(i < ihi && kmers[i] < query[j]){
        i++;
      }
      if(i < ihi && kmers[i] == query[j]){
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
  if(updating){
    // Same totals as adding the genomes one by one.
    if(addedGenomes > 0){
      addNGenomes(addedGenomes);
      addSumfreq(addedSum);
    }
  }else{
    // Nothing trained since loading: write the class as it is.
    newKmers.assign(1, vector<T>());
    newFreqcnt_lg.assign(1, vector<double>());
    newKmers[0].swap(kmers);
    newFreqcnt_lg[0].swap(freqcnt_lg);
  }
  updating = false;
  double ngenomes_v = getNGenomes_lg(), sumfreq_v = getSumFreq_lg();
  uint64_t n = 0, m = trainedGenomes.size();
  for(size_t r = 0; r < newKmers.size(); r++){
    n += newKmers[r].size();
  }

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
    // Ranges are in k-mer order, so writing them one after another gives the
    // whole sorted array.
    for(size_t r = 0; r < newKmers.size(); r++){
      out.write((const char*)newKmers[r].data(), newKmers[r].size() * sizeof(T));
    }
    for(size_t r = 0; r < newFreqcnt_lg.size(); r++){
      out.write((const char*)newFreqcnt_lg[r].data(),
                newFreqcnt_lg[r].size() * sizeof(double));
    }
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
  // The class is unloaded after saving; drop the written arrays now.
  vector<vector<T> >().swap(newKmers);
  vector<vector<double> >().swap(newFreqcnt_lg);
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
  vector<vector<pair<T, int> > >().swap(pendingRanges);
  vector<vector<T> >().swap(newKmers);
  vector<vector<double> >().swap(newFreqcnt_lg);
  splitters.clear();
  updating = false;
  trainedGenomes.clear();

  isLoaded = false;
}

template <class T>
bool Class<T>::loaded(){
  return isLoaded;
}
