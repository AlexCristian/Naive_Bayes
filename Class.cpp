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

template <class T> int Class<T>::io_slots = 0;
template <class T> int Class<T>::io_busy = 0;
template <class T> mutex Class<T>::io_lock;
template <class T> condition_variable Class<T>::io_free;

template <class T>
void Class<T>::setIoSlots(int slots){
  io_slots = slots;
}

template <class T>
Class<T>::IoSlot::IoSlot(){
  unique_lock<mutex> lk(io_lock);
  io_free.wait(lk, [](){ return io_slots <= 0 || io_busy < io_slots; });
  io_busy++;
}

template <class T>
Class<T>::IoSlot::~IoSlot(){
  {
    lock_guard<mutex> lk(io_lock);
    io_busy--;
  }
  io_free.notify_one();
}

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
uint64_t Class<T>::queuedKmerEstimate(){
  uint64_t total = 0;
  queue<Genome*> copy = queuedGenomes;
  while(!copy.empty()){
    total += copy.front()->estimatedKmers();
    copy.pop();
  }
  return total;
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
  if(version < 2 || version > SAVEFILE_VERSION){
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
  IoSlot slot;
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
  if(!in || version < 2 || version > SAVEFILE_VERSION || k != kmerSize
     || (canonical != 0) != CANONICAL_KMERS){
    // checkSavefileHeader() reports mismatches with a clear message.
    checkSavefileHeader(source_path);
    savefileError(source_path, "corrupt header");
  }
  string rest;
  getline(in, rest); // end of the header line
  if(version == 4){
    readCompactBody(in);
  }else if(version == 3){
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

  readGenomeNames(in);
}

template <class T>
void Class<T>::readGenomeNames(std::istream& in){
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
void Class<T>::readCompactBody(std::istream& in){
  // See save() for the layout.
  double ngenomes_v, sumfreq_v;
  uint64_t n, ndict, key_bytes;
  uint8_t width;
  in.read((char*)&ngenomes_v, sizeof(ngenomes_v));
  in.read((char*)&sumfreq_v, sizeof(sumfreq_v));
  in.read((char*)&n, sizeof(n));
  in.read((char*)&ndict, sizeof(ndict));
  if(!in || ndict > n + 1){
    savefileError(savefile, "truncated or corrupt (header values)");
  }
  ngenomes_lg = make_pair(ngenomes_v, true);
  sumfreq_lg = make_pair(sumfreq_v, true);

  vector<double> dict(ndict);
  in.read((char*)dict.data(), ndict * sizeof(double));
  in.read((char*)&width, sizeof(width));
  in.read((char*)&key_bytes, sizeof(key_bytes));
  if(!in || (width != 1 && width != 2 && width != 4)){
    savefileError(savefile, "truncated or corrupt (value table)");
  }
  vector<unsigned char> keys(key_bytes);
  in.read((char*)keys.data(), key_bytes);
  vector<unsigned char> idx(n * width);
  in.read((char*)idx.data(), n * width);
  if(!in){
    savefileError(savefile, "truncated or corrupt (k-mer table)");
  }

  kmers.resize(n);
  freqcnt_lg.resize(n);
  const unsigned char *p = keys.data(), *end = p + key_bytes;
  T prev = 0;
  for(uint64_t i = 0; i < n; i++){
    // Gap to the previous k-mer, 7 bits per byte, low bits first.
    T gap = 0;
    int shift = 0;
    while(true){
      if(p == end || shift > 63){
        savefileError(savefile, "truncated or corrupt (k-mer stream)");
      }
      unsigned char b = *p++;
      gap |= (T)(b & 0x7f) << shift;
      if(!(b & 0x80)){
        break;
      }
      shift += 7;
    }
    prev += gap;
    kmers[i] = prev;
    uint32_t j = 0;
    memcpy(&j, &idx[i * width], width); // little-endian
    if(j >= ndict){
      savefileError(savefile, "truncated or corrupt (value index)");
    }
    freqcnt_lg[i] = dict[j];
  }
  readGenomeNames(in);
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
   * Savefile format (version 4):
   * NBSAVE <version> <k-mer size> <canonical 0/1>\n      (text)
   * ngenomes_lg sumfreq_lg n ndict                       (double, double, uint64, uint64)
   * dict[ndict]                                          (double: distinct log counts)
   * width key_bytes                                      (uint8: 1, 2 or 4; uint64)
   * key stream: gap from the previous k-mer, 7 bits per byte, low first
   * index[n]: position of each k-mer's log count in dict (width bytes each)
   * m, then m genome file names, one per line
   * Binary fields are little-endian (the machine's order on x86/ARM).
   * Version 3 (still readable) stored kmers[n] and freqcnt_lg[n] raw.
   *
   * Written to a temp file and renamed over the old savefile, so a crash or
   * full disk mid-write never leaves a truncated savefile behind.
   */
  // Encode before taking an I/O slot, so only the write waits on the disk.
  // K-mers are sorted, so each is stored as its gap to the previous one, 7
  // bits per byte. Log counts take few distinct values (one per total
  // count), so each is stored as an index into a table of exact doubles.
  vector<unsigned char> keys;
  keys.reserve(n * 3);
  vector<double> dict;
  vector<uint32_t> idx;
  idx.reserve(n);
  {
    unordered_map<uint64_t, uint32_t> where;
    T prev = 0;
    for(size_t r = 0; r < newKmers.size(); r++){
      for(size_t i = 0; i < newKmers[r].size(); i++){
        T gap = newKmers[r][i] - prev;
        prev = newKmers[r][i];
        while(gap >= 0x80){
          keys.push_back((unsigned char)(gap | 0x80));
          gap >>= 7;
        }
        keys.push_back((unsigned char)gap);
        uint64_t bits = serializeDouble(newFreqcnt_lg[r][i]);
        typename unordered_map<uint64_t, uint32_t>::iterator it = where.find(bits);
        if(it == where.end()){
          it = where.insert(make_pair(bits, (uint32_t)dict.size())).first;
          dict.push_back(newFreqcnt_lg[r][i]);
        }
        idx.push_back(it->second);
      }
    }
  }
  uint8_t width = dict.size() <= 256 ? 1 : (dict.size() <= 65536 ? 2 : 4);
  vector<unsigned char> idx_bytes(idx.size() * width);
  for(size_t i = 0; i < idx.size(); i++){
    memcpy(&idx_bytes[i * width], &idx[i], width); // little-endian
  }
  uint64_t ndict = dict.size(), key_bytes = keys.size();

  path tmp_path = path(destination_path.native() + ".tmp");
  {
    IoSlot slot;
    std::ofstream out(tmp_path.native(), ios::binary);
    out<<"NBSAVE "<<SAVEFILE_VERSION<<" "<<kmerSize<<" "<<CANONICAL_KMERS<<"\n";
    out.write((const char*)&ngenomes_v, sizeof(ngenomes_v));
    out.write((const char*)&sumfreq_v, sizeof(sumfreq_v));
    out.write((const char*)&n, sizeof(n));
    out.write((const char*)&ndict, sizeof(ndict));
    out.write((const char*)dict.data(), ndict * sizeof(double));
    out.write((const char*)&width, sizeof(width));
    out.write((const char*)&key_bytes, sizeof(key_bytes));
    out.write((const char*)keys.data(), key_bytes);
    out.write((const char*)idx_bytes.data(), idx_bytes.size());
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
