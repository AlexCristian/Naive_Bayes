//
//  Genome.cpp
//  NB_C++
//
//  Created by Alexandru Cristian on 01/05/2017.
//
//

#include "Genome.hpp"
#include <stdexcept>
#include <zlib.h>
#include <algorithm>
#include <cctype>
#include <iterator>

bool Genome::STORE_ALL_NUMERATORS = false;
unsigned int Genome::KMER_SIZE = 0;
bool Genome::USE_PRIOR = true;

bool Genome::ScoreLess::operator()(const score& a, const score& b) const {
  if(a.first != b.first){
    return a.first < b.first;
  }
  // Equal scores: the lexicographically smaller class ID ranks higher.
  return a.second->getId() > b.second->getId();
}

Genome::Genome(path _kmr_path, path _sequence_path){
  kmr_path = _kmr_path;
  sequence_path = _sequence_path;
}

Genome::~Genome(){
  unload();
}

path Genome::getKmrPath(){
  return kmr_path;
}

void Genome::loadSequence(){
  // Remove this code, it never gets used & was never meant for production anyway
  std::ifstream in(sequence_path.native());

  sequence = new string;

  // Read the first line first (genome ID, other info)
  in>>getSequence();

  // Discard and overwrite with actual sequence.
  in>>getSequence();

  in.close();
  sequenceLoaded = true;
}

template <class F>
void Genome::parseKmrFile(F add){
  // Reads the whole file and parses "<k-mer> <count>" pairs by hand, about
  // 1.5x faster than stream extraction.
  // One read() into a pre-sized buffer: filling it through
  // istreambuf_iterator copies a character at a time and costs ~6x more.
  std::ifstream in(kmr_path.native(), ios::binary | ios::ate);
  if(!in.good()){
    cerr<<"[ERROR] Cannot read k-mer file "<<kmr_path.native()<<"\n";
    exit(1);
  }
  string buf((size_t)in.tellg(), '\0');
  in.seekg(0);
  in.read(&buf[0], buf.size());
  if(!in){
    cerr<<"[ERROR] Cannot read k-mer file "<<kmr_path.native()<<"\n";
    exit(1);
  }
  // Same characters as isspace() in the C locale, without the locale lookup.
  auto space = [](char c){ return c == ' ' || (c >= '\t' && c <= '\r'); };
  const char *p = buf.data(), *end = p + buf.size();
  while(true){
    while(p < end && space(*p)) p++;
    if(p == end) break;
    const char *kmer_start = p;
    while(p < end && !space(*p)) p++;
    size_t len = p - kmer_start;
    while(p < end && space(*p)) p++;
    if(p == end || !isdigit((unsigned char)*p)){
      cerr<<"[ERROR] Malformed k-mer count in "<<kmr_path.native()<<"\n";
      exit(1);
    }
    long long count = 0;
    while(p < end && isdigit((unsigned char)*p)){
      count = count * 10 + (*p - '0');
      p++;
    }
    if(len > 32 || (KMER_SIZE != 0 && len != KMER_SIZE)){
      cerr<<"[ERROR] k-mer \""<<string(kmer_start, len)<<"\" in "
          <<kmr_path.native()<<" does not match the k-mer size ("<<KMER_SIZE
          <<", max 32).\n";
      exit(1);
    }
    kmer_t kmer = 0;
    bool valid = true;
    for(size_t i = 0; valid && i < len; i++){
      kmer<<=2;
      switch(kmer_start[i]){
        case 'A': case 'a': break;
        case 'C': case 'c': kmer+=1; break;
        case 'G': case 'g': kmer+=2; break;
        case 'T': case 't': kmer+=3; break;
        // Ambiguous bases (N, IUPAC codes) carry no k-mer identity; skip.
        default: valid = false;
      }
    }
    if(valid){
      add(kmer, (int)count);
    }
  }
}

bool Genome::isFasta(){
  string name = kmr_path.filename().native();
  if(name.size() > 3 && name.compare(name.size() - 3, 3, ".gz") == 0){
    name = name.substr(0, name.size() - 3);
  }
  string ext = path(name).extension().native();
  return ext == ".fa" || ext == ".fasta" || ext == ".fna";
}

uint64_t Genome::estimatedKmers(){
  boost::system::error_code ec;
  uintmax_t size = file_size(kmr_path, ec);
  if(ec){
    return 0;
  }
  if(isFasta()){
    // About one k-mer per base; gzip packs DNA roughly 3.5x.
    string name = kmr_path.filename().native();
    bool gz = name.size() > 3 && name.compare(name.size() - 3, 3, ".gz") == 0;
    return gz ? size * 7 / 2 : size;
  }
  return size / (KMER_SIZE + 4); // one "<k-mer> <count>" line each
}

void Genome::countFastaKmers(vector<pair<kmer_t, int> >& out){
  out.clear();
  if(KMER_SIZE == 0){
    cerr<<"[ERROR] k-mer size must be set to count k-mers from FASTA.\n";
    exit(1);
  }
  // gzread() reads both gzipped and plain files.
  gzFile in = gzopen(kmr_path.native().c_str(), "rb");
  if(in == NULL){
    cerr<<"[ERROR] Cannot read FASTA file "<<kmr_path.native()<<"\n";
    exit(1);
  }
  const unsigned int k = KMER_SIZE;
  const kmer_t mask = (k == 32) ? ~(kmer_t)0 : (((kmer_t)1 << (2 * k)) - 1);
  const int rc_shift = 2 * (k - 1);
  vector<kmer_t> codes;
  kmer_t fwd = 0, rev = 0;
  unsigned int valid = 0;     // bases since the last break
  bool header = false, line_start = true;
  vector<char> buf(1 << 20);
  int got;
  while((got = gzread(in, buf.data(), buf.size())) > 0){
    for(int i = 0; i < got; i++){
      char c = buf[i];
      if(line_start && c == '>'){
        header = true;
        valid = 0;            // k-mers never span two sequences
      }
      line_start = (c == '\n');
      if(header){
        if(c == '\n'){
          header = false;
        }
        continue;
      }
      kmer_t b;
      switch(c){
        case 'A': case 'a': b = 0; break;
        case 'C': case 'c': b = 1; break;
        case 'G': case 'g': b = 2; break;
        case 'T': case 't': b = 3; break;
        case '\n': case '\r': case ' ': case '\t': continue;
        default: valid = 0; continue; // N and other codes break k-mers
      }
      fwd = ((fwd << 2) | b) & mask;
      rev = (rev >> 2) | ((3 - b) << rc_shift);
      if(++valid >= k){
        // With A<C<G<T coded 0-3, numeric order is lexicographic order, so
        // the smaller code is jellyfish's canonical k-mer.
        codes.push_back(fwd < rev ? fwd : rev);
      }
    }
  }
  int err;
  const char* msg = gzerror(in, &err);
  if(got < 0 || (err != Z_OK && err != Z_STREAM_END)){
    cerr<<"[ERROR] Cannot read FASTA file "<<kmr_path.native()<<": "<<msg<<"\n";
    exit(1);
  }
  gzclose(in);

  sort(codes.begin(), codes.end());
  for(size_t i = 0; i < codes.size(); ){
    size_t j = i;
    while(j < codes.size() && codes[j] == codes[i]){
      j++;
    }
    out.push_back(make_pair(codes[i], (int)(j - i)));
    i = j;
  }
}

void Genome::loadKmerCounts(){
  /*
   * Several classifier threads can ask for the same read at once. Re-check the
   * flag under the lock so only the first one loads, and the others wait for
   * it to finish instead of reloading into a fresh map.
   */
  lock_guard<mutex> guard(loadKmersLock);
  if(kmersLoaded){
    return;
  }

  kmer_counts = new unordered_map<kmer_t, int>;
  if(isFasta()){
    vector<pair<kmer_t, int> > counts;
    countFastaKmers(counts);
    for(size_t i = 0; i < counts.size(); i++){
      (*kmer_counts)[counts[i].first] += counts[i].second;
    }
  }else{
    parseKmrFile([this](kmer_t kmer, int count){ (*kmer_counts)[kmer] += count; });
  }
  kmersLoaded = true;
}

void Genome::loadSortedKmerCounts(vector<pair<kmer_t, int> >& out){
  if(isFasta()){
    countFastaKmers(out);
    return;
  }
  out.clear();
  parseKmrFile([&out](kmer_t kmer, int count){ out.push_back(make_pair(kmer, count)); });
  sort(out.begin(), out.end());
  // Combine repeats (jellyfish output has none, but hand-made files might).
  size_t w = 0;
  for(size_t r = 0; r < out.size(); r++){
    if(w > 0 && out[w-1].first == out[r].first){
      out[w-1].second += out[r].second;
    }else{
      out[w++] = out[r];
    }
  }
  out.resize(w);
}

void Genome::indexKmers(const vector<kmer_t>& batch_kmers){
  unordered_map<kmer_t, int>& counts = getKmerCounts();
  batch_index.clear();
  batch_count.clear();
  batch_index.reserve(counts.size());
  batch_count.reserve(counts.size());
  for(unordered_map<kmer_t, int>::iterator freq = counts.begin();
    freq != counts.end(); freq++){
      batch_index.push_back(lower_bound(batch_kmers.begin(), batch_kmers.end(),
                                        freq->first) - batch_kmers.begin());
      batch_count.push_back(freq->second);
    }
  unload();
}

void Genome::unload(){
  if(kmersLoaded){
    delete kmer_counts;
    kmersLoaded = false;
  }

  if(sequenceLoaded){
    delete sequence;
    sequenceLoaded = false;
  }
}

unordered_map<kmer_t, int>& Genome::getKmerCounts(){
  if(!kmersLoaded){
    loadKmerCounts();
  }
  return *kmer_counts;
}

string& Genome::getSequence(){
  return *sequence;
}

void Genome::scoreClass(Class<kmer_t>* cl, const vector<double>& batch_lg){
  accumulator_set<double, features<tag::sum_kahan> > current, sumfrq;

  if(USE_PRIOR){
    current(cl->getNGenomes_lg());
  }
  ostringstream strs;
  if(NB::debug_flag == NB::Debug::LOG_ALL){
    strs<<"("<<cl->getId()<<"): "<<sum_kahan(current);
  }

  // Same terms in the same order as the original per-k-mer lookups, so the
  // Kahan sums (and predictions) are bit-identical.
  for(size_t i = 0; i < batch_index.size(); i++){
      double lg = batch_lg[batch_index[i]];
      sumfrq(batch_count[i]);
      current(batch_count[i] * lg);
      if(NB::debug_flag == NB::Debug::LOG_ALL){
        strs<<" + "<<batch_count[i]<<" * "<<lg;
      }
  }
  current(-sum_kahan(sumfrq) * cl->getSumFreq_lg());
  if(NB::debug_flag == NB::Debug::LOG_ALL){
    strs<<" - "<<sum_kahan(sumfrq)<<" * "<<cl->getSumFreq_lg()<<" = ";
    strs<<sum_kahan(current)<<"\n";
  }

  numeratorAccess.lock();

  if(!Genome::STORE_ALL_NUMERATORS){
    double candidateNumerator = sum_kahan(current);
    if(maximumNumeratorClass == NULL
      || ScoreLess()(make_pair(maximumNumerator, maximumNumeratorClass),
                     make_pair(candidateNumerator, cl))){
        maximumNumerator = candidateNumerator;
        maximumNumeratorClass = cl;
      }
  }else{
    numerator.push(make_pair(sum_kahan(current), cl));
    if(NB::debug_flag == NB::Debug::LOG_ALL){
      cout<<strs.str();
    }
  }

  numeratorAccess.unlock();
}

Genome::score Genome::getMaximum(){
  numeratorAccess.lock();
  Genome::score sc = make_pair(maximumNumerator, maximumNumeratorClass);
  numeratorAccess.unlock();

  return sc;
}

Genome::pqueue Genome::getConfidences(){
  accumulator_set<double, features<tag::sum_kahan> >  p_denominator(1);
  double max, denominator;

  pqueue num_cpy = numerator;

  max = num_cpy.top().first;
  vector<score> cache; cache.push_back(num_cpy.top());
  num_cpy.pop();

  while(!num_cpy.empty()){
    p_denominator(exp(num_cpy.top().first - max));
    cache.push_back(num_cpy.top());
    num_cpy.pop();
  }
  denominator = max + log(sum_kahan(p_denominator));

  pqueue confidence;
  for(vector<score>::iterator num=cache.begin();
    num != cache.end(); num++){
    double confidence_lg = num->first - denominator;
    confidence.push(make_pair(confidence_lg, num->second));
    //processFreqs(getKmerCounts(), num->second, confidence_lg);
  }

  return confidence;
}

int Genome::size(){
  return getSequence().size();
}

long long int Genome::computeAlignmentScore(Genome* seq){
  // The Smith-Waterman scoring (distance & gap penalty) was never written.
  // Fail loudly instead of running a loop over uninitialized state.
  throw logic_error("Genome::computeAlignmentScore is not implemented");
}

string Genome::charAt(int pos){
  if(!sequenceLoaded){
    loadSequence();
  }
  return string(1, getSequence()[pos]);
}

long long int Genome::max(long long int a,
                  long long int b,
                  long long int c,
                  long long int d){
  if(a>=b && a>=c && a>=d){
    return a;
  }else if(b>=a && b>=c && b>=d){
    return b;
  }else if(c>=b && c>=a && c>=d){
    return c;
  }else{
    return d;
  }
}
