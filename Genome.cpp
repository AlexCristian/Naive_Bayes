//
//  Genome.cpp
//  NB_C++
//
//  Created by Alexandru Cristian on 01/05/2017.
//
//

#include "Genome.hpp"
#include <stdexcept>

bool Genome::STORE_ALL_NUMERATORS = false;
unsigned int Genome::KMER_SIZE = 0;

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

  std::ifstream in(kmr_path.native());

  kmer_counts = new unordered_map<kmer_t, int>;

  kmer_t kmer; int kmer_count; string kmer_s;
  while(in>>kmer_s>>kmer_count){
    if(kmer_s.size() > 32 || (KMER_SIZE != 0 && kmer_s.size() != KMER_SIZE)){
      cerr<<"[ERROR] k-mer \""<<kmer_s<<"\" in "<<kmr_path.native()
          <<" does not match the k-mer size ("<<KMER_SIZE<<", max 32).\n";
      exit(1);
    }
    kmer = 0;
    bool valid = true;
    for(string::iterator it = kmer_s.begin(); valid && it != kmer_s.end(); it++){
      kmer<<=2;
      switch(*it){
        case 'A':
        case 'a':
          kmer+=0;
          break;
        case 'C':
        case 'c':
          kmer+=1;
          break;
        case 'G':
        case 'g':
          kmer+=2;
          break;
        case 'T':
        case 't':
          kmer+=3;
          break;
        default:
          // Ambiguous bases (N, IUPAC codes) carry no k-mer identity; skip.
          valid = false;
      }
    }
    if(valid){
      (*kmer_counts)[kmer]+=kmer_count;
    }
  }
  in.close();
  kmersLoaded = true;
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

void Genome::computeClassificationNumerator(Class<kmer_t>* cl){
  accumulator_set<double, features<tag::sum_kahan> > current, sumfrq;
  
  current(cl->getNGenomes_lg());
  ostringstream strs;
  if(NB::debug_flag == NB::Debug::LOG_ALL){
    strs<<"("<<cl->getId()<<"): "<<sum_kahan(current);
  }

  unordered_map<kmer_t, int>& counts = getKmerCounts();
  for(unordered_map<kmer_t, int>::iterator freq = counts.begin();
    freq != counts.end(); freq++){
      sumfrq(freq->second);
      current(freq->second * cl->getFreqCount_lg(freq->first));
      if(NB::debug_flag == NB::Debug::LOG_ALL){
        strs<<" + "<<freq->second<<" * "<<cl->getFreqCount_lg(freq->first);
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
