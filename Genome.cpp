//
//  Genome.cpp
//  NB_C++
//
//  Created by Alexandru Cristian on 01/05/2017.
//
//

#include "Genome.hpp"

bool Genome::STORE_ALL_NUMERATORS = false;

Genome::Genome(path _kmr_path, path _sequence_path){
  kmr_path = _kmr_path;
  sequence_path = _sequence_path;
}

Genome::~Genome(){
  unload();
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
  if(!loadKmersLock.try_lock()){
    /*
     * Another thread is trying to load this at the same time,
     * wait for it to finish and then return.
     */
    loadKmersLock.lock();
    loadKmersLock.unlock();

    // Since we're here, another thread has already loaded this data. Return.
    return;
  }
  std::ifstream in(kmr_path.native());

  kmer_counts = new unordered_map<int, int>;

  int kmer, kmer_count; string kmer_s;
  while(in>>kmer_s>>kmer_count){
    kmer = 0;
    for(string::iterator it = kmer_s.begin(); it != kmer_s.end(); it++){
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
      }
    }
    (*kmer_counts)[kmer]=kmer_count;
  }
  in.close();
  kmersLoaded = true;
  loadKmersLock.unlock();
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

unordered_map<int, int>& Genome::getKmerCounts(){
  if(!kmersLoaded){
    loadKmerCounts();
  }
  return *kmer_counts;
}

string& Genome::getSequence(){
  return *sequence;
}

void Genome::computeClassificationNumerator(Class<int>* cl){
  accumulator_set<double, features<tag::sum_kahan> > current, sumfrq;
  
  current(cl->getNGenomes_lg());
  ostringstream strs;
  if(NB::debug_flag == NB::Debug::LOG_ALL){
    strs<<"("<<cl->getId()<<"): "<<sum_kahan(current);
  }

  for(unordered_map<int, int>::iterator freq = getKmerCounts().begin();
    freq != getKmerCounts().end(); freq++){
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
      || maximumNumerator < candidateNumerator){
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
  Genome *a, *b;
  if(seq->size() < getSequence().size()){
    a=seq;
    b=this;
  }else{
    a=this;
    b=seq;
  }

  vector<vector<long long int> > cost(2);
  int v_len = a->size();
  for(int i=0; i < v_len; i++){
    cost[0].push_back(0);
    cost[1].push_back(0);
  }

  vector<int> cmax(v_len), c_ind(v_len);

  int rep_n = b->size(), lmax, l_ind, highscore, hs_ind=-1;
  for(int i=0; i < rep_n; i++){
    for(int j=1; i < v_len; j++){
      throw exception();
      cost[1][j] = max(cost[0][j-1] /*+ distance(a->charAt(i), b->charAt(j))*/,
                       cmax[j] /*- penalty(j - c_ind[j])*/,
                       (l_ind != -1) ? lmax /*- penalty(j - l_ind)*/ : 0,
                       0);
      if(cost[1][j] > lmax || l_ind == -1){
        lmax = cost[1][j];
        l_ind = j;
      }
      if(cost[1][j] > cmax[j] /*- penalty(j - c_ind[j])*/){
        cmax[j] = cost[1][j];
        c_ind[j] = j;
      }
      if(cost[1][j] > highscore || hs_ind == -1){
        hs_ind = 1;
        highscore = cost[1][j];
      }
    }
    cost[0].swap(cost[1]);
  }

  return highscore;
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
