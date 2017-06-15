//
//  NB.cpp
//  NB_C++
//
//  Created by Alexandru Cristian on 06/04/2017.
//
//

#include "NB.hpp"

NB::Debug NB::debug_flag = NB::Debug::NO_LOG;

NB::NB(int _kmer_size, path _save_dir, int _nthreads)  {
  nthreads = _nthreads;
  debug_flag = Debug::NO_LOG;
  kmer_size = _kmer_size;
  save_dir = _save_dir;

  load();
}

NB::~NB() {
  for(unordered_map<string, Class<int>* >::iterator iter = classes.begin();
    iter != classes.end(); iter++)
      delete iter->second;
}

void NB::load(){
  vector<path> res = Diskutil::getItemsInDir(save_dir);
  for(vector<path>::iterator iter = res.begin(); iter != res.end(); iter++){

    if(!Diskutil::hasFileExtension(*iter, Diskutil::SAVE_FILE_EXT)){
      continue;
    }

    string filename = iter->filename().native();
    string cls_s = filename.substr(0,filename.rfind('-'));
    addClass(new Class<int>(cls_s, kmer_size, *iter));
  }
}

void NB::save(){
  for(unordered_map<string, Class<int>* >::iterator iter = classes.begin();
      iter != classes.end(); iter++){
        if(iter->second->loaded()){
          iter->second->save();
        }
      }
}

void NB::addClass(Class<int>* cl){
    classes[cl->getId()] = cl;
}

Class<int>* NB::getClass(string cl_id){
  if(classes.find(cl_id) == classes.end()){
    return NULL;
  }
  return classes[cl_id];
}

int NB::getKmerSize(){
  return kmer_size;
}

path NB::getSavedir(){
  return save_dir;
}

void NB::addClassToUpdateQueue(Class<int>* cl){
  classesToProcess.push(cl);
}

void NB::processClassUpdates(){
  progress = 1;
  total = classesToProcess.size();

  if(NB::debug_flag == NB::Debug::LOG_SOME){
    cout<<"Started training on batch...\n";
    cout.flush();
  }

  int c_nthreads = nthreads;
  if(nthreads > classes.size()){
    c_nthreads = classes.size();
  }

  vector<thread*> workers(c_nthreads);
  for(int i=0; i < c_nthreads; i++){
    workers[i] = new thread(&NB::trainThreadController, this);
  }
  for(int i=0; i < c_nthreads; i++){
    workers[i]->join();
    delete workers[i];
  }
}

void NB::trainThreadController(){
  bool exitCondition = false;
  while(!exitCondition){
    classQueueAccess.lock();
    if(classesToProcess.empty()){
      exitCondition = true;
      classQueueAccess.unlock();
      continue;
    }

    Class<int>* cl = classesToProcess.front();
    classesToProcess.pop();

    if(NB::debug_flag == NB::Debug::LOG_SOME
       || NB::debug_flag == NB::Debug::LOG_ALL){
      cout<<"("<<progress++<<"/"<<total<<") Training genomes for class ";
      cout<<cl->getId()<<"\n";
      cout.flush();
    }
    classQueueAccess.unlock();

    cl->load();
    cl->addGenomesInQueue();
    cl->save();
    cl->unload();
  }
}

void NB::classifyThreadController(vector<Genome*>& reads){
  bool exitCondition = false;
  while(!exitCondition){
    classQueueAccess.lock();
    if(classesToProcess.empty()){
      exitCondition = true;
      classQueueAccess.unlock();
      continue;
    }

    Class<int>* cl = classesToProcess.front();
    classesToProcess.pop();

    if(NB::debug_flag == NB::Debug::LOG_SOME){
      cout<<"("<<progress++<<"/"<<total<<") Computing confidence for class ";
      cout<<cl->getId()<<"\n";
      cout.flush();
    }
    classQueueAccess.unlock();

    cl->load();
    cl->computeBatchNumerators(reads);
    cl->unload();
  }
}

void NB::classify(vector<Genome*> reads){
  progress = 1;
  total = classes.size();

  if(NB::debug_flag == NB::Debug::LOG_SOME
     || NB::debug_flag == NB::Debug::LOG_ALL){
    cout<<"Started classifying...\n";
    cout.flush();
  }

  for(unordered_map<string, Class<int>* >::iterator iter = classes.begin();
    iter != classes.end(); iter++){
      Class<int> *cl = iter->second;
      classesToProcess.push(cl);
    }

  int c_nthreads = nthreads;
  if(nthreads > classes.size()){
    c_nthreads = classes.size();
  }
  vector<thread*> workers(c_nthreads);
  for(int i=0; i < c_nthreads; i++){
    workers[i] = new thread(&NB::classifyThreadController, this, ref(reads));
  }
  for(int i=0; i < c_nthreads; i++){
    workers[i]->join();
    delete workers[i];
  }

  for(vector<Genome*>::iterator iter = reads.begin();
    iter != reads.end(); iter++){
      (*iter)->unload();
    }
}

int NB::min(int a, int b, int c){
  if(a<=b && a<=c){
    return a;
  }else if(b<=a && b<=c){
    return b;
  }else{
    return c;
  }
}

int NB::levenshteinDistance(string src, string tgt){
  vector<vector <int> > dist;

  int s_len = src.length(), t_len = tgt.length(), elem_cost;

  // Initialize
  dist.emplace_back(s_len+1);
  dist.emplace_back(s_len+1);
  for(int j=0; j <= s_len; j++){
    dist[0][j] = j;
  }

  for(int i=1; i <= t_len; i++){
    int adj_i = 1;
    dist[adj_i][0] = i;
    for(int j=1; j <= s_len; j++){
      if(src[j] == tgt[i])
        elem_cost = 0;
      else
        elem_cost = 1;

      dist[adj_i][j] = min(dist[adj_i-1][j] + 1,
        dist[adj_i][j-1] + 1,
        dist[adj_i-1][j-1] + elem_cost);
    }

    dist[0].swap(dist[1]);
  }

  return dist[0][s_len];
}
