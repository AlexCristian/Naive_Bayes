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
  for(unordered_map<string, Class<kmer_t>* >::iterator iter = classes.begin();
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
    Class<kmer_t>* cl = new Class<kmer_t>(cls_s, kmer_size, *iter);
    // Catch a savedir trained with another k or counting mode up front,
    // instead of silently scoring against mismatched k-mers.
    cl->checkSavefileHeader(*iter);
    addClass(cl);
  }
}

void NB::save(){
  for(unordered_map<string, Class<kmer_t>* >::iterator iter = classes.begin();
      iter != classes.end(); iter++){
        if(iter->second->loaded()){
          iter->second->save();
        }
      }
}

void NB::addClass(Class<kmer_t>* cl){
    classes[cl->getId()] = cl;
}

Class<kmer_t>* NB::getClass(string cl_id){
  if(classes.find(cl_id) == classes.end()){
    return NULL;
  }
  return classes[cl_id];
}

size_t NB::numClasses(){
  return classes.size();
}

int NB::getKmerSize(){
  return kmer_size;
}

path NB::getSavedir(){
  return save_dir;
}

bool NB::isClassQueuedForUpdate(Class<kmer_t>* cl){
  return classesQueuedForUpdate.count(cl) != 0;
}

bool NB::addClassToUpdateQueue(Class<kmer_t>* cl){
  // Two queue entries for one class would let two threads train it at once.
  if(!classesQueuedForUpdate.insert(cl).second){
    return false;
  }
  classesToProcess.push(cl);
  return true;
}

void NB::queueAllClassesForUpdate(){
  for(unordered_map<string, Class<kmer_t>* >::iterator iter = classes.begin();
    iter != classes.end(); iter++){
      addClassToUpdateQueue(iter->second);
    }
}

void NB::setWorkMemory(uint64_t bytes){
  work_memory = bytes;
}

int NB::splitCount(size_t items, size_t per_task){
  // Enough tasks to keep every thread busy (a few per thread, to even out
  // uneven pieces), but none smaller than per_task items.
  size_t n = items / per_task;
  size_t most = 4 * (size_t)std::max(nthreads, 1);
  return (int)std::max((size_t)1, std::min(n, most));
}

/*
 * Training runs on a task pool. Each class in the batch is a job:
 *   load the class -> parse each new genome (one task per genome) ->
 *   commit each k-mer range (one task per range) -> save and unload.
 * Idle threads open the next job, largest first, while the memory budget
 * allows, so big classes are split across all threads and small ones fill
 * the gaps. The saved values don't depend on the order tasks run in.
 */
void NB::processClassUpdates(){
  progress = 1;
  total = classesToProcess.size();

  if(NB::debug_flag == NB::Debug::LOG_SOME){
    cout<<"Started training on batch...\n";
    cout.flush();
  }

  struct Job {
    Class<kmer_t>* cl;
    uint64_t work, bytes;
  };
  vector<Job> jobs;
  while(!classesToProcess.empty()){
    Class<kmer_t>* cl = classesToProcess.front();
    classesToProcess.pop();
    uint64_t saved = 0, new_kmers = cl->queuedKmerEstimate();
    boost::system::error_code ec;
    uintmax_t size = file_size(cl->getSavefilePath(), ec);
    if(!ec){
      saved = size;
    }
    // In memory: the class twice during the commit (old and new arrays),
    // plus the new counts, at ~16 bytes per k-mer.
    jobs.push_back({cl, saved + 16 * new_kmers, 2 * saved + 16 * new_kmers + 1});
  }
  classesQueuedForUpdate.clear();
  sort(jobs.begin(), jobs.end(),
       [](const Job& a, const Job& b){ return a.work > b.work; });

  TaskPool pool(nthreads, work_memory);
  size_t next = 0;
  pool.run(
    [&](){ return next < jobs.size() ? jobs[next].bytes : 0; },
    [&](uint64_t& bytes) -> function<void()> {
      if(next >= jobs.size()){
        return function<void()>();
      }
      Job job = jobs[next++];
      bytes = job.bytes;
      if(NB::debug_flag == NB::Debug::LOG_SOME
         || NB::debug_flag == NB::Debug::LOG_ALL){
        cout<<"("<<progress++<<"/"<<total<<") Training genomes for class ";
        cout<<job.cl->getId()<<"\n";
        cout.flush();
      }
      return [this, job, &pool](){ trainClass(job.cl, job.bytes, pool); };
    });
}

void NB::trainClass(Class<kmer_t>* cl, uint64_t bytes, TaskPool& pool){
  cl->load();
  vector<Genome*> todo = cl->takeGenomesToTrain();

  size_t estimate = cl->numKmers();
  for(size_t i = 0; i < todo.size(); i++){
    estimate += todo[i]->estimatedKmers();
  }
  size_t R = splitCount(estimate, 1 << 20);

  // A new class has no k-mers to split on yet: parse one genome first and
  // split on its k-mers instead.
  shared_ptr<vector<Genome*> > rest = make_shared<vector<Genome*> >(todo);
  if(cl->numKmers() < R && !rest->empty()){
    vector<pair<kmer_t, int> > sample;
    rest->front()->loadSortedKmerCounts(sample);
    cl->beginUpdate(R, &sample);
    cl->addGenomeCounts(sample);
    rest->erase(rest->begin());
  }else{
    cl->beginUpdate(R, NULL);
  }

  auto commit = [cl, bytes, &pool](){
    size_t ranges = cl->numRanges();
    shared_ptr<atomic<size_t> > left = make_shared<atomic<size_t> >(ranges);
    for(size_t r = 0; r < ranges; r++){
      pool.submit([cl, r, left, bytes, &pool](){
        cl->commitRange(r);
        if(--(*left) == 0){
          cl->save();
          cl->unload();
          pool.release(bytes);
        }
      });
    }
  };

  if(rest->empty()){
    commit();
    return;
  }
  shared_ptr<atomic<size_t> > left = make_shared<atomic<size_t> >(rest->size());
  for(size_t i = 0; i < rest->size(); i++){
    Genome* g = (*rest)[i];
    pool.submit([cl, g, left, commit](){
      vector<pair<kmer_t, int> > counts;
      g->loadSortedKmerCounts(counts);
      cl->addGenomeCounts(counts);
      if(--(*left) == 0){
        commit();
      }
    });
  }
}

void NB::forEachRead(vector<Genome*>& reads, function<void(Genome*)> fn){
  int c_nthreads = std::max(1, std::min((int)reads.size(), nthreads));
  vector<thread> workers;
  for(int t = 0; t < c_nthreads; t++){
    workers.emplace_back([&reads, &fn, t, c_nthreads](){
      for(size_t i = t; i < reads.size(); i += c_nthreads){
        fn(reads[i]);
      }
    });
  }
  for(size_t t = 0; t < workers.size(); t++){
    workers[t].join();
  }
}

/*
 * Classification runs on the same kind of task pool. Each class is a job:
 *   load the class -> look up its log counts for slices of the batch's
 *   k-mer list (one task per slice) -> score chunks of reads (one task per
 *   chunk) -> unload.
 * Every thread can work on the same loaded class, because nothing in it
 * changes while classifying. Each read is still summed in the same order,
 * so scores are the same as before.
 */
void NB::classify(vector<Genome*> reads){
  progress = 1;
  total = classes.size();

  if(NB::debug_flag == NB::Debug::LOG_SOME
     || NB::debug_flag == NB::Debug::LOG_ALL){
    cout<<"Started classifying...\n";
    cout.flush();
  }

  // Load every read once, then list the batch's distinct k-mers (sorted) and
  // point each read's k-mers into that list.
  forEachRead(reads, [](Genome* g){ g->getKmerCounts(); });
  shared_ptr<vector<kmer_t> > batch = make_shared<vector<kmer_t> >();
  vector<kmer_t>& batch_kmers = *batch;
  for(vector<Genome*>::iterator iter = reads.begin();
    iter != reads.end(); iter++){
      unordered_map<kmer_t, int>& counts = (*iter)->getKmerCounts();
      for(unordered_map<kmer_t, int>::iterator k = counts.begin();
        k != counts.end(); k++){
          batch_kmers.push_back(k->first);
        }
    }
  sort(batch_kmers.begin(), batch_kmers.end());
  batch_kmers.erase(unique(batch_kmers.begin(), batch_kmers.end()),
                    batch_kmers.end());
  forEachRead(reads, [&batch_kmers](Genome* g){ g->indexKmers(batch_kmers); });

  struct Job {
    Class<kmer_t>* cl;
    uint64_t bytes;
  };
  vector<Job> jobs;
  size_t B = batch_kmers.size();
  for(unordered_map<string, Class<kmer_t>* >::iterator iter = classes.begin();
    iter != classes.end(); iter++){
      uint64_t saved = 0;
      boost::system::error_code ec;
      uintmax_t size = file_size(iter->second->getSavefilePath(), ec);
      if(!ec){
        saved = size;
      }
      jobs.push_back({iter->second, saved + B * sizeof(double) + 1});
    }
  sort(jobs.begin(), jobs.end(),
       [](const Job& a, const Job& b){ return a.bytes > b.bytes; });

  int slices = splitCount(B, 1 << 16);
  int chunks = splitCount(reads.size(), 256);
  TaskPool pool(nthreads, work_memory);
  size_t next = 0;
  pool.run(
    [&](){ return next < jobs.size() ? jobs[next].bytes : 0; },
    [&](uint64_t& bytes) -> function<void()> {
      if(next >= jobs.size()){
        return function<void()>();
      }
      Job job = jobs[next++];
      bytes = job.bytes;
      if(NB::debug_flag == NB::Debug::LOG_SOME){
        cout<<"("<<progress++<<"/"<<total<<") Computing confidence for class ";
        cout<<job.cl->getId()<<"\n";
        cout.flush();
      }
      return [job, batch, &reads, slices, chunks, &pool](){
        Class<kmer_t>* cl = job.cl;
        cl->load();
        size_t B = batch->size();
        shared_ptr<vector<double> > lg = make_shared<vector<double> >(B);
        shared_ptr<atomic<int> > lookups = make_shared<atomic<int> >(slices);
        for(int s = 0; s < slices; s++){
          size_t qlo = B * s / slices, qhi = B * (s + 1) / slices;
          pool.submit([=, &reads, &pool](){
            cl->lookupFreqCounts_lg(*batch, qlo, qhi, *lg);
            if(--(*lookups) > 0){
              return;
            }
            shared_ptr<atomic<int> > scored = make_shared<atomic<int> >(chunks);
            for(int c = 0; c < chunks; c++){
              size_t rlo = reads.size() * c / chunks;
              size_t rhi = reads.size() * (c + 1) / chunks;
              pool.submit([=, &reads, &pool](){
                for(size_t i = rlo; i < rhi; i++){
                  reads[i]->scoreClass(cl, *lg);
                }
                if(--(*scored) == 0){
                  cl->unload();
                  pool.release(job.bytes);
                }
              });
            }
          });
        }
      };
    });

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
      if(src[j-1] == tgt[i-1])
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
