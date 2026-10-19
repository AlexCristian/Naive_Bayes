//
//  TaskPool.hpp
//  NB_C++
//
//  A fixed set of worker threads that run small tasks, and open new jobs
//  (classes to train or classify) whenever a worker would otherwise be idle
//  and the memory budget allows.
//

#ifndef TaskPool_hpp
#define TaskPool_hpp

#include <condition_variable>
#include <cstdint>
#include <deque>
#include <functional>
#include <mutex>
#include <thread>
#include <vector>
using namespace std;

class TaskPool {
public:
  /**
   * Tries to start the next job. Called by an idle worker, with the pool's
   * lock held: it must not call submit() or release().
   * @param bytes Set to the memory the job will hold until it calls release().
   * @return      The job's first task, or an empty function when no jobs are
   *              left.
   */
  typedef function<function<void()>(uint64_t& bytes)> Opener;

  /**
   * @param nthreads   Number of worker threads.
   * @param budget     Memory (bytes) that open jobs may hold at once. A job
   *                   larger than the budget still runs, but alone.
   */
  TaskPool(int nthreads, uint64_t budget);

  /**
   * Runs until every job has been opened and every task has finished.
   * @param peek   Returns the memory the next job needs, or 0 if none is left.
   * @param open   Starts the next job (see Opener).
   */
  void run(function<uint64_t()> peek, Opener open);

  /** Queues a task; callable from inside tasks. */
  void submit(function<void()> task);

  /** Returns memory a job no longer holds; callable from inside tasks. */
  void release(uint64_t bytes);

private:
  int nthreads;
  uint64_t budget, used = 0;
  mutex m;
  condition_variable cv;
  deque<function<void()> > tasks;
  int busy = 0;
  bool exhausted = false;
  function<uint64_t()> peek;
  Opener open;

  void worker();
};

inline TaskPool::TaskPool(int _nthreads, uint64_t _budget)
  : nthreads(_nthreads < 1 ? 1 : _nthreads), budget(_budget) {}

inline void TaskPool::submit(function<void()> task){
  {
    lock_guard<mutex> lk(m);
    tasks.push_back(move(task));
  }
  cv.notify_one();
}

inline void TaskPool::release(uint64_t bytes){
  {
    lock_guard<mutex> lk(m);
    used -= bytes;
  }
  cv.notify_all();
}

inline void TaskPool::run(function<uint64_t()> _peek, Opener _open){
  peek = _peek;
  open = _open;
  exhausted = false;
  used = 0;
  vector<thread> workers;
  for(int i = 0; i < nthreads; i++){
    workers.emplace_back(&TaskPool::worker, this);
  }
  for(size_t i = 0; i < workers.size(); i++){
    workers[i].join();
  }
}

inline void TaskPool::worker(){
  unique_lock<mutex> lk(m);
  while(true){
    // Tasks of jobs already open come first, so they finish and free memory.
    if(!tasks.empty()){
      function<void()> task = move(tasks.front());
      tasks.pop_front();
      busy++;
      lk.unlock();
      task();
      lk.lock();
      busy--;
      cv.notify_all();
      continue;
    }
    // Otherwise open another job, if one is left and fits the budget. A job
    // that doesn't fit waits until others release memory (or runs alone).
    if(!exhausted){
      uint64_t need = peek();
      if(need == 0){
        exhausted = true;
        cv.notify_all();
      }else if(used == 0 || used + need <= budget){
        uint64_t bytes = 0;
        function<void()> first = open(bytes);
        if(!first){
          exhausted = true;
          cv.notify_all();
        }else{
          used += bytes;
          busy++;
          lk.unlock();
          first();
          lk.lock();
          busy--;
          cv.notify_all();
        }
        continue;
      }
    }
    if(exhausted && tasks.empty() && busy == 0){
      cv.notify_all();
      return;
    }
    cv.wait(lk);
  }
}

#endif /* TaskPool_hpp */
