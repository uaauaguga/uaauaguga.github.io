---
layout: post
title:  "Parallel Programming in Practice"
date:   2021-06-15 13:21:52 +0800
usemathjax: true
categories: jekyll update
---

## Some basic concepts

- [process](https://en.wikipedia.org/wiki/Process_(computing))
- [thread](https://en.wikipedia.org/wiki/Thread_(computing)) 
- [Cooperative multitasking](https://en.wikipedia.org/wiki/Cooperative_multitasking)
- mutex (互斥锁),用于加锁
- conditional variable (条件变量)，用于等待
- semaphore (信号量)，用于限制资源使用



## Simultaneously excute multiple programs

### GNU [parallel](https://www.gnu.org/software/parallel/)

- <https://www.gnu.org/software/parallel/>

```bash
parallel -j 4 echo "hell await {1} {2}" ::: A B C ::: 1 2 3
#hell await A 1
#hell await A 2
#hell await A 3
#hell await B 1
#hell await B 2
#hell await B 3
#hell await C 1
#hell await C 2
#hell await C 3
```

- Here is some examples <https://www.gnu.org/software/parallel/parallel.html#examples>
- A example for running a bash loop in parallel

```bash
#!/bin/bash
for a in 3 7 11 15 19 23;do
for b in 0.0 0.2 0.4 0.6 0.8 1.0 ;do
  sem -j 4  "scripts.py -a $a -b $b > log/${a}-${b}.log"
done
done
sem --wait
# sem is an alias for parallel --semaphore
```

## GNU [make](https://www.gnu.org/software/make/manual/make.html)

```text
#cat sample_ids.txt
0001
0002
0003
0004
```

```makefile
# the != operator read variables from bash command
idx != cat sample_ids.txt 
# pattern substitution
tgt = $(idx:%=count/%.txt)

.PHONY: all

# inform make not to remove intermediate files
.SECONDARY: data/%.txt

all: $(tgt)

# count occurence of each character in the random string
count/%.txt: data/%.txt
	scripts/count-character.py -i $^ -o $@

# generate some random strings
data/%.txt:
	openssl rand -base64 100000 > $@

```

- `make --jobs 4` will spawn 4 process for parallel processing the pipeline

## [snakemake](https://snakemake.readthedocs.io/en/stable/)

- Inspired by gnu make, but more firendly
- Better support for cluster environment
- Native support for python scripting
- Support for multiple named wildcards



## Paralleling in programming language
  
### C++

- [pthread](https://man7.org/linux/man-pages/man7/pthreads.7.html): gnu support for multithreading in C
- [openmp](https://www.openmp.org/): higher level paralleling
- [std::thread](https://en.cppreference.com/w/cpp/thread/thread): multithreading feature added in C++ 11

- Discussion about pros and cons of these three methods
  - <https://stackoverflow.com/questions/3949901/pthreads-vs-openmp>
  - <https://stackoverflow.com/questions/13134186/c11-stdthreads-vs-posix-threads>
  - <https://www.zhihu.com/question/36236334>

- meme, bwa and STAR use pthread, kraken2 use openmp, seems few uses thread in C++ 11 ...

- Introductions for [pthread](https://man7.org/linux/man-pages/man7/pthreads.7.html)
  - <https://www.cs.cmu.edu/afs/cs/academic/class/15492-f07/www/pthreads.html>
  - <https://hpc-tutorials.llnl.gov/posix/>
  - <https://github.com/hailinzeng/Programming-POSIX-Threads>


- An example for `pthread` with mutex

```c++
#include<stdio.h>
#include<pthread.h>

#define NTHREADS 10
static pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
static int glob = 0;

void * func(void * args){
    printf("Thread id is: %ld\n", pthread_self());
    pthread_mutex_lock(&mutex);
    glob ++;
    pthread_mutex_unlock(&mutex);
}

int main(int argc, char* argv[]){
    pthread_t threads[NTHREADS];
    for(unsigned int i=0;i<NTHREADS;i++){
        pthread_create(&threads[i],NULL,func,NULL);
    }
    for(unsigned int j=0;j<NTHREADS;j++){
        pthread_join(threads[j],NULL);
    }  
    printf("The global variable is now: %d\n",glob);
    return 0;
}
```

- An example for `pthread` with conditional variable

```c++
#include<stdio.h>
#include<stdlib.h>
#include<pthread.h>

pthread_mutex_t count_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t condition_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_cond_t condition_cond = PTHREAD_COND_INITIALIZER;

static int glob = 0;

void * func1(void * args){
  while(1){
    printf("func 1 try to aquire lock\n");
    pthread_mutex_lock(&condition_mutex);
    printf("func 1 aquired lock\n");
    while(glob>=3&&glob<=7){
        printf("thread 1 start waiting ...\n");
        pthread_cond_wait(&condition_cond,&condition_mutex);
        printf("thread 1 stop waiting.\n");
    }
    pthread_mutex_lock(&count_mutex);
    glob++;
    printf("func 1 update to glbal variable to: %d\n",glob);
    pthread_mutex_unlock(&count_mutex);
    pthread_mutex_unlock(&condition_mutex);
    printf("func 1 release lock\n");
    if(glob>=10)return NULL;
  }
}

void * func2(void * args){
  while(1){
    printf("func 2 try to aquire lock\n");
    pthread_mutex_lock(&condition_mutex);
    printf("func 2 aquired lock\n");
    if(glob<3||glob>7){
        printf("thread 2 is sending signal...\n");
        pthread_cond_signal(&condition_cond);
    }
    
    pthread_mutex_lock(&count_mutex);
    glob++;
    printf("func 2 update to glbal variable to: %d\n",glob);
    pthread_mutex_unlock(&count_mutex);
    pthread_mutex_unlock(&condition_mutex);
    printf("func 2 release lock\n");
    if(glob>=10)return NULL;
  }
}

int main(int argc,char * argv[]){
  // [3,10] will be printed out by function 1
  pthread_t t1,t2;
  pthread_create(&t1,NULL,&func1,NULL);
  pthread_create(&t2,NULL,&func2,NULL);
  pthread_join(t1,NULL);
  pthread_join(t2,NULL);
  return 0;
}
```


### python

- Due to the global interpreter lock (GIL), pyhton is lack in native support for multithreading. 
- multithreading in python is not real multithreading, but it still helps in IO bounded problems.

- CPU bounded problem
  - multiprocessing
- IO bounded problem
  - multithreading
  - asyncio

- <https://realpython.com/python-concurrency/>

- OO style multi[threading](https://docs.python.org/3/library/threading.html)

```python
#!/usr/bin/env python
from threading import Thread
from time import sleep
from random import randint
class MThread(Thread):
    def __init__(self,thread_id,delay):
        Thread.__init__(self)
        self.thread_id = thread_id
        self.delay = delay
    def run(self):
        print(f"{self.thread_id} start")
        sleep(self.delay)
        print(f"{self.thread_id} stop")
def main():
    threads = []
    for i in range(10):
        t = MThread(thread_id=i,delay=randint(1,5))
        threads.append(t)
        t.start()
    print("Join threads ...")
    for i in range(10):
        print(f"to join: {i}")
        threads[i].join() # block the program until this thread is finished
        print(f"{i} joined .")
if __name__ == "__main__":
    main()
```

- Multithreading with Lock

```python
#!/usr/bin/env python
from threading import Lock, Thread
glob = 0
lock = Lock()
def increment():
    global glob
    for i in range(1000000):
        lock.acquire()
        glob += 1
        lock.release()
def decrement():
    global glob
    for i in range(1000000):
        lock.acquire()
        glob -= 1
        lock.release()
def main():
    threads = []
    for i in range(10):
        if i%2==0:
            target = increment
        else:
            target = decrement
        t = Thread(target=target)
        threads.append(t)
    for i in range(10):
        threads[i].start()
    for i in range(10):
        threads[i].join()
    print(glob)
if __name__ == "__main__":
    main()

```

- A multithreading thread example with semaphore

```python
#!/usr/bin/env python
from threading import Thread, Semaphore, active_count
from time import sleep
from random import randint
n_active = 0
semaphore = Semaphore(5)
def func(i):
    global n_active
    print(f"thread {i} is accquiring a semphore ...")
    semaphore.acquire()
    print(f"thread {i} accquired a semphore ...")
    n_active += 1
    sleep(randint(1,3))
    semaphore.release()
    n_active -= 1
def main():
    threads = []
    for i in range(20):
        t = Thread(target=func,args=(i,))
        threads.append(t)
    for i in range(20):
        threads[i].start()
    for i in range(20):
        print(f"the number of activated threads is {n_active}")
        threads[i].join()
if __name__ == "__main__":
    main()
```

- API for [multiprocessing](https://docs.python.org/3/library/multiprocessing.html) is almost same as multithreading

```python
#!/usr/bin/env python
from multiprocessing import Process
from time import sleep
from random import randint
def function(i):
    print(f"{i} start")
    sleep(randint(1,5))
    print(f"{i} stop")
def main():
    processes = []
    for i in range(10):
        p = Process(target=function,args=(i,))
        processes.append(p)
        p.start()
    print("Join threads ...")
    for i in range(10):
        print(f"to join: {i}")
        processes[i].join() # block the program until this thread is finished
        print(f"{i} joined .")
if __name__ == "__main__":
    main()
```

- Using `ProcessPoolExecutor` in [concurrent.futures](https://docs.python.org/3/library/concurrent.futures.html)

```python
#!/usr/bin/env python
from concurrent.futures import ProcessPoolExecutor, as_completed
import time
def count(number) :
    for i in range(0,10000000):
        i=i+1
    return i * number
def main():
    global numbers
    numbers = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    start = time.time()
    with ProcessPoolExecutor(max_workers=5) as executor:
        futures = [executor.submit(count, number) for number in numbers]
        for future in as_completed(futures):
            print(future.result())
    end = time.time()
    print(f"{end-start} seconds passed.")
if __name__ == "__main__":
    main()
```

- For more examples, see [python-parallel-programmning-cookbook](https://python-parallel-programmning-cookbook.readthedocs.io/zh_CN/latest/index.html)

- Parallel loading of gzipped files
  - see <https://www.blopig.com/blog/2016/08/processing-large-files-using-python/>
  - See cutadapt's implementation <https://github.com/marcelm/cutadapt/blob/main/src/cutadapt/pipeline.py>

### R
  - The parallel package <https://stat.ethz.ch/R-manual/R-devel/library/parallel/doc/parallel.pdf>
  - See this tutorial <https://bookdown.org/rdpeng/rprogdatascience/parallel-computation.html>
