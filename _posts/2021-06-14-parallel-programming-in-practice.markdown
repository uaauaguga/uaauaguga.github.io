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

- [Preemption](https://en.wikipedia.org/wiki/Preemption_%28computing%29)


## External tools for paralleling

### GNU's parallel utility

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

```

### Paralleling in workflow manager

- gnu make

- [snakemake](https://snakemake.readthedocs.io/en/stable/)


## Paralleling in programming language
  
- C++

- python

  - IO bounded problem
    - multithreading
    - asyncio
  - CPU bounded problem
    - multiprocessing

<https://realpython.com/python-concurrency/>

<https://joblib.readthedocs.io/en/latest/>


  - Parallel loading of gzipped files
    - see <https://www.blopig.com/blog/2016/08/processing-large-files-using-python/>
    - See cutadapt's implementation <https://github.com/marcelm/cutadapt/blob/main/src/cutadapt/pipeline.py>

- R


## Inter process communication

- The named pipe or fifo
- If several processses simultaneously write to fifo, there is a chance that the data get interleaved

- <https://unix.stackexchange.com/questions/68146/what-are-guarantees-for-concurrent-writes-into-a-named-pipe>


