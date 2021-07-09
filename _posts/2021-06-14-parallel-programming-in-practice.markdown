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

<https://www.gnu.org/software/parallel/>

### Paralleling in workflow manager

- gnu make

- snakemake


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


- R


## Inter process communication

- The named pipe or fifo
- If several processses simultaneously write to fifo, there is a chance that the data get interleaved

- <https://unix.stackexchange.com/questions/68146/what-are-guarantees-for-concurrent-writes-into-a-named-pipe>



## Example
- Parallel loading <https://github.com/marcelm/cutadapt/blob/main/src/cutadapt/pipeline.py>