## Simulate reactivity with SHAKER
- Preparation of SHAPE reacitivities with reference structure
  - PATTERNA set
  - Shaker set
  - Weeks 2008 set


- Enter SHAKER environment
  - `conda activate shaker-env-py27`
  ```{bash}
  scripts/shaker-loo.py --input data/reactivity/patteRNA-weeks.txt --performance test/reactivity-simulation/shaker-performance.txt --reactivity test/reactivity-simulation/shaker-loo.shape
  ```
  - Under infernal env
  ```{bash} 
  scripts/fit-genextreme.py --dataset data/reactivity/patteRNA-weeks.txt --statistics patteRNA-weeks-fit-statistics.txt --model patteRNA-weeks-genextreme.pkl
  ```
