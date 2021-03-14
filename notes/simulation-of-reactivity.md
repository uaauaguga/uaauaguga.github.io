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
  scripts/fit-genextreme.py -d data/reactivity/patteRNA-weeks.txt -s data/reactivity/patteRNA-weeks-genextreme-statistics.txt -m data/reactivity/patteRNA-weeks-genextreme-model.pkl
  scripts/simulate-genextreme.py -i test/predicted/RNAstructure-no-shape.dot -o test/predicted/RNAstructure-no-shape.txt -m data/reactivity/patteRNA-weeks-genextreme-model.pkl
  ```
