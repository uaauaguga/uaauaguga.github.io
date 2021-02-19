## Simulate reactivity with SHAKER
- Preparation of SHAPE reacitivities with reference structure
  - PATTERNA set
  - Shaker set
  - Weeks 2008 set


- Enter SHAKER environment
  - `conda activate shaker-env-py27`
  ```{bash} 
  scripts/shaker-fit.py --dot-bracket data/ShaKer/data/RNA16.dbn --reactivity data/ShaKer/data/RNA16.react -m RNA16.pkl
  scripts/simulate-reactivity.py -i data/ShaKer/data/RNA20.dbn -m RNA16.pkl -o RNA20-prected.react
   ```
