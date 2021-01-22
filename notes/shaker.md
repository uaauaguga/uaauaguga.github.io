### Install shaker for reactivity simulation
```{bash}
conda create -n shaker-env-py27 python=2.7
conda activate shaker-env-py27
conda install -c bioconda viennarna
pip install xgboost seaborn tabulate toolz
## The following two package get error if use pip to install
conda install osqp
conda install scs
pip install git+https://github.com/smautner/EDeN.git --user
pip install shaker-rna
```
