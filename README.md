# A More Practical Algorithm for Weighted First-Order Model Counting with Linear Order Axiom


## Installation
Install the package:
```
$ pip install -e .
```

## Run the experiments

You can find the sentence files at the directory 'models'.

### Three-way split
Run the following command:
```
$ python wfomc_fo2/ecai_wfomc.py -i models/three-way.wfomcs -a inc
$ python wfomc_fo2/ecai_wfomc.py -i models/three-way.wfomcs -a rec
```

### Immediate predecessor
Run the following command:
```
$ python wfomc_fo2/ecai_wfomc.py -i models/predecessor.wfomcs -a inc
$ python wfomc_fo2/ecai_wfomc.py -i models/predecessor.wfomcs -a rec
```

### BA(3) model
Run the following command:
```
$ python wfomc_fo2/ecai_wfomc_cc.py -i models/BA3.wfomcs -a inc
$ python wfomc_fo2/ecai_wfomc_cc.py -i models/BA3.wfomcs -a rec
```

### Smoking-drinking friends
Run the following command:
```
$ python wfomc_fo2/ecai_wfomc.py -i models/smoke.mln -a inc
$ python wfomc_fo2/ecai_wfomc.py -i models/smoke.mln -a rec_real
```

### Watts-Strogatz model
Run the following command:
```
$ python wfomc_fo2/ecai_wfomc.py -i models/WS.wfomcs -a inc
$ python wfomc_fo2/ecai_wfomc.py -i models/WS.wfomcs -a rec
```