## ALPIDE-analysis
Framework for analysis of R3B & FRS EC experimental data involving the ALPIDE detectors.
To start, clone the repo and run make in the main directory.
Executables existing so far are
- clusterise
- analyse
```sh
./clusterize --file=FILE_IN.root --output=FILE_OUT.root --veto=[N](1) --firstEvent=[fE](0) --max-events=[mE](0) --dets=[{x_i}]()
file=/PATH/TO/file.root     Select file.
--first-event=N             Start from N-th event. Default 0. 
--max-events=N              Specify maximum number of events. Default all entries.
--veto=N                    Only consider clusters with size>N.
--output=/PATH/TO/OUT.root  Write clusterised data in this file. Default $(inputName)_coarse_cl.root.
                            Specify output file name. Default same as input file with cl suffix.
--dets=[d1,d2,..]           Condition to only write events which have clusters every specified detectors.
      =all				    Equivalent to dets=1, 2, ... ALPIDE_NUM. Every event must contain a cluster in all detectors.
--help                      Print this message to stdout. 

The exe will cluster all the hits in a (selected) format and write an output root file.
A cluster is represented as a tuple<float,float,uint> 
corresponding to (meanX, meanY, N) of the cluster. N is the size of the cluster.
```
```sh
./analyse --file=FILE_IN.root --output=FILE_OUT.root --veto=[N](1) --firstEvent=[fE](0) --max-events=[mE](0) --dets=[{x_i}]()
file=/PATH/TO/FILE/file.root
--first-event=N             Start from N-th event. Default 0. 
--max-events=N	            Specify maximum number of events. Default all entries.
--raw=X,Y                   Plots only raw correlations of AlpideX:AlpideY.
--raw=X                     Plots the raw hitmap of AlpideX.
--cal                       Will do the calibration process.
--help                      Print this message to stdout. 
``` 
# ALPIDE-analysis
