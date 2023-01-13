## ALPIDE-analysis
Framework for analysis of R3B & FRS EC experimental data involving the ALPIDE detectors.
To start, clone the repo and run make in the main directory.
Executables existing so far are
- clusterise
- analyse

Default arguments are in () brackets, and user arguments are in [].
```sh
./clusterize --file=FILE_IN.root --output=FILE_OUT.root --veto=[N](1) --firstEvent=[fE](0) --max-events=[mE](-1) --dets=[{x_i}]()

--file=FILE_IN.root         ..Select file.
--first-event=N             ..Start from N-th event. Default 0. 
--max-events=N              ..Specify maximum number of events. Default all entries.
--veto=N                    ..Only consider clusters with size>N. Default 1.
--output=/PATH/TO/OUT.root  ..Write clusterised data in this file. Default $(FILE_IN)_cl.root.
--dets=[d1,d2,..]           ..Condition to only write events which have clusters in every specified detectors.
      =all				    ..Equivalent to dets=1,2,..,ALPIDE_NUM. Every event must contain a cluster in all detectors.
--help                      ..Print this message to stdout.

The exe will cluster all the hits in a (selected) format and write an output root file.
A cluster is represented as a (float,float,uint)
corresponding to (meanX,meanY,N) of the cluster. N is the size of the cluster.
```

```sh
./calibrate --file=FILE_IN.root --output=FILE_OUT.root --firstEvent=[fE](0) --max-events=[mE](-1)
		\
--file=FILE_IN.root         ..Input file.
--first-event=N             ..Start from N-th event. Default 0.
--max-events=N              ..Specify maximum number of events. Default all entries.
--output=FILE_OUT.root		..Specify output file name. Default same as input file with 'calib' suffix.
--help                      ..Print this message to stdout.

The exe will calibrate the detectors (dXi,dYi) relative to ALPIDE1 and write an output root file.
Calibration of (dx, dy) for each alpide with reference to the first ALPIDE is stored in the
aCol,bCol,aRow,bRow and corresponding sigma branches. 'a' is the slope of the fit line and 'b' the offset,
in units of rows & columns.
```

```sh
./analyse --file=FILE_IN.root --output=FILE_OUT.root --firstEvent=[fE](0) --max-events=[mE](0) --dets=[{x_i}]()
--file=inputName.root        --Input file.
--first-event=N              ..Start from N-th event. Default 0 
--max-events=N	             ..Specify maximum number of events. Default all entries.

--hitmap=X                   ..Plots the hitmap of AlpideX. Works for files containing raw or clustered data.
--hitmap		             ..Plots the hitmaps of all Alpides. Works for files containing raw or clustered data.
--track			             ..Do the tracking.
	--cal=calFile.root       ..Pass calibration file to the tracker (OPTIONAL).
	--save=saveFile.root     ..Save the tracks into a rootfile. If not specified, only plots data.
--help                       ..Print this message to stdout.

This exe can plot the hitmaps of all or specified detector(s). Works with either raw or clustered data.
To do tracking pass --track --cal=<calFile>.root --save=<saveFile>.root flags. Works only for input files with clustered data.
This will create a rootfile with (X,Y,Z) branches containing (calibrated) x,y,z hit positions that the algorithm recognized 
form a track. Track can only be formed with 3 or more hit points.
``` 
# ALPIDE-analysis
