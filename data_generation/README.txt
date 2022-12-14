* This is a rough draft README specifying how to create input data for the cgm_world quasar 
  absorption line neural network. *

-----------------------------------------------------------------------------------------------------

1.  Copy program files 

	Navigate to /fs1/project/cgm_world/code/python and copy the two program files located there 
	into your directory of choice. 

2.  Run 'cubegen.py' 

	This generates our data cube by setting our parameters. It creates two files: labels.txt 
	(each line has the velocity, column density, and b parameter) and noise_log.txt (contains 
	the signal-to-noise ratio and two negative integers used to set random seed noise character-
	istics in specsynth). The command is:

	apptainer exec --nv "<PATH>" python cubegen.py <NUM>
		- <PATH> is the path to your .sif file
		- <NUM> is the number of examples you want to generate

3.  Copy supporting input files

	Navigate to /fs1/project/cgm_world/code/fortan/input_files/ and copy the following files into
	your local directory: ions.table (specifies what sigma level to detect each transition), 
	sysanal.inp (points to atomic data), and zabs.dat (gives the system redshift). 

	You must adjust the "R" factor in sysanal.inp. This defines the resolution of the spectro-
	graph's instrumental spread function.

	sysanal.inp is currently configured to mimic HIRES and UVES.

4.  Run 'specgen_1cloud_same_init.py'

	This code reads in the text file, generates two spectrum files, runs sysanal to analyze 
	where the absorption regions are, and runs rebin spec on both spectrum files to clean the 
	continuum and re-bin it onto the particular velocity grid that is used. The command is:

	apptainer exec --nv "<PATH>" python specgen_1cloud_same_init.py
		- <PATH> is the path to your .sif file
