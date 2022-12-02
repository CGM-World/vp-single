HOW TO RUN FILE CONVERSION

This code exists to convert files from .txt to .h5 format. In order to run it, you will need to:

- Edit the 'conversion.py' file. 
	- Put the path to the directory where your file lives as the 'datadir' variable.
	- Put your file name (without the extension) as the 'fi' variable.

	* Note: The program currently only works for numerical data but can easily be changed *

- Run the program. In order to use the needed libraries, we must utilize the apptainer. The way to run the script, then, is with this command in your terminal:
	- sbatch conversion_script.sh	 
