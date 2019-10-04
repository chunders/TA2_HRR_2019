# Inteferometry Analysis
**Author:** Chris Underwood  
**Email:** christopher.underwood@york.ac.uk

Interferometry analysis takes an interferogram and returns the density.  
This code has been written in: `python 3`


## Code Structure
1. Import data to numpy array.
2. Get reference
	* Rotate fringes to vertical and create reference if required
	* Load reference file.
3. Remove background, and crop to ROI.
5. Recreate phase.
6. Unwrap phase.
7. Abel transform and return density.

## Code structure
The code has been seperated into seperate files for each of the different stages.
The main file to run is:
`Interferometry_Extraction.py`

There are the following classes:

* `backgroundRemover_class.py`
* `createPhase_class.py`
* `loadDataToNumpy_class.py`
* `test_createPhase_class.py`
* `createDensity_class.py`
* `createRefenceFromImage_class.py`
* `test_createDensity_class.py`
* `unwrapPhase_class.py`
* `fourier_mask_class.py`

The documentation is in each file.

## Interactive notebook to get fitting parameters
A jupyter notebook file exists to find the ideal parameters for extraction.
This will save the required parameters to the database (see next section).


### Experimental Data Database
There is the option of an experimental data base of the paramaters required for extraction.
Then to run in batch mode the extraction parameters can be loaded from this database.

