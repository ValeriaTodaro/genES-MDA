*****************************************************************************************************************************************
							 University of Parma, Italy
								    &
						   Universitat Politènica de València, Spain

							     Valeria Todaro 
							       May 2022               
*****************************************************************************************************************************************


genES-MDA is an open-source Python software package for the solution of inverse
problems by means of the Ensemble Smoother with Multiple Data Assimilation (ES-MDA).
It can be used with any Python-based platform and does not require an installation. 

This directory contains the genES-MDA software and the Examples folder.
All files in the genES-MDA package need to be downloaded to the same folder preserving the directory hierarchy.


genES-MDA contains the following files:

	ESMDA.py	 : main module, it does not require modifications by the user.
	Mod.py 		 : subordinate module, it requires modifications by the user.
	InputSettings.py : subordinate module, it requires modifications by the user.
	Tools 		 : Python package it does not require modifications by the user.
	Model		 : folder that contains the forward model files, it requires modifications by the user.
	Obs.txt  	 : guidelines to set up the text input file that contains observations.	
	Par.txt 	 : guidelines to set up the text input file that contains parameters.
	Ens.txt		 : guidelines to set up the text input file that contains the initial ensemble of parameters (optional).
	Errors.txt	 : guidelines to set up the text input file that contains the ensemble of measurement errors (optional).
	R.txt 		 : guidelines to set up the text input file that contains the covariance matrix of measurement errors (optional).

Examples contains the following case study applications:
	Case1	: solution of the reverse flow routing problem for the estimation of the inflow hydrograph to a river reach based on 
                  observed water levels in a downstream section.
	Case2	: estimation of a hydraulic conductivity field using piezometric observations.
	Case3	: identification of the release history of contaminant spills in an aquifer based on observed concentration data.

To execute one of the application example, replace the file contained in genESMDA with those reported in Case1, Case2 or Case3 folders.
