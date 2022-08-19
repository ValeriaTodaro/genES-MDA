
						 University of Parma, Italy
							     &
					  Universitat Politènica de València, Spain

						      Valeria Todaro 
						         May 2022               


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
	Case1	: solution of the reverse flow routing problem for the estimation of the inflow hydrograph to a river reach based on observed water levels in a downstream section. 
                  Case 1 uses HEC-RAS v.5.0.7 as forward model and it must be installed on the computer. 
		  HEC-RAS v.5.0.7 is freely downloadable from the USACE Hydrologic Engineering Center's website (https://www.hec.usace.army.mil/software/hec-ras/the).

	Case2	: estimation of a hydraulic conductivity field using piezometric observations. 
		  Case 2 uses MODFLOW2005 as forward model and its excutable file must be present in the Model folder. 
		  For your convinience, the executable file mf2005dbl.exe is already provided in the Model folder. 
		  
	Case3	: identification of the release history of contaminant spills in an aquifer based on observed concentration data. 
		  Case 3 uses MODFLOW2005 and MT3D-USGS as forward models and their excutable files must be present in the Model folder. 
		  For your convinience, the executable files mf2005dbl.exe and mt3d-usgs_1.1.0_64.exe are already provided in the Model folder. 

To execute one of the application example, the python command ES_MDA shoul be run using the default settings.


For more information, and to cite this tool, please refer to the following paper:

Todaro, V.; D'Oria, M.; Tanda, M.G.; Gómez-Hernández, J.J. genES-MDA: a generic open-source software package to solve inverse problems via the Ensemble	
Smoother with Multiple Data Assimilation. Computers & Geosciences, 2022, 167, 105210. https://doi.org/10.1016/j.cageo.2022.105210.
