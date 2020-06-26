This folder contains example data and MATLAB scripts to create a correlation matrix 
Files in data folder (unzip before use!!):
    - Atlas_inPax.nii = atlas created by app
    - Atlas.txt = text file combining numeric values of atlas with ROI names
    - FunctionalData.nii = example 4D dataset of functional-MRI dataset
    - AnatomicalData.nii = T2-weighted anatomical image (not needed; for 'illustration' purposes only)

The script "do_all.m" combines all necessary steps to extract ROI time courses, calculate and plot the correlation matrix.
Change the MATLAB directory to this folder and run "do_all".

In ExtractTimecourse.m the atlas is converted to the resolution of the functional data (named 'atlas_func.nii').
