%% define files
% functional data set
Pfunc = [pwd filesep 'data' filesep 'FunctionalData.nii']
% atlas
Patlas = [pwd filesep 'data' filesep 'Atlas_inPax.nii']
% atlas txt file
Ptxt = [pwd filesep 'data' filesep 'Atlas.txt']

%% extract time courses
tc_struc=ExtractTimecourse(Patlas,Ptxt,Pfunc);

%% calculate correlation matrix
[cormat, names]=Calc_Cormat(tc_struc);

%% plot correlation matrix
Plot_cormat(cormat{1},names)