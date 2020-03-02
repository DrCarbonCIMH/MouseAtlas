%% old non-working
% matlabbatch{1, 1}.spm.util.reorient.srcfiles = {'/data/Markus/TestAllen/Test1.nii'};
% matlabbatch{1, 1}.spm.util.reorient.transform.transM = [1 0 0 -5
%                                                         0 1 0 -30
%                                                         0 0 1 -55
%                                                         0 0 0 1];
% matlabbatch{1, 1}.spm.util.reorient.transform.transprm = '<UNDEFINED>';
% matlabbatch{1, 1}.spm.util.reorient.transform.transF = '<UNDEFINED>';
% matlabbatch{1, 1}.spm.util.reorient.prefix = 'ax';
% matlabbatch{2, 1}.spm.spatial.normalise.write.subj.def = {'/data/matlab/ms_AllenMouse_v3/20171030_Allen2Pax/y_shift2pax.nii'};
% % matlabbatch{2, 1}.spm.spatial.normalise.write.subj.resample(1) = cfg_dep('Reorient Images: Reoriented Images', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
% matlabbatch{2, 1}.spm.spatial.normalise.write.woptions.bb = [-50 -95 -80
%                                                              50 60 10];
% matlabbatch{2, 1}.spm.spatial.normalise.write.woptions.vox = [0.5 0.5 0.5];
% matlabbatch{2, 1}.spm.spatial.normalise.write.woptions.interp = 0;
% matlabbatch{2, 1}.spm.spatial.normalise.write.woptions.prefix = 'p';

function [matlabbatch] = ms_allen2pax_batch(inputF)
%% Shift
inputV = spm_vol( inputF );
inputP = spm_imatrix( inputV(1).mat );
shiftMatrix = spm_matrix( [-5, -30, -55] ); % right -5 forward -30 up -55
[~, defReorient] = cfg_util( 'harvestdef', 'spm.util.reorient' );

defReorient.srcfiles = {inputF};
defReorient.transform.transM = shiftMatrix;
defReorient.prefix = 'ax';

%% OldNorm
deformationFile = 'shift2pax_sn_20171221.mat';
[~, defOldnormWrite] = cfg_util( 'harvestdef', 'spm.tools.oldnorm.write' );

defOldnormWrite.subj.matname = {[fileparts(which('ms_allen2pax_batch')) filesep deformationFile]};
defOldnormWrite.subj.resample = ...
    cfg_dep('Reorient Images: Reoriented Images', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));

defOldnormWrite.roptions.preserve = 0;
defOldnormWrite.roptions.bb = [-50 -95 -80; 50 60 10];
defOldnormWrite.roptions.vox = round( abs( inputP(7:9) ) * 10 ) / 10; % use the original spacing
defOldnormWrite.roptions.interp = 0;
defOldnormWrite.roptions.wrap = [0 0 0];
defOldnormWrite.roptions.prefix = 'p';

%% Batch

matlabbatch = cell( 2, 1 );
matlabbatch{1}.spm.util.reorient = defReorient;
% matlabbatch{2}.spm.spatial.normalise.write = defNormalizeWrite;
matlabbatch{2}.spm.tools.oldnorm.write = defOldnormWrite;
