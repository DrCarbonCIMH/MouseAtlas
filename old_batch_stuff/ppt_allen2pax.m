function ppt_allen2pax( inputF, inputInterp )
% Transform image from the Allen mouse space to the Paxinos mouse space
% FORMAT ppt_allen2pax( P, interp )
% P            - a char array of input image filenames
% interp   - interpolation hold (see spm_slice_vol)
%                 [defaults (missing or empty) to 4 - 4th Degree B-Spline]
% The output file is in the same folder as the input file, with prefix 'pax'


%% Control

switch nargin
    case 1
        inputInterp = 4;
end

shiftMatrix = spm_matrix( [-5, -30, -55] ); % right -5 forward -30 up -55
d=fileparts(which('ppt_allen2pax'));
deformationFile = [d filesep 'y_shift2pax.nii'];

inputV = spm_vol( inputF );
inputP = spm_imatrix( inputV(1).mat );


%% Shift

[~, defReorient] = cfg_util( 'harvestdef', 'spm.util.reorient' );

defReorient.srcfiles = {inputF};
defReorient.transform.transM = shiftMatrix;
defReorient.prefix = 'ax';


%% Normalize

[~, defNormalizeWrite] = cfg_util( 'harvestdef', 'spm.spatial.normalise.write' );

defNormalizeWrite.subj.def = {deformationFile};
defNormalizeWrite.subj.resample = ...
    cfg_dep('Reorient Images: Reoriented Images', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));

defNormalizeWrite.woptions.bb = [-50 -95 -80; 50 60 10];
%         defNormalizeWrite.woptions.vox = [1 1 1];
defNormalizeWrite.woptions.vox = round( abs( inputP(7:9) ) * 10 ) / 10; % use the original spacing
defNormalizeWrite.woptions.interp = inputInterp;
defNormalizeWrite.woptions.prefix = 'p';


%% Batch

matlabbatch = cell( 2, 1 );
matlabbatch{1}.spm.util.reorient = defReorient;
matlabbatch{2}.spm.spatial.normalise.write = defNormalizeWrite;

strBatch = gencode( matlabbatch );

oneBatchFile = strrep( inputF, '.nii', '_allen2pax_batch.m' );
fid = fopen( oneBatchFile, 'w' );
fprintf( fid, '%s\n', strBatch{:} );
fclose( fid );

spm_jobman( 'run', oneBatchFile );


%% Cleanup

[inputPath, inputName, inputExt] = spm_fileparts( inputF );

reorientFile = fullfile( inputPath, [defReorient.prefix, inputName, inputExt] );

delete( reorientFile )
