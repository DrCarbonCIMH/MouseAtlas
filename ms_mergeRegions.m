function [ output_args ] = ms_mergeRegions( atlas, toMerge, newRegion )
%MS_MERGEREGIONS Merge brain regions of an existing atlas
%   Usage: ms_mergeRegions( atlas, toMerge, newRegion )
% atlas: nii-file of the atlas
% toMerge: cellarray containing the region numbers to be merged
% newRegion: nx2 cellarray containing the number and name of the new ('merged') brain region
%
% Example: 
% toMerge = {[477 549 512]; [313 698]};
% newRegion = {12, 'newRegion'; 18, 'newRegion2'};
% will combine the regions 477, 549, and 512 to 'newRegion' with number 12
% and regions 313 and 698 to 'newRegion2' with number 18
% If you do not assign any numbers to the new region it will get the first number
% of the corresponding merged regions!
% newRegion = {'newRegion'; 'newRegion2'}; -> regions 477, 549, and 512 to
% 'newRegion' with number 477 AND regions 313 and 698 to 'newRegion2' with
% number 313
    
[d, name, ~] = fileparts(atlas);
%% get the text file
txtFile = [d filesep name '.txt'];
fid=fopen(txtFile, 'r'); C=textscan(fid, '%s %f'); fclose(fid); % C{1} are the names; C{2} the numbers


%% do some checking
if size(newRegion,2)==1 % that means no number is assigned to the new regions
    for ix=1:length(newRegion)
        tmp{ix,1} = toMerge{ix}(1); tmp{ix,2}=newRegion{ix};
    end
    newRegion=tmp;
else
    for ix=1:size(newRegion,1) % check if the number of newRegion is already in the atlas
        if any(C{2}==newRegion{ix,1}); warndlg(['The new Region ' newRegion{ix,2} ' has the same number as ' C{1}{C{2}==newRegion{ix,1}} '. I hope that is okay.']); end
    end
end
for ix=1:length(toMerge)
    fprintf('Regions ');
    for iy=1:length(toMerge{ix})
        if ~any(C{2}==toMerge{ix}(iy)); errordlg(['There is no Regionnumber ' num2str(toMerge{ix}(iy))]); return; end
        fprintf('%s ', C{1}{C{2}==toMerge{ix}(iy)});
    end
    fprintf('will be combined to %s (%i)\n', newRegion{ix,2}, newRegion{ix,1})
end
%% change the matrix and remove the entries of the txt-file
V=spm_vol(atlas);
Mtx=spm_read_vols(V);

for ix=1:length(toMerge)
    for iy=1:length(toMerge{ix})
        Mtx(Mtx==toMerge{ix}(iy))=newRegion{ix,1};
        C{1}(C{2}==toMerge{ix}(iy))=[];
        C{2}(C{2}==toMerge{ix}(iy))=[];
    end
end

%% write the nii
Vout=V;
Vout.fname = [d filesep name '_merged.nii'];
fprintf('Writing new atlas %s \n', Vout.fname);
spm_write_vol(Vout, Mtx);

%% write the txt file
txtFileOut=[d filesep name '_merged.txt'];
for ix=1:size(newRegion,1); C{1}(end+1)=newRegion(ix,2); C{2}(end+1)=newRegion{ix,1}; end
fprintf('Writing new txt-file %s \n', txtFileOut);
fid=fopen(txtFileOut,'w');
for ix=1:length(C{1})-1
    fprintf(fid, '%s\t%i\n',C{1}{ix},C{2}(ix));
end
fprintf(fid, '%s\t%i',C{1}{ix+1},C{2}(ix+1));
fclose(fid);

end

