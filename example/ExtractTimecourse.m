function tc_struc=ExtractTimecourse(Patlas,Ptxt,P)

% extract regional timecourses for ROIS defined by Ptxt and atlas from .nii-files 
%
% input: Ptxt = path-string for txt-file with roi-names
%           P = cell array of path-strings for nii-files to analize
%           atlas = path of atlas region definition nii-file
%
% output: tcdata = cell array containing regional timecourses 
%          names = vector of roi names


if nargin<1
    Patlas=spm_select(1,'image','Select atlas structure file');
end    

if nargin<2
	[Ptxt,sts]=spm_select(1,'any','Select Textfile',[],pwd,'.*.txt');
end
    
if nargin<3
	[Pchar,sts]=spm_select(inf,'image','Select Images');
    P={Pchar};
end
fprintf('Extracting time courses\n')
% Make a cellarray out of the Filenames
if ~iscell(P), P=cellstr(P);end
nsubs=length(P);
regfile=read_t2s(Ptxt);
nreg=size(regfile,1);

% read region names
for lnum=1:nreg
     tline=strtrim(regfile(lnum,:));
     %format tline: 'RegionName \t #  '
     cix=regexp(tline,'\t');
     region(lnum).name=strtrim(tline(1:cix-1));
     %Kick out ',' that sometimes appear in names
     region(lnum).name=strrep(region(lnum).name,',','');
     %Kick out '-' that sometimes appear in names
     region(lnum).name=strrep(region(lnum).name,'-','_');
     %Kick out '''  that sometimes appear in names
     region(lnum).name=strrep(region(lnum).name,'''', '');
     try
        eval(['nums=[' tline(cix:end) '];']);
     catch
         keyboard
     end
     region(lnum).nums=nums;
end


% Loop over functional Files (per subject)
for ns=1:nsubs
    %Create atlas image in space of functional image:
    Pcur=P{ns};
    [fpath fname ext]=fileparts(Pcur);
    Pcur=[fpath filesep fname '.nii'];  
    V1=spm_vol([Pcur ',1']);
    V2=spm_vol(Patlas);
    Vi=[V1 V2];
    Vo=V1;
    % Atlas in functional resolution is written out in current directory.
    Vo.fname=[pwd filesep 'atlas_func.nii'];
    Vatlas=spm_imcalc(Vi,Vo,'i2',{0, 0, 0});        
    Vcur=spm_vol([fpath filesep fname '.nii']);
    nimg=size(Vcur,1); %Number of repetitions in functional file
    
    % regional partition according to atlas
    
    tcourse=zeros(nreg,nimg);
    mtx_func=spm_read_vols(Vcur);
    funcsize=size(mtx_func);
    mtx_func=reshape(mtx_func,prod(funcsize(1:3)),funcsize(4));
    mtx_atl=spm_read_vols(Vatlas);
    clear V1 V2 Vcur
    fprintf('%s: ',fname);
    for nr=1:nreg
        fprintf('%d ',nr);
        if 1==1
            rmask=ismember(mtx_atl,region(nr).nums);
            %Uncomment following to write out single regional masks in
            %functional resolution
            %Vout=Vatlas;
            %Vout.fname=[region(nr).name '.nii'];
            %spm_write_vol(Vout,rmask);
            indx=find(rmask);
            size(indx,1);
            roidat=mtx_func(indx,:);
            nanIndex = any( isnan( roidat ), 2 ); %Kick out any NANs that might be in the region
            meantime = squeeze( mean( roidat(~nanIndex, :) ) );
            try
                tc_struc(ns).(region(nr).name)=meantime;
            catch
                keyboard
            end
            clear roidat;
        end
    end
    fprintf('\n');
end
