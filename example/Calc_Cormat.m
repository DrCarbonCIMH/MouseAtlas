function [cormat, names]=Calc_Cormat(tc_struc)

% calculate set of windowed connectivity matrices
%
% tc_struc = structure with regional timecourses (output of ExtractTimecourse.m)

Nsub=size(tc_struc,2);
names=fieldnames(tc_struc);
Nreg=length(names);
Nimg=length(tc_struc(1).(names{1}));

for jsub=1:Nsub
    tc_mat=zeros(Nreg,Nimg);
    for jreg=1:Nreg
        tc_mat(jreg,:)=tc_struc(jsub).(names{jreg});
    end
    cormat{jsub}=corrcoef(tc_mat', 'rows', 'pairwise' );
end
