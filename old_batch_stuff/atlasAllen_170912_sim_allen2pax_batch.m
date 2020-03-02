matlabbatch{1, 1}.spm.util.reorient.srcfiles = {'/home/lei.zheng/work/output/20171030_Allen2Pax/atlasAllen_170912_sim.nii'};
matlabbatch{1, 1}.spm.util.reorient.transform.transM = [1 0 0 -5
                                                        0 1 0 -30
                                                        0 0 1 -55
                                                        0 0 0 1];
matlabbatch{1, 1}.spm.util.reorient.transform.transprm = '<UNDEFINED>';
matlabbatch{1, 1}.spm.util.reorient.transform.transF = '<UNDEFINED>';
matlabbatch{1, 1}.spm.util.reorient.prefix = 'ax';
matlabbatch{2, 1}.spm.spatial.normalise.write.subj.def = {'y_shift2pax.nii'};
matlabbatch{2, 1}.spm.spatial.normalise.write.subj.resample(1) = cfg_dep('Reorient Images: Reoriented Images', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
matlabbatch{2, 1}.spm.spatial.normalise.write.woptions.bb = [-50 -95 -80
                                                             50 60 10];
matlabbatch{2, 1}.spm.spatial.normalise.write.woptions.vox = [0.5 0.5 0.5];
matlabbatch{2, 1}.spm.spatial.normalise.write.woptions.interp = 4;
matlabbatch{2, 1}.spm.spatial.normalise.write.woptions.prefix = 'p';
