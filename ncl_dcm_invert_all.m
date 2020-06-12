function ncl_dcm_invert_all(F)
% Load DCM files 
%--------------------------------------------------------------------------
fs      = filesep;
D       = spm_eeg_load(F.spm_file); 
conds   = condlist(D);  
clear DCM
for c = 1:length(conds)
    T = load([F.outp fs 'DCM' fs 'DCM_' conds{c}]);
    DCM(c) = T.DCM; 
end

% Invert individual DCMs without altered priors
%--------------------------------------------------------------------------
clear ICM
for c = 1:length(conds)
    if ~strcmp(conds{c}, 'Control baseline')
        TMP = load([F.outp fs 'DCM' fs 'DCM_Control baseline']);
        DCM(c).xY.scale = TMP.DCM.xY.scale;
    end
    DCM(c).M.pE.L = 1000; 
    ICM(c) = ncl_spm_dcm_csd(DCM(c)); 
end

% Plot initial model fits
%--------------------------------------------------------------------------
clf
cols = ['k', 'b', 'r']; 
for c = 1:length(conds)
    plot(log(abs(ICM(c).Hc{1})), 'color', cols(c)), hold on
    plot(log(abs(ICM(c).xY.y{1})), ':', 'color', cols(c));  hold on
end

