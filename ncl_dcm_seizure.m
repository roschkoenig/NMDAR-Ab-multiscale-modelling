fs  = filesep;
CNT = load([F.outp fs 'DCM' fs 'DCM_Control baseline']);
NMD = load([F.outp fs 'DCM' fs 'Model_comparison_NMDAR.mat']); 
SZR = load([F.outp fs 'DCM' fs 'DCM_NMDAR-Ab positive seizure.mat']); 

%% Reinvert seizure effect using NMD average as priors
%--------------------------------------------------------------------------
szr         = rmfield(SZR.DCM, 'name'); 
szr.M.pE    = NMD.AVG.Ep; 
szr         = ncl_spm_dcm_csd(szr); 


% Estimate main effect of seizure using PEB
%--------------------------------------------------------------------------

% Apply the seizure effect to control setup - stability analysis 
%--------------------------------------------------------------------------
