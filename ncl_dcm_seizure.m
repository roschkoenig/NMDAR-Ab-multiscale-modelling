function ncl_dcm_seizure(F)

fs  = filesep;

% Load control DCMs (initial inversion)
%--------------------------------------------------------------------------
CNT = load([F.outp fs 'DCM' fs 'ICM_Control baseline']);
CNT = rmfield(CNT.DCM, 'name');

% Load NMDA-baseline DCMs (from Bayesian model comparison)
%-------------------------------------------------------------------------
NMD = load([F.outp fs 'DCM' fs 'Model_comparison_NMDAR.mat']); 
clear FE
for n = 1:length(NMD.P), FE(n) = NMD.P{n}.F; end
[val winner] = max(FE); 
NMD = NMD.P{winner}; 

% Load (uninverted) seizure DCMs 
%--------------------------------------------------------------------------
SZR = load([F.outp fs 'DCM' fs 'DCM_NMDAR-Ab positive seizure.mat']); 
SZR = rmfield(SZR.DCM, 'name');

% Reinvert seizure effect using NMD average as priors
%==========================================================================
SZR.M.pE    = NMD.Ep; 
SZR         = ncl_spm_dcm_csd(SZR); 

% Estimate main effect of seizure using PEB
%==========================================================================
PCM{1,1} = NMD; 
PCM{2,1} = SZR; 
M.X(:,1) = [1, 1];      % group mean
M.X(:,2) = [0, 1];      % main effect of seizure 
M.Xnames = {'Mean', 'Seizure'};

[PEB, PCM] = spm_dcm_peb(PCM, M, {'T', 'G', 'S'}); 
REB        = spm_dcm_peb_bmc(PEB); 

% Extract the main effect of seizure
%--------------------------------------------------------------------------
pid = find(abs(REB.Ep(length(REB.Pnames)+1:length(REB.Pnames)*2)) > 0.001);
Sfx.Pnames = REB.Pnames(pid);
Sfx.pid    = REB.Pind(pid);
Sfx.vec    = REB.Ep(length(REB.Pnames)+pid); 

% Pack up for export
%--------------------------------------------------------------------------
clear AllCM
AllCM{1} = CNT; 
AllCM{2} = NMD; 
AllCM{3} = SZR;

save([F.outp fs 'DCM' fs 'Ictogenesis'], 'Sfx', 'REB', 'AllCM')
