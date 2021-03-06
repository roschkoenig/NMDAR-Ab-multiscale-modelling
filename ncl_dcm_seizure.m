function ncl_dcm_seizure(F)
dopeb   = 0; 
fs      = filesep;

% Load control DCMs (initial inversion)
%--------------------------------------------------------------------------
CNT = load([F.outp fs 'DCM' fs 'ICM_Control baseline']);
CNT = rmfield(CNT.DCM, 'name');

% Load NMDA-baseline DCMs (from Bayesian model comparison)
%-------------------------------------------------------------------------
NMD = load([F.outp fs 'DCM' fs 'Model_comparison_NMDAR.mat']);
AVG = NMD.AVG; 
selection = 'winner';   % 'average' or 'winner' or 'full'

if strcmp(selection, 'average')
    pC  = NMD.P{1}; 
    NMD = NMD.AVG; 
    NMD.M.pC = pC; 
    
elseif strcmp(selection, 'winner')
    clear FE
    for n = 1:length(NMD.P), FE(n) = NMD.P{n}.F; end
    [val winner] = max(FE); 
    NMD = NMD.P{winner};
    
elseif strcmp(selection, 'full')
    NMD = NMD.P{1}; 
end

% Load (uninverted) seizure DCMs 
%--------------------------------------------------------------------------
SZR = load([F.outp fs 'DCM' fs 'DCM_NMDAR-Ab positive seizure.mat']); 
SZR = rmfield(SZR.DCM, 'name');

% Reinvert seizure effect using NMD average as priors
%==========================================================================
SZR.M.pE    = AVG.Ep; 
SZR         = ncl_spm_dcm_csd(SZR); 

if dopeb 
    
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

else    
    Sfx.pid = [spm_fieldindices(CNT.Ep, 'G'); 
               spm_fieldindices(CNT.Ep, 'T'); 
               spm_fieldindices(CNT.Ep, 'S')]; 
    nmdvec  = spm_vec(NMD.Ep); 
    szrvec  = spm_vec(SZR.Ep); 
    Sfx.vec = szrvec(Sfx.pid)-nmdvec(Sfx.pid);

end

% Pack up for export
%--------------------------------------------------------------------------
clear AllCM
AllCM{1} = CNT; 
AllCM{2} = NMD; 
AllCM{3} = SZR;

save([F.outp fs 'DCM' fs 'Ictogenesis'], 'Sfx', 'AllCM')
