function ncl_dcm_specify(F) 
clear DCM
fs      = filesep; 

% Extract specifics from SPM file of LFPs
%==========================================================================
D       = spm_eeg_load(F.spm_file); 
Fs      = fsample(D); 
smpls   = size(D,2); 
timax   = time(D); 
conds   = condlist(D); 

for cond_no = 1:length(conds)
% cond_no = find(strcmp(conds, 'Control baseline'));
% % 'Control baseline', 
% % 'NMDAR-Ab positive baseline', 
% % 'NMDAR-Ab positive seizure'

% Define frequency band of interest
%--------------------------------------------------------------------------
% 'broad', 'delta', 'theta', 'alpha', 'beta', 'gamma'

disp(['Working on condition: ' conds{cond_no}]); 
fband   = 'broad';
switch fband 
    case 'broad',   Fbp = [1 60];
    case 'delta',   Fbp = [1 4];
    case 'theta',   Fbp = [4 8];
    case 'alpha',   Fbp = [8 12];
    case 'beta',    Fbp = [12 30];
    case 'gamma',   Fbp = [30 60];
    case 'high gamma',  Fbp = [60 100];
    case 'alpha2gamma', Fbp = [8 100];
end

% Set up DCM details
%--------------------------------------------------------------------------
% Set up DCM details
%--------------------------------------------------------------------------
DCM                  = []; 
DCM.xY.Dfile         = F.spm_file; 
DCM.options.analysis = 'CSD';       % cross-spectral density 
DCM.options.model    = 'CMC';      	% structure canonical microcircuit (for now)
DCM.options.spatial  = 'LFP';           
DCM.options.Tdcm     = [timax(1) timax(end)] * 1000;     

DCM.options.Fdcm    = Fbp;          % frequency range  
DCM.options.D       = 1;         	% frequency bin, 1 = no downsampling
DCM.options.Nmodes  = 8;          	% cosine reduction components used 
DCM.options.h       = 0;         	% no hanning 
DCM.options.trials  = cond_no;   	% index of ERPs within file

DCM.Sname           = chanlabels(D);
DCM.M.Hz            = DCM.options.Fdcm(1):0.5:DCM.options.Fdcm(2);
DCM.xY.Hz           = DCM.M.Hz;

% Define different connection types
%==========================================================================
Fd      = zeros(size(D,1));    % Forward (empty)
Bd      = zeros(size(D,1));    % Backward (empty)
Sf      = zeros(size(D,1));   for s = 1:length(Sf); Sf(s,s) = 1; end

% Define model arcitecture (A), conditional effects (B) and input (C) 
%==========================================================================
DCM.A{1}    =   Fd + Bd;
DCM.A{2}    =   Fd + Bd;
DCM.A{3}    =   Sf; 
DCM.B       = 	{}; 
DCM.C       =   sparse(length(DCM.A{1}),0); 

% Reorganise model parameters in specific structure
%==========================================================================
DCM.M.dipfit.Nm    = DCM.options.Nmodes;
DCM.M.dipfit.model = DCM.options.model;
DCM.M.dipfit.type  = DCM.options.spatial;

DCM.M.dipfit.Nc    = size(D,1);
DCM.M.dipfit.Ns    = length(DCM.A{1});

% Prepare and load data
%=========================================================================
DCM  = spm_dcm_erp_data(DCM);                   % data
DCM  = spm_dcm_erp_dipfit(DCM, 1);       

% Define priors
%==========================================================================
% Load standard neural priors
%--------------------------------------------------------------------------
[pE,pC]  = ncl_spm_dcm_neural_priors(DCM.A,DCM.B,DCM.C,DCM.options.model);  
[pE,pC]  = spm_L_priors(DCM.M.dipfit,pE,pC);
[pE,pC]  = spm_ssr_priors(pE,pC);

DCM.M.pE   = pE;
DCM.M.pC   = pC;

% Load pre-specified noise covariance from baseline
%==========================================================================
if ~strcmp(conds{cond_no}, 'Control baseline')
for b = 1:max(length(DCM.B),1) 
    TMP = load([F.outp fs 'DCM' fs 'DCM_Control baseline.mat']); 
    DCM.xY.mar{b}.noise_cov = TMP.DCM.xY.mar{b}.noise_cov; 
end 
end

% Calculate Data features
%==================================================================
DCM = ncl_spm_dcm_csd_data(DCM);
DCM.options.DATA = 0;

% Save extracted data features and DCM structure
%--------------------------------------------------------------------------
DCM.name = [F.outp fs 'DCM' fs 'DCM_' conds{cond_no}]; 
save(DCM.name, 'DCM'); 

end