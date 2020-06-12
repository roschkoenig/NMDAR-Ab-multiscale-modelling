fs      = filesep;

% Load control
%--------------------------------------------------------------------------
control = dir([F.outp fs 'DCM' fs 'DCM_Control*']);
CNT     = load([control.folder fs control.name]);
CNT     = CNT.DCM; 

% Load NMDAR-baseline
%--------------------------------------------------------------------------
nmdar   = dir([F.outp fs 'DCM' fs 'DCM_NMDAR-Ab positive s*']); 
NMD     = load([nmdar.folder fs nmdar.name]); 
NMD     = NMD.DCM; 

% Switch off variation in observation parameters
%--------------------------------------------------------------------------
Ep      = CNT.Ep; 
Cp      = spm_unvec(diag(CNT.Cp), CNT.M.pC);

pE      = Ep;               % prior assumption: parameters same as control
pC      = NMD.M.pC;         % free up parameters as for initial inversio
pC.L    = pC.L / 1000;      % Switch off lead field variation
pC.J    = pC.J / 1000;      % Switch off population contribution variation

% Scale neuronal parameters by observed microscale behaviour
%==========================================================================
% Excitatory connectivity 

% Set up model space to explain NMDAR-Ab effect
%--------------------------------------------------------------------------
