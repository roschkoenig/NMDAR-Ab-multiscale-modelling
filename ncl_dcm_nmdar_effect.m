function ncl_dcm_nmdar_effect(C,F) 
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
ci        = find(strcmp({C.name}, 'Control'));
ni        = find(strcmp({C.name}, 'NMDAR-Ab positive')); 

% Set up model space
%--------------------------------------------------------------------------
mspace( 1).pops = {'deep', 'superficial'}; mspace( 1).pars = {'amplitude', 'time constant', 'variance'}; 
mspace( 2).pops = {'deep', 'superficial'}; mspace( 2).pars = {'amplitude', 'time constant'};            
mspace( 3).pops = {'deep', 'superficial'}; mspace( 3).pars = {'amplitude', 'variance'};                
mspace( 4).pops = {'deep', 'superficial'}; mspace( 4).pars = {'time constant', 'variance'};        
mspace( 5).pops = {'deep', 'superficial'}; mspace( 5).pars = {'amplitude'};                          
mspace( 6).pops = {'deep', 'superficial'}; mspace( 6).pars = {'time constant'};                         
mspace( 7).pops = {'deep', 'superficial'}; mspace( 7).pars = {'variance'};                             
mspace( 8).pops = {'deep'};                mspace( 8).pars = {'amplitude', 'time constant', 'variance'};              
mspace( 9).pops = {'deep'};                mspace( 9).pars = {'amplitude', 'time constant'};                            
mspace(10).pops = {'deep'};                mspace(10).pars = {'amplitude', 'variance'};                                 
mspace(11).pops = {'deep'};                mspace(11).pars = {'time constant', 'variance'};                             
mspace(12).pops = {'deep'};                mspace(12).pars = {'amplitude'};                                             
mspace(13).pops = {'deep'};                mspace(13).pars = {'time constant'};                                         
mspace(14).pops = {'deep'};                mspace(14).pars = {'variance'};                                              
mspace(15).pops = {'superficial'};         mspace(15).pars = {'amplitude', 'time constant', 'variance'};         
mspace(16).pops = {'superficial'};         mspace(16).pars = {'amplitude', 'time constant'};                    
mspace(17).pops = {'superficial'};         mspace(17).pars = {'amplitude', 'variance'};                         
mspace(18).pops = {'superficial'};         mspace(18).pars = {'time constant', 'variance'};                     
mspace(19).pops = {'superficial'};         mspace(19).pars = {'amplitude'};                                      
mspace(20).pops = {'superficial'};         mspace(20).pars = {'time constant'};                                  
mspace(21).pops = {'superficial'};         mspace(21).pars = {'variance'};                                      
mspace(22).pops = {};                      mspace(22).pars = {};                   

clear P
for m = 1:length(mspace)
this_pE = pE; 

% Excitatory connectivity 
%--------------------------------------------------------------------------
% Custom Mapping 
%--------------------------------------------------------------------------
% G(:,1)  ss -> ii (+ve rec )  4    Excitatory
% G(:,2)  dp -> ii (+ve rec )  2    Excitatory
% G(:,3)  ss -> sp (+ve rec )  4    Excitatory
% G(:,4)  sp -> ss (-ve rec )  4    Inhibitory
% G(:,5)  ii -> ss (-ve rec )  4    Inhibitory
% G(:,6)  ii -> dp (-ve rec )  2    Inhibitory
% G(:,7)  ss -> ss (-ve self)  4    Gain
% G(:,8)  ii -> ii (-ve self)  4    Gain
% G(:,9)  sp -> sp (-ve self)  4    Gain
% G(:,10) dp -> dp (-ve self)  1    Gain

if ~isempty(find(strcmp(mspace(m).pars, 'amplitude')))
    pid       = [];
    amp_ratio = log(C(ni).amplitude / C(ci).amplitude);
    if ~isempty(find(strcmp(mspace(m).pops, 'deep'))),        pid = [pid(:); 1,2];  end
    if ~isempty(find(strcmp(mspace(m).pops, 'superficial'))), pid = [pid(:); 3];    end
    this_pE.G(pid) = this_pE.G(pid) + amp_ratio; 
end

% Time constants
%--------------------------------------------------------------------------
% populations
% T(1)  ss   
% T(2)  sp  receives excitation
% T(3)  ii  receives excitation
% T(4)  dp

if ~isempty(find(strcmp(mspace(m).pars, 'time constant')))
    pid         = []; 
    tim_ratio   = log(C(ni).halflife / C(ci).halflife);
    if ~isempty(find(strcmp(mspace(m).pops, 'deep'))),        pid = [pid(:); 3];  end
    if ~isempty(find(strcmp(mspace(m).pops, 'superficial'))), pid = [pid(:); 2];    end
    this_pE.T(pid) = this_pE.T(pid) + tim_ratio; 
end

% Set up model space to explain NMDAR-Ab effect
%--------------------------------------------------------------------------
if ~isempty(find(strcmp(mspace(m).pars, 'variance')))
    sig_ratio   = log(C(ni).sigma / C(ci).sigma); 
    this_pE.S   = this_pE.S + sig_ratio; 
end

% Pack up and reinvert
%--------------------------------------------------------------------------
disp(['Inverting model ' num2str(m) ' of ' num2str(length(mspace))]); 
P{m}            = rmfield(NMD, 'name'); 
P{m}            = rmfield(P{m}, 'F'); 
P{m}.M.pE       = this_pE; 
P{m}.nmda_mod   = mspace(m);
P{m}            = ncl_spm_dcm_csd(P{m});

end

save([F.outp fs 'DCM' fs 'Model_comparison_NMDAR.mat'], 'P', 'mspace')