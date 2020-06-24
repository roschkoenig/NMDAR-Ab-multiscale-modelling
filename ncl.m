fs  = filesep ;
F   = ncl_housekeeping;                 % Specify filenames and path  
C   = ncl_preproc_loadfiles(F);         % Loads original datafiles
C   = ncl_preproc_epspseg(C);           % Segments out EPSCs
C   = ncl_preproc_featextract(C);       % Calculates some statistics based on EPSCs
F   = ncl_preproc_segment(C,F);         % Segment out relevant data segments

%% DCM part
% ncl_dcm_specify(F);                     % Specify DCMs and estimate CSDs 
% ncl_dcm_invert_all(F);                  % Initial inversion of all DCMs 
% ncl_dcm_nmdar_effect(C,F);              % Evaluate microscale explanation of NMDAR-Ab effect
% ncl_dcm_bmc(F); 
ncl_dcm_seizure(F); 
% ncl_dcm_simulate(F); 