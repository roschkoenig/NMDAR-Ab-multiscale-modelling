fs  = filesep ;
F   = ncl_housekeeping;                 % Specify filenames and path  
C   = ncl_preproc_loadfiles(F);         % Loads original datafiles
C   = ncl_preproc_epspseg(C);           % Segments out EPSCs
C   = ncl_preproc_featextract(C);       % Calculates some statistics based on EPSCs
F   = ncl_preproc_segment(C,F);         % Segment out relative data segments

%% DCM part
ncl_dcm_specify(F); 
