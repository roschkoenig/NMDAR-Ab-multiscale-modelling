function F = ncl_housekeeping()
    fs = filesep; 
    F.base = '/Users/roschkoenig/Dropbox/Research/1906 NMDAR Ab Cells';
    F.code = [F.base fs '01_Code'];
    F.data = [F.base fs '02_Data'];
    F.outp = [F.base fs '03_Output'];
    F.addpaths = {
        '/Users/roschkoenig/Dropbox/Research/0002 Tools/Matlab/spm',
        '/Users/roschkoenig/Dropbox/Research/0002 Tools/Matlab/abfload',
        '/Users/roschkoenig/Dropbox/Research/0002 Tools/Matlab/cbrewer',
        F.code
    };

    F.abf_Fs = 10000; 

    for a = 1:length(F.addpaths), addpath(F.addpaths{a}); end
    
    spm('defaults', 'eeg'); 