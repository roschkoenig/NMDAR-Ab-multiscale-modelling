function C = ncl_preproc_loadfiles(F)
fs      = filesep; 
doplot  = 0; 
Fbs     = {}; % {[49,51],[99,101],[149,151]};

% Find EDF files
%--------------------------------------------------------------------------
clear C
C(1).name       = 'Control';
C(1).path       = [F.data fs 'EEG'];
C(1).edfs       = {'M1548954677.edf', 8;...
                   'M1548958277.edf', 8;...
                   'M1548961877.edf', 8;...
                   'M1541793130.edf', 6;...
                   'M1541753111.edf', 6;...
                   'M1541709443.edf', 6 };
C(1).edf_Fs     = 512; 
C(1).abfs       = {'18124005.abf'};
C(1).abf_Fs     = F.abf_Fs; 

C(2).name       = 'NMDAR-Ab positive';
C(2).path       = [F.data fs 'EEG'];
C(2).edfs       = {'M1548954677.edf', 5; ...
                   'M1548958277.edf', 5; ...
                   'M1548961877.edf', 5};
               
C(2).abfs       = {'18201000.abf'};
C(2).abf_Fs     = F.abf_Fs; 

% Load edf files 
%--------------------------------------------------------------------------
for c = 1:length(C)
    
C(c).cell   = abfload([F.data fs 'whole cell patch clamp recordings' fs C(c).abfs{1}]); 
C(c).cell   = C(c).cell(1:10*60*C(c).abf_Fs);   % Limit to the first ten minutes
C(c).dat    = []; 
disp(['Found ' num2str(length(C(c).edfs)) ' files for ' C(c).name]); 

for e = 1:length(C(c).edfs)
    hdr         = ft_read_header([C(c).path fs C(c).edfs{e}]); 
    C(c).edf_Fs = hdr.Fs; 
    
    if size(C(c).edfs,2) > 1    % Signal channels predefined
    %----------------------------------------------------------------------
        ci  = find(strcmp(hdr.label, num2str(C(c).edfs{e,2})));
        dat = ft_read_data([C(c).path fs C(c).edfs{e,1}]);
        id  = C(c).edfs{e,2}; 
        
    else    % Automatically find the signal channels 
    %----------------------------------------------------------------------
        dat = ft_read_data([C(c).path fs C(c).edfs{e}]);
        % Find channels with high continuous variability (these are the signal
        % channels)
        %--------------------------------------------------------------------------
        clear medabsdiff
        for dr = 1:size(dat,1)
            medabsdiff(dr) = median(abs(diff(dat(dr, :)))); 
        end
        id = find(medabsdiff > 10);

        % Plot selected channels for quick sanity check
        %--------------------------------------------------------------------------
        if doplot
        figure
        subplot(1+length(id), 1, 1)
        scatter(1:size(dat,1), medabsdiff)
        for i = 1:length(id)
            subplot(1+length(id), 1, i+1)
            plot(dat(id(i),1:10000))
            ylim([4000,6000])
        end
        end
    end
    
    if e > 1, minlength   = min(size(dat,2), size(C(c).dat,2));
              C(c).dat    = [C(c).dat(:,1:minlength); dat(id,1:minlength)]; 
    else,     C(c).dat    = dat(id,:); end
    
end
for f = 1:length(Fbs)
    C(c).dat = ft_preproc_bandstopfilter(C(c).dat, C(c).edf_Fs, Fbs{f});
end
end