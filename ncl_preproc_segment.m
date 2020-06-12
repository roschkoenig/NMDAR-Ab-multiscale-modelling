function F = ncl_preproc_segment(C,F)
fs = filesep; 
%% Plot z-scored segments
%--------------------------------------------------------------------------
z_cutoff    = 5.5; 
timewin     = 45; 
bl_win      = 30 * C(1).edf_Fs; 
sz_win      = 30 * C(1).edf_Fs; 

% Control means data
%--------------------------------------------------------------------------
cmean = mean(C(1).dat(:)); 
cstd  = std(C(1).dat(:)); 
Sz = {};
Bl = {}; 

for c = 1:length(C)
Fs = C(c).edf_Fs; 
disp(['Processing ' C(c).name]);

for d = 1:size(C(c).dat,1)
    disp(['Dataset ' num2str(d) ' of ' num2str(size(C(c).dat,1))])
    
    % Z-score based on control data and find excessive amplitude regions
    %----------------------------------------------------------------------
    zdat = (C(c).dat(d,:) - cmean)/cstd; 
    excess  = find(abs(zdat) > z_cutoff); 
    classfd = ones(length(C(c).dat),1); 
    for e = 1:length(excess)
        srt = max(1, excess(e)-timewin*Fs);
        stp = min(length(C(c).dat), excess(e)+timewin*Fs); 
        classfd(srt:stp) = zeros(length(srt:stp),1); 
    end    
    
    % Plot colour coded by their classification
    %----------------------------------------------------------------------
    yes = find(classfd);
    no  = find(~classfd);
    t   = [1:length(C(c).dat)]/(60*Fs);
    
    subplot(2,1,c)
        plot(t(no), zdat(no) + d*10, 'color', [.5 .5 .5]), hold on
        plot(t(yes), zdat(yes) + d*10, 'color', 'k'), hold on
        ylim([0 (d+1)*10])
    
    % Segment out baseline and seizure segments
    %---------------------------------------------------------------------
    boundaries = [0, find(diff(classfd))', length(classfd)]; 
    for b = 1:length(boundaries)-1
        if classfd(boundaries(b)+1) == 1  
            newid       = size(Bl,1)+1;
            Bl{newid,1} = C(c).dat(d, boundaries(b)+1:boundaries(b+1)); 
            Bl{newid,2} = c; 
            Bl{newid,3} = C(c).name; 
        elseif strcmp(C(c).name, 'NMDAR-Ab positive')
            newid       = size(Sz,1) + 1;
            Sz{newid,1} = C(c).dat(d, boundaries(b)+1:boundaries(b+1));
            Sz{newid,2} = c;
            Sz{newid,3} = C(c).name; 
        end
    end
end
end

% Cut out segments of data that can be used for further analysis
%--------------------------------------------------------------------------
segs    = []; 
conds   = {}; 

% Take all chunks from the Baseline 
%--------------------------------------------------------------------------
for b = 1:length(Bl)    
    Nwin    = floor(length(Bl{b,1}) / bl_win);
    for n = 1:Nwin
        segs    = [segs; Bl{b,1}((n-1)*bl_win+1 : n*bl_win)];
        conds   = {conds{:}, [Bl{b,3} ' baseline']};
    end
end

% Take sequence of windows from seizure (preserve order)
%--------------------------------------------------------------------------
for s = 1:length(Sz)
    Nwin    = floor(length(Sz{s,1}) / sz_win);
    for n = 1:Nwin
        segs    = [segs; Sz{s,1}((n-1)*sz_win+1 : n*sz_win)]; 
        conds   = {conds{:}, ['NMDAR-Ab Seizure Window ' num2str(n, '%03d')]}; 
    end
end

% Manual correction
%--------------------------------------------------------------------------
if size(segs,1) == 863
    disp('I have removed previously selected artefact segments'); 
    ind = ones(size(segs,1), 1); 
    ind([297, 362, 364, 365, 476, 477, 533]) = 0; 
    segs    = segs(find(ind),:);
    conds   = conds(find(ind));
    
else
    figure
    disp('The number of segments is not what I expected - please check manually');
    for s = 1:10:size(segs,1)
        clf
        for k = 1:10
        if s+k <= size(segs,1)
            plot(segs(s+k,:)+k*400, 'k', 'linewidth', 0.5); hold on
        end
        end
        title(s)
        pause
    end
end


% Plot segments as fourier spectra
%--------------------------------------------------------------------------
close all
cnt     = find(~cellfun(@isempty, strfind(conds, 'Control')));
nmd     = find(~cellfun(@isempty, strfind(conds, 'NMDAR-Ab positive')));
szr     = find(~cellfun(@isempty, strfind(conds, 'Seizure'))); 
szrt    = unique(conds(szr)); 

for t = 1:length(szrt)
    szrt{t} = find(~cellfun(@isempty, strfind(conds, szrt{t}))); 
end

nyquist     = floor(C(1).edf_Fs/2); 
Nffsmpls    = floor(size(segs,2)/2); 
frqs        = linspace(0, nyquist, Nffsmpls); 

fcnt        = log(mean(abs(fft(segs(cnt,:)')),2)); 
fnmd        = log(mean(abs(fft(segs(nmd,:)')),2)); 
fszr        = log(mean(abs(fft(segs(szr,:)')),2)); 
fszrt       = {};
for t = 1:length(szrt)
    fszrt{t}       = log(mean(abs(fft(segs(szrt{t},:)')),2)); 
end

subplot(2,1,1)
    plot(frqs(2:end), fcnt(2:Nffsmpls),'k'), hold on
    plot(frqs(2:end), fnmd(2:Nffsmpls),'b')
    plot(frqs(2:end), fszr(2:Nffsmpls),'r')
    xlabel('Frequency'); 
    ylabel('Log-power');
    title('Average power spectrum by condition'); 
    xlim([0 50])

subplot(2,1,2)
    reds = cbrewer('seq', 'Reds', length(fszrt)); 
    for t = 1:length(fszrt)
        fullfft = fszrt{t}(2:Nffsmpls);
        fullfrq = frqs(2:end);
        stp     = .5;
        fr      = [];
        ft      = [];
        k       = 0;
        for b = 0:stp:fullfrq(end)-1
            k  = k + 1;
            id = intersect(find(fullfrq >= b), find(fullfrq < b+stp));
            fr(k)   = mean(fullfrq(id)); 
            ft(k)   = mean(fullfft(id));
        end
       
        plot(fr, ft, 'color', reds(t,:)), hold on
    end
    colormap(reds), caxis([0, length(fszrt)*0.5]), cb = colorbar();
    xlim([0 50]); 
    xlabel('Frequency');
    ylabel('Log-power'); 
    title('Transition into seizure'); 
    ylabel(cb, 'Minutes from onset'); 

% Save everything in an SPM file
%--------------------------------------------------------------------------
for s = 1:size(segs,1)
    ftdata.trial{s} = segs(s,:); 
    ftdata.time{s}  = linspace(0, (size(segs,2)-1)/C(1).edf_Fs, size(segs,2));
end
ftdata.label = 'LFP';
try mkdir([F.outp fs 'DCM']);  end
D = spm_eeg_ft2spm(ftdata, [F.outp fs 'DCM' fs 'LFP_seg']); 
for c = 1:length(conds)  
    if strfind(conds{c}, 'Seizure'), cond = 'NMDAR-Ab positive seizure';
    else,                            cond = conds{c}; end 
    D = conditions(D, c, cond); 
end
save(D)

F.spm_file = [F.outp fs 'DCM' fs 'LFP_seg']; 
