function C = ncl_preproc_epspseg(C)

for c = 1:length(C)
Fs   = C(c).abf_Fs; 
filt = ft_preproc_bandpassfilter(C(c).cell', Fs, [0.3, 2000]);

% Find peaks 
[val loc, w, p] = findpeaks(-filt, 'MinPeakProminence', 20, 'MinPeakHeight', 0);
loc             = loc(val < 200); % This will exclude action potentials 


% Segment out EPSCs
%--------------------------------------------------------------------------
win       = [-0.03, 0.06] * Fs; 
bl        = [0, 0.02] * Fs + 1; 
segs = [];
tims = []; 
for l = 1:length(loc)
    if (loc(l)+win(1) > 0) && (loc(l)+win(2) < length(filt))    % window fully within dataset
    if range(filt(loc(l)+win(1):loc(l)+win(2))) < 250           % Exclude high amplitude artefacts
        seg = filt(loc(l)+win(1):loc(l)+win(2));  
        seg = seg - mean(seg(bl(1):bl(end))); 
        if isempty(find((seg) > 20))
            segs = [segs; seg];
            tims = [tims, loc(l)/Fs/60];  
        end
    end % if 
    end % if 
end % loop through individually detected EPSPs
C(c).segs = segs; 
C(c).tims = tims; 
end
