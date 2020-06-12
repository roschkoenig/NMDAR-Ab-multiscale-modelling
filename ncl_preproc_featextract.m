function C = ncl_preproc_featextract(C)

for c = 1:length(C)
    Fs = C(c).abf_Fs; 
    
    % Calculate peak amplitude
    %----------------------------------------------------------------------
    mDec = mean(C(c).segs); 
    [peak, peaksample] = min(mDec); 
    C(c).amplitude     = abs(peak); 
    
    % Calculate half life (proportional to time constant)
    %----------------------------------------------------------------------    
    halfpksmpl         = find(mDec < peak/2); 
    C(c).halflife      = (halfpksmpl(end) - peaksample) / Fs;
    
    % Calculate variability
    %----------------------------------------------------------------------
    peaks       = min(C(c).segs');
    [cnts bins] = histcounts(log(-peaks), 100); 
    freqs       = cnts / size(C(c).segs,1); 
    cfrqs       = cumsum(freqs);
    inflex      = find(cfrqs >= 0.5);
    inflex      = inflex(1);
    bins        = bins-(bins(inflex)+bins(inflex+1))/2; 
    
    dfrqs       = diff(cfrqs);
    sigma       = dfrqs(inflex); 
    C(c).sigma  = sigma;
    plot(bins(1:end-1), cfrqs), hold on
    plot(bins(2:end-1), dfrqs); 
end