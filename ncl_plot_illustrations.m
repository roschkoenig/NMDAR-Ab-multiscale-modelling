
for s = 1:2 
    subplot(2,1,s), plot(C(s).segs(1:floor(size(C(s).segs,1)/100):size(C(s).segs,1),:)'); ylim([-200 50]); hold on
                    plot(mean(C(s).segs), 'linewidth', 2); 
end

%%
for c = 1:2
    peak_seg    = 276:350; 
    peaks       = abs(min(C(c).segs(:,peak_seg)'));
    [cnts bins] = histcounts(log(peaks), 100); 
    freqs       = cnts / size(C(c).segs,1); 
    cfrqs       = cumsum(freqs);

    lrange = linspace(-6, 6, 1001);
    srange = exp(linspace(-4, 2, 1000)); 
    F_m     = @(beta,xdata) 1 ./ (1 + exp(-(xdata-beta(1))./beta(2)));
    LSE     = zeros(length(lrange), length(srange)); 

    for l = 1:length(lrange)
    for s = 1:length(srange)
        beta_0      = [lrange(l), srange(s)];
        x           = bins(1:end-1); 
        y           = cfrqs; 
        LSE(l,s)    = mean((y - F_m(beta_0, x)).^2);
    end
    end
    subplot(2,2,c), 
        imagesc(LSE); hold on
        set(gca, 'Ydir', 'normal'); 

        [l,s] = find(LSE == min(LSE(:)));
        scatter(s, l, 'r', 'filled'); 

    subplot(2,2,c+2)
        scatter(x, cfrqs, 'k', 'filled'); hold on
        plot(linspace(0,6,100),F_m([lrange(l), srange(s)], linspace(0,6,100)), 'b'); 
        ylim([0 1]);
        xlim([0 6]); 
         
end


%% Plot parameterised sigmoid
%--------------------------------------------------------------------------
F_h = @(params, t) params(1) * params(2) * t.* exp(-params(2)*t);
params = [1 1/10];
t   = 1:0.2:100; 
plot(F_h(params, t))

%% Plot parameters of model fits
%--------------------------------------------------------------------------
I   = load([F.outp fs 'DCM' fs 'Ictogenesis']); 
CNT = I.AllCM{1};
NMD = I.AllCM{2}; 
SZR = I.AllCM{3}; 

ids = [spm_fieldindices(CNT.Ep, 'G'); 
       spm_fieldindices(CNT.Ep, 'T'); 
       spm_fieldindices(CNT.Ep, 'S')]; 

CNT.vec = spm_vec(CNT.Ep); 
CNT.vec = CNT.vec(ids);

NMD.vec = spm_vec(NMD.Ep); 
NMD.vec = NMD.vec(ids); 
   