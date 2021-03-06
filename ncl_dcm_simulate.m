function ncl_dcm_simulate(F)
fs  = filesep;
I   = load([F.outp fs 'DCM' fs 'Ictogenesis']); 
CNT = I.AllCM{1};
NMD = I.AllCM{2}; 
SZR = I.AllCM{3}; 

% % Simulate baseline (sanity check
% % --------------------------------------------------------------------------
% y = ncl_spm_csd_mtf(CNT.Ep, CNT.M, CNT.xU); 
% plot(log(abs(y{1}))); hold on
% plot(log(abs(CNT.xY.y{1})))

% Microscale differences
%--------------------------------------------------------------------------
clear EG IG
cnt_pars        = spm_vec(CNT.Ep); 
nmd_pars        = spm_vec(NMD.Ep); 
szr_pars        = spm_vec(SZR.Ep); 

% Epileptogenicity vector (1) reduced
%--------------------------------------------------------------------------
EG(1).id           = [spm_fieldindices(CNT.Ep, 'G(1:3)'); 
                      spm_fieldindices(CNT.Ep, 'T(2:3)'); 
                      spm_fieldindices(CNT.Ep, 'S')]; 
EG(1).vec          = zeros(size(cnt_pars)); 
EG(1).vec(EG.id)   = nmd_pars(EG.id) - cnt_pars(EG.id);

% Epiletogenicity vector (2) full
%--------------------------------------------------------------------------
EG(2).vec          =  nmd_pars - cnt_pars; %- EG(1).vec;

% Ictogenicity vector
%--------------------------------------------------------------------------
% for i = 1:length(I.Sfx.Pnames)
%     IG(1).id(i) = spm_fieldindices(CNT.Ep, I.Sfx.Pnames{i}); 
% end
% IG(1).vec          = zeros(size(cnt_pars)); 
% IG(1).vec(IG.id)   = I.Sfx.vec; 
IG(2).vec          = szr_pars - nmd_pars; 


% Simulation
%==========================================================================
% Set up simulation settings
%--------------------------------------------------------------------------
e_effect = 2;   % 1 - 'reduced' or 2 - 'full'
i_effect = 2;   % 1 - 'reduced' or 2 - 'full'

steps  = 20; 
escale = linspace(0, 1, steps);   
iscale = linspace(0, 1, steps); 

% Run simulation
%--------------------------------------------------------------------------
clear sims ratio
for e = 1:steps
    disp(['Simulating step ' num2str(e) ' of ' num2str(steps)]); 
for i = 1:steps 
    simpars         = spm_unvec(cnt_pars + escale(e) * EG(e_effect).vec + iscale(i) * IG(i_effect).vec, CNT.Ep); 
    y               = ncl_spm_csd_mtf(simpars, CNT.M, CNT.xU); 
    sims{e,i}       = y{1}; 
end
end

%% Plot results
%--------------------------------------------------------------------------
% Calculate spectral patterns
%--------------------------------------------------------------------------
frqhi        = intersect(find(CNT.M.Hz >= 0), find(CNT.M.Hz < 8));
frqlo        = intersect(find(CNT.M.Hz >= 20), find(CNT.M.Hz < 60)); 
paramspace   = zeros(size(sims)); 
hi           = zeros(size(sims)); 
lo           = zeros(size(sims)); 

for r = 1:size(sims,1)     % epileptogenicity
for c = 1:size(sims,2)     % ictogenicity
    paramspace(r,c) = classification(sims{r,c}, {CNT, NMD, SZR});
    hi(r,c)         = mean(log(abs(sims{r,c}(frqhi))));
    lo(r,c)         = mean(log(abs(sims{r,c}(frqlo))));
end
end

ratio = hi ./ lo; 

% Calculate difference function
%--------------------------------------------------------------------------
clear diff
for r = 2:(size(sims,1)-1)
for c = 2:(size(sims,2)-1)
    % Average of all neighbours
    %----------------------------------------------------------------------
    neighbours = {sims{r,c-1}, sims{r,c+1}};
    nsum = zeros(size(neighbours{1}));
    for n = 1:length(neighbours)
        nsum = nsum + abs(log(neighbours{n})); 
    end
    nsum = nsum / n;
    diff(r-1,c-1) = log(mean((log(abs(sims{r,c})) - nsum).^2)); 
end
end

fn = 0; 
close all
figure(fn+1),   fn = fn+1;  
    set(gcf, 'color', 'w', 'position', [200,800,400,400])
    cols = cbrewer('qual', 'Dark2', 3); 
    plot(log(abs(CNT.xY.y{1})), '--', 'color', cols(3,:), 'linewidth', 1.5), hold on
    plot(log(abs(CNT.Hc{1})), 'color', cols(3,:), 'linewidth', 1.5),
    plot(log(abs(NMD.xY.y{1})), '--', 'color', cols(3,:), 'linewidth', 1.5), hold on
    plot(log(abs(NMD.xY.y{1})), 'color', cols(1,:), 'linewidth', 1.5), 
    plot(log(abs(SZR.xY.y{1})), '--', 'color', cols(3,:), 'linewidth', 1.5), hold on
    plot(log(abs(SZR.xY.y{1})), 'color', cols(2,:), 'linewidth', 1.5); 

figure(fn+1),   fn = fn+1;  
    set(gcf, 'color', 'w', 'position', [600,800,400,400])
    colormap(flip(cbrewer('div', 'Spectral', 100))); 
    imagesc((ratio)); axis square
    title('Spectral features across simualted parameter space'); 
    set(gca, 'Ytick', [1,size(ratio,1)], 'YTickLabel', {'control', 'NMDAR-Ab +ive'});
    set(gca, 'Xtick', [1,size(ratio,2)], 'XTickLabel', {'interictal', 'ictal'}); 
    ytickangle(90); 
    cb = colorbar();
    ylabel(cb, 'High/Low frequency amplitude ratio');
    
figure(fn+1), fn = fn+1;  
    set(gcf, 'color', 'w', 'position', [1000,800,400,400])
    cm  = cbrewer('qual', 'Set2', 3);
    colormap(cm([3 1 2],:)); 
    imagesc(paramspace);    axis square
    title('State classification'); 
    set(gca, 'Ytick', [1,size(ratio,1)], 'YTickLabel', {'control', 'NMDAR-Ab +ive'});
    set(gca, 'Xtick', [1,size(ratio,2)], 'XTickLabel', {'interictal', 'ictal'}); 
    ytickangle(90)
    cb = colorbar;
    set(cb, 'ytick', [1 2 3], 'yticklabel', {'CNT', 'NMD', 'SZR'}); 
    
figure(fn+1), fn = fn+1;  
    set(gcf, 'color', 'w', 'position', [1000,300,400,400])
    imagesc(diff);  axis square
    colormap(flip(cbrewer('div', 'RdGy', 100))); 
    cb = colorbar;
    ylabel(cb, 'Sensitivity to parameter change'); 
    set(gca, 'Ytick', [1,size(ratio,1)-2], 'YTickLabel', {'control', 'NMDAR-Ab +ive'});
    set(gca, 'Xtick', [1,size(ratio,2)-2], 'XTickLabel', {'interictal', 'ictal'});
    ytickangle(90); 
    
figure(fn+1), fn = fn+1;  
    set(gcf, 'color', 'w', 'position', [200,300,800,400]); 
    Hz  = CNT.M.Hz; 
    subplot(1,2,1)
        cols = cbrewer('qual', 'Dark2', 3); 
        gris = (cbrewer('seq', 'YlOrRd', size(sims,2))); 
        cnt  = log(abs(CNT.xY.y{1})); 
        for s =1:5:size(sims,2), plot(Hz,log(abs(sims{1,s}))-cnt, 'color', gris(s,:)); hold on; end
        plot(Hz,log(abs(CNT.xY.y{1}))-cnt, '--', 'color', cols(3,:), 'linewidth', 2); hold on
        plot(Hz,log(abs(SZR.xY.y{1}))-cnt, '--','color', cols(2,:), 'linewidth', 2); 
        title('Sensitivity of normal circuit to ictogenic parameter changes'); 
        ylim([-2.5 2.5]); 
        
    subplot(1,2,2)
        cols = cbrewer('qual', 'Dark2', 3); 
        gris = (cbrewer('seq', 'YlOrRd', size(sims,2))); 
        nmd  = log(abs(NMD.xY.y{1})); 
        for s =1:3:size(sims,2), plot(Hz,log(abs(sims{end,s}))-nmd, 'color', gris(s,:)); hold on; end
        plot(Hz,log(abs(NMD.xY.y{1}))-nmd, '--', 'color', cols(1,:), 'linewidth', 2); hold on
        plot(Hz,log(abs(SZR.xY.y{1}))-nmd, '--', 'color', cols(2,:), 'linewidth', 2); 
        title('Sensitivity of NMDAR-Ab circuit to ictogenic parameter changes'); 
        ylim([-2.5, 2.5])
        
figure(fn+1), fn = fn+1;  
    set(gcf, 'color', 'w', 'position', [1400,300,400,800]); 
    plthis = [ceil(0.1*size(sims,1)), ceil(0.5*size(sims,1)), ceil(0.9*size(sims,1))]; 

    for p = 1:3
    plid = plthis(p);
    separat = 0; 
    subplot(3,1,p)
    for s = 1:size(sims,2)
        colid = [3 1 2];
        plot(log(abs(sims{plid,s})) - log(abs(sims{1,1})) +s*separat, 'color', cols(colid(paramspace(plid,s)),:))
        ylim([-2 2.5]); 
        hold on
    end
    end
end

% Local functions
%-=========================================================================
function cl = classification(y, dcms)
    for d = 1:length(dcms)
        y_sc = abs(log(y));              % - mean(abs(log(y{1})));
        t_sc = abs(log(dcms{d}.xY.y{1}));   % - mean(abs(log(dcms{d}.xY.y{1})));
        mse(d) = mean( [y_sc - t_sc].^2 );
    end
    [val cl] = min(mse); 
end


% sim_pars            = cnt_pars; 
% sim_pars(egen_id)   = cnt_pars(egen_id) + dif_pars(egen_id); 
% sim_pars            = spm_unvec(sim_pars, CNT.Ep); 



% All Empirical differences 

% %%
% CNT = load('ICM_Control baseline.mat'); 
% NMD = load('ICM_NMDAR-Ab positive baseline.mat'); 
% SZR = load('ICM_NMDAR-Ab positive seizure.mat'); 
% 
% subplot(2,1,1) 
%     plot(log(abs(CNT.DCM.xY.y{1})), 'k'); hold on
%     plot(log(abs(NMD.DCM.xY.y{1})), 'b'); hold on
%     plot(log(abs(SZR.DCM.xY.y{1})), 'r'); hold on
% 
% 
% subplot(2,1,2)
% 
% for p = 1:length(P)
%     plot(log(abs(P{p}.xY.y{1}))); hold on
% end