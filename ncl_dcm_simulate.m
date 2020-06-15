fs  = filesep;
I   = load([F.outp fs 'DCM' fs 'Ictogenesis']); 
CNT = I.AllCM{1};
NMD = I.AllCM{2}; 
SZR = I.AllCM{3}; 

% Simulate baseline (sanity check
%--------------------------------------------------------------------------
% y = ncl_spm_csd_mtf(CNT.Ep, CNT.M, CNT.xU); 
% plot(log(abs(y{1}))); hold on
% plot(log(abs(CNT.xY.y{1})))

%% Microscale differences
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
EG(2).vec          = nmd_pars - cnt_pars;     

% Ictogenicity vector
%--------------------------------------------------------------------------
for i = 1:length(I.Sfx.Pnames)
    IG(1).id(i) = spm_fieldindices(CNT.Ep, I.Sfx.Pnames{i}); 
end
IG(1).vec          = zeros(size(cnt_pars)); 
IG(1).vec(IG.id)   = I.Sfx.vec; 
IG(2).vec          = szr_pars - cnt_pars; 


% Simulation
%==========================================================================
% Set up simulation settings
%--------------------------------------------------------------------------
e_effect = 2;   % 1 - 'reduced' or 2 - 'full'
i_effect = 2;   % 1 - 'reduced' or 2 - 'full'

steps  = 50; 
escale = linspace(0, 1, steps);   
iscale = linspace(0, 1, steps); 
frqhi  = intersect(find(CNT.M.Hz >= 0), find(CNT.M.Hz < 8));
frqlo  = intersect(find(CNT.M.Hz >= 8), find(CNT.M.Hz < 60)); 

clear sims ratio


% Run simulation
%--------------------------------------------------------------------------
for e = 1:steps
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
for r = 2:(size(sims,1)-1)
for c = 2:(size(sims,2)-1)
    % Average of all neighbours
    %----------------------------------------------------------------------
%     neighbours = {sims{r-1,c-1}, sims{r-1,c}, sims{r-1,c+1}, ...
%                   sims{r,c-1}, sims{r,c+1}, ...
%                   sims{r+1,c-1}, sims{r+1,c}, sims{r+1,c+1}};
    neighbours = {sims{r,c-1}, sims{r,c+1}};
    nsum = zeros(size(neighbours{1}));
    for n = 1:length(neighbours)
        nsum = nsum + abs(log(neighbours{n})); 
    end
    nsum = nsum / n;
    diff(r-1,c-1) = log(mean((log(abs(sims{r,c})) - nsum).^2)); 
end
end


subplot(2,2,1), imagesc(log(ratio)); axis square
subplot(2,2,2), imagesc(diff);  axis square

subplot(2,2,3), imagesc(paramspace);    axis square
subplot(2,2,4), 
    cols = cbrewer('seq', 'YlOrRd', size(diff,1));
    for d = 1:size(diff,1)  % epileptogenicity
        plot(exp(diff(d,:)), 'color', cols(d,:));  hold on
        axis square
        pause
    end


%% Local functions
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