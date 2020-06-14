function ncl_dcm_bmc(F)
fs  = filesep;
load([F.outp fs 'DCM' fs 'Model_comparison_NMDAR']); 
CNT = load([F.outp fs 'DCM' fs 'DCM_Control baseline.mat']); 
CNT = CNT.DCM; 

% Perform individual model Bayesian model comparison
%--------------------------------------------------------------------------
clear FE 
for p = 1:length(P)
    FE(p)       = P{p}.F;
    names{p}    = ''; 
    for k = 1:length(mspace(p).pars)
        names{p} = [names{p}, mspace(p).pars{k} ' ']; 
    end
end

figure(1)
subplot(2,1,1),     bar(FE - min(FE)); xlabel('Models'); 
                    set(gca, 'Xtick', 1:length(FE)); ylabel('Free energy difference'); 
subplot(2,1,2),     bar(spm_softmax(FE')); xlabel('Models'); ylabel('Posterior probability'); 
                    set(gca, 'Xtick', 1:length(FE))
                    set(gca, 'Xticklabel', names)
                    xtickangle(90)

% Family-wise comparison
%--------------------------------------------------------------------------
fams = [1 1 1 1 1 1 1 2 2 2 2 2 2 2 3 3 3 3 3 3 3 4]; % Family by effect location
clear famFE
for k = unique(fams)
    id = find(fams == k); 
    famFE(k) = mean(FE(id)); 
end

% Plot family wise comparison
%--------------------------------------------------------------------------
figure(2)
subplot(2,2,[1,3]); 
    bar(famFE - min(famFE))
    [val winner] = max(famFE); 
    title('Bayesian model comparison - Family-wise');
    xlabel('Model families');   ylabel('Free energy difference');

% Bayesian model averaging
%--------------------------------------------------------------------------
id = find(fams == winner);   
AVG = spm_dcm_bma(P(id)); 

cnt = [CNT.Ep.G(1:3), CNT.Ep.T(2:3), CNT.Ep.S];
prd = [P{1}.M.pE.G(1:3), P{1}.M.pE.T(2:3), P{1}.M.pE.S]; 
emp = [AVG.Ep.G(1:3), AVG.Ep.T(2:3), AVG.Ep.S]; 
Cp  = [AVG.Cp.G(1:3), AVG.Cp.T(2:3), AVG.Cp.S]; 

subplot(2,2,2), bar(prd - cnt); ylim([-1.5 1.5]); title('Prediction from microscale');
subplot(2,2,4), bar(emp - cnt); ylim([-1.5 1.5]); title('Estimation with microscale-informed priors');
set(gca, 'Xtick', [1:6], 'XtickLabels', {'G_1', 'G_2', 'G_3', 'T_2', 'T_3', 'S'})

save([F.outp fs 'DCM' fs 'Model_comparison_NMDAR'], 'P', 'mspace', 'AVG'); 
