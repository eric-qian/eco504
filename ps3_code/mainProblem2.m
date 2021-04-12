% mainProblem2.m  Main file for problem 2
% TODO: SIMULATION OF LOG CONSUMPTION TENDS TO INFINITY
%% Preliminaries
clear
clc

addpath('functions/')

% Structure giving specifications
Spec = struct;

figPath = 'figures/';
figSize = [6,3];

% Simulation settings
nSim     = 11000;
nDiscard = 1000;


%% Baseline specification
j=1;
Spec(j).lab = 'Baseline';

% Parameters
Spec(j).phi     = 0;  % Borrowing limit
Spec(j).gamma   = 2;  % CRRA parameter
Spec(j).beta    = 0.95; 
Spec(j).r       = 0.02;
Spec(j).Na      = 1000;  
Spec(j).epsConv = 5e-7;

% Income process parameters
% w_t = wBar + rho w_t-1 + eps_t, eps_t~N(0, vareps)
Spec(j).vareps  = 0.06;                 % Variance of epsilon of w process
Spec(j).m       = 4;                    % grid is 0 +/- m*sigma^2_w
Spec(j).rho     = 0.9;                  % AR coefficient of w process
Spec(j).Nw      = 5;                    % Grid size

%% Alternative specifications

% Changing gamma, gamma = 1
j             = 2;
Spec(j)       = Spec(1);
Spec(j).lab   = '\gamma=1';
Spec(j).gamma = 1;

% j             = 7;
% Spec(j)       = Spec(1);
% Spec(j).lab   = '\gamma=1';
% Spec(j).gamma = 1;


% Changing gamma, gamma=5
j             = 3;
Spec(j)       = Spec(1);
Spec(j).lab   = '\gamma=5';
Spec(j).gamma = 5;

% Changing vareps, vareps = .01
j              = 4;
Spec(j)        = Spec(1);
Spec(j).lab    = '\sigma^2_\varepsilon=0.01';
Spec(j).vareps = 0.01;

% Changing vareps, vareps = .12
j              = 5;
Spec(j)        = Spec(1);
Spec(j).lab    = '\sigma^2_\varepsilon=0.12';
Spec(j).vareps = 0.12;


% Changing natural borrowing limit
j              = 6;
Spec(j)        = Spec(1);
Spec(j).lab    = '\phi=-y_{min}/r';
Spec(j).phi    = 'NBL';


%% Estimate models

% Compare solution methods
compare = 0;

if compare == 1
    j=1;
    
    Spec_j       = Spec(j);
    Spec_j.Na    = 600;  % Make the grid smaller so that the problem is tractable
    Spec_j.gamma = 1;

    
    SpecEG  = Spec_j;
    SpecVFI = Spec_j;
    
    SpecEG.Res   = solveConsSaving(SpecEG);
    SpecVFI.Res  = solveConsSaving(SpecVFI, 'VFI');
    
    aPolEG  = SpecEG.Res.aPol;
    aPolVFI = SpecVFI.Res.aPol;
    aGrid   = SpecEG.Res.aGrid;
    
    cPolEG   = SpecEG.Res.cPol;
    cPolVFI = SpecVFI.Res.cPol;
   
    
    close all
    subplot(2,1,1)
    plot(aGrid, cPolEG(:, 3))
    hold on
    plot(aGrid, cPolVFI(:, 3))
    title('c policy')
    
    subplot(2,1,2)
    plot(aGrid, aPolEG(:, 3))
    hold on
    plot(aGrid, aPolVFI(:, 3))
    title('a policy')
    legend({'EG', 'VFI'})
    
    rng(1)
    ResEG.Sim  = simConsSaving(SpecEG, nSim, nDiscard);
    rng(1)
    ResVFI.Sim = simConsSaving(SpecVFI, nSim, nDiscard);

    
    
    close all
    subplot(2,1,1)
    plot(ResEG.Sim.aPath)
    title('a path')
    hold on
    plot(ResVFI.Sim.aPath)
    formatFig(figSize)
%    xlim([0, 100])
    box on;grid on;
        subplot(2,1,2)
    plot(ResEG.Sim.cPath)
    hold on
    plot(ResVFI.Sim.cPath)
    title('c path')
    formatFig(figSize)
%    xlim([0, 100])

    
    legend({'EG', 'VFI'})
    saveas(gcf, [figPath 'VFI_EG_compare.png'])        
    
end

for j = 1:length(Spec)
    Spec_j      = Spec(j);
    Spec(j).Res = solveConsSaving(Spec_j, 'VFI');
end

%% Simulate model


for jSpec = 1:length(Spec)
    Spec_j      = Spec(jSpec);
    rng(1)
    Spec(jSpec).Sim = simConsSaving(Spec_j, nSim, nDiscard);
end


%% Look at .0001 

Spec1 = Spec(1);
Spec2 = Spec(2);


aPol1 = Spec1.Res.aPol;
aPol2 = Spec2.Res.aPol;

diff = (aPol2-aPol1) ./ aPol1 * 100;
ksdensity(diff(:))
close all
yyaxis left
plot(Spec2.Sim.aPath)
yyaxis right
plot(Spec2.Sim.yPath)
xlim([0, 250])



%% Part b

plotIdx  = [2, 1, 3];

cPlot    = []; 
aPlot    = [];

for j =plotIdx
    cPlot = [cPlot Spec(j).Sim.cPath];
    aPlot = [aPlot Spec(j).Sim.aPath];
    
end

sigmaVec = std(cPlot, 0, 1, 'omitnan');
plotLabs = {['\gamma=1 (\sigma=' num2str(sigmaVec(1), 3) ')'], ...
            ['\gamma=3 (\sigma=' num2str(sigmaVec(2), 3) ')'],...
            ['\gamma=5 (\sigma='  num2str(sigmaVec(3), 3) ')']};

close all
plot(cPlot())
xlim([0, nSim-nDiscard])
legend(plotLabs, 'Location', 'southoutside', 'Orientation', 'horizontal')
formatFig(figSize)
saveas(gcf, [figPath 'p2_partb.png'])

%% Part c
plotIdx  = [4, 1, 5];
plotLabs = {'\sigma_\epsilon^2=0.01', '\sigma_\epsilon^2=0.06', '\sigma_\epsilon^2=0.12'};

cPlot    = []; 
yPlot    = []; 

for j =plotIdx
    cPlot = [cPlot Spec(j).Sim.cPath];
    yPlot = [yPlot Spec(j).Sim.yPath(:)];
    
end

close all
sPlot = ones(size(cPlot))- cPlot./yPlot;
plot(sPlot)
xlim([0, nSim-nDiscard])
legend(plotLabs, 'Location', 'southoutside', 'Orientation', 'horizontal')
saveas(gcf, [figPath 'p2_partc.png'])


%% Part d
plotIdx  = [1, 6];

cPlot    = []; 
yPlot    = []; 
aPlot    = [];

for j =plotIdx
    cPlot = [cPlot Spec(j).Sim.cPath];
    aPlot = [aPlot Spec(j).Sim.aPath];    
    yPlot = [yPlot Spec(j).Sim.yPath'];
end

muVec = mean(cPlot, 'omitnan');
plotLabs = {['Baseline (\mu=' num2str(muVec(1), 3) ')'], ...
    ['NBL (\mu=' num2str(muVec(2),3) ')' ]};

close all
plot(cPlot)
xlim([0, nSim-nDiscard])
legend(plotLabs, 'Location', 'southoutside', 'Orientation', 'horizontal')
saveas(gcf, [figPath 'p2_partd.png'])


temp1=Spec(plotIdx(1)).Res.cPol
temp2=Spec(plotIdx(2)).Res.cPol

close all
plot(Spec(1).Res.aGrid, temp1(:,3,:))
hold on
plot(Spec(6).Res.aGrid, temp2(:,3,:))
legend({'BL', 'NBL'})


subplot(2,1,1)
plot(aPlot(:, 2) - aPlot(:, 1))
xlim([0, nSim])
subplot(2,1,2)
plot(cPlot(:, 2) - cPlot(:, 1))
xlim([0, nSim])

% Check budget constraint
j=1;
plot(Spec(j).Sim.cPath(1:end-1) + Spec(j).Sim.aPath(2:end) - (1+Spec(j).r)*Spec(j).Sim.aPath(1:end-1) - Spec(j).Sim.yPath(1:end-1)')
hold on
j=6
plot(Spec(j).Sim.cPath(1:end-1) + Spec(j).Sim.aPath(2:end) - (1+Spec(j).r)*Spec(j).Sim.aPath(1:end-1) - Spec(j).Sim.yPath(1:end-1)')


legend(plotLabs, 'Location', 'southoutside', 'Orientation', 'horizontal')
%% Part e


dl_c    = log(cPlot(2:end, :)) -     log(cPlot(1:end-1, :));
eps     = log(yPlot(2:end, 1)) - 0.9*log(yPlot(1:end-1, 1));  % Same shocks for both
vcov0   = cov(dl_c(:, 1), eps, 'omitrows');
vcovNBL = cov(dl_c(:, 2), eps, 'omitrows');

% close all
% scatter(dl_c(:, 1), eps)
% hold on
% scatter(dl_c(:, 2), eps)
% xlabel('\Delta c_t')
% ylabel('\epsilon_t')
% legend({'Borrowing limit', 'No borrowing limit'})




clc
disp('Borrowing limit (psi):')
1 -  vcov0(1,2)  ./ var(eps) 
disp('No borrowing limit (psi):')
1 - vcovNBL(1,2) ./ var(eps)



