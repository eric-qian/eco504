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
Spec(j).Na      = 20000;  
Spec(j).epsConv = 1e-6;

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
Spec(j).gamma = 1.01;

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
    
    Spec_j = Spec(j);
    ResEG  = solveConsSaving(Spec_j);
    ResVFI = solveConsSaving(Spec_j, 'VFI');
    pctDiff = (ResVFI.aPol - ResEG.aPol) ./ ResEG.aPol;
    nanmean(pctDiff(:)*100)
    nanmedian(pctDiff(:)*100)
    min(pctDiff(:)*100)
end


for j = 1:length(Spec)
    Spec_j      = Spec(j);
    Spec(j).Res = solveConsSaving(Spec_j);
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

sPlot = ones(size(cPlot))- cPlot./yPlot;
plot(sPlot)
xlim([0, nSim-nDiscard])
legend(plotLabs, 'Location', 'southoutside', 'Orientation', 'horizontal')
saveas(gcf, [figPath 'p2_partc.png'])


%% Part d
plotIdx  = [1, 6];

cPlot    = []; 
yPlot    = []; 

for j =plotIdx
    cPlot = [cPlot Spec(j).Sim.cPath];
    yPlot = [yPlot Spec(j).Sim.yPath'];
end

muVec = mean(cPlot, 'omitnan');
plotLabs = {['Baseline (\mu=' num2str(muVec(1), 3) ')'], ...
    ['NBL (\mu=' num2str(muVec(2),3) ')' ]};

plot(cPlot)
xlim([0, nSim-nDiscard])
legend(plotLabs, 'Location', 'southoutside', 'Orientation', 'horizontal')
saveas(gcf, [figPath 'p2_partd.png'])




%% Part e


dl_c    = log(cPlot(2:end, :)) -     log(cPlot(1:end-1, :));
eps     = log(yPlot(2:end, 1)) - 0.9*log(yPlot(1:end-1, 1));  % Same shocks for both
vcov0   = cov(dl_c(:, 1), eps, 'omitrows');
vcovNBL = cov(dl_c(:, 2), eps, 'omitrows');

clc
1 -  vcov0(1,2)  ./ var(eps, 'omitnan')
1 - vcovNBL(1,2) ./ var(eps, 'omitnan')



