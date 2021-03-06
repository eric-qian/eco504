% mainProblem2.m  Main file for problem 2
%% Preliminaries
clear
clc

addpath('functions/')

% Structure giving specifications
Spec = struct;


Nw = 5;  % # Grid points for income shock
%Nw = 101;  % # Grid points for income shock, big grid

% Make figures
figPath = ['figures_Nw=' num2str(Nw) '/'];
mkdir(figPath)
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
Spec(j).epsConv = 5e-7;
Spec(j).method = 'EGM';

% Income process parameters
% w_t = wBar + rho w_t-1 + eps_t, eps_t~N(0, vareps)
Spec(j).vareps  = 0.06;                 % Variance of epsilon of w process
Spec(j).m       = 3;                    % grid is 0 +/- m*sigma^2_w
Spec(j).rho     = 0.9;                  % AR coefficient of w process
Spec(j).Nw      = Nw;                    % Grid size

%% Alternative specifications

% Changing gamma, gamma = 1
j             = 2;
Spec(j)       = Spec(1);
Spec(j).lab   = '\gamma=1';
Spec(j).gamma = 1.0;


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


% Imposing the natural borrowing limit
j              = 6;
Spec(j)        = Spec(1);
Spec(j).lab    = '\phi=-y_{min}/r';
Spec(j).phi    = 'NBL';


%% Estimate models

for j = 1:length(Spec)
    Spec_j      = Spec(j);
    Spec(j).Res = solveConsSaving(Spec_j, Spec_j.method);
end


%% Simulate model time paths


for jSpec = 1:length(Spec)
    Spec_j      = Spec(jSpec);
    rng(1)
    Spec(jSpec).Sim = simConsSaving(Spec_j, nSim, nDiscard);
end


%% Part a: Plot policy function

figure

surf(repmat(Spec(1).Res.aGrid', 1, Spec(1).Nw),...
    repmat(Spec(1).Res.yGrid, Spec(1).Na, 1), Spec(1).Res.aPol)
shading interp
xlabel('a')
ylabel('y')
zlabel('a''')
box on
grid on
resizeFig(figSize)
saveas(gcf, [figPath 'p2_parta.png'])

%% Part b

plotIdx  = [2, 1, 3];

cPlot     = []; 
aPlot     = [];


for j =plotIdx
    cPlot = [cPlot Spec(j).Sim.cPath];
    aPlot = [aPlot Spec(j).Sim.aPath];
    
end

sigmaVec = std(cPlot);
plotLabs = {['\gamma=1 (\sigma=' num2str(sigmaVec(1), 3) ')'], ...
            ['\gamma=2 (\sigma=' num2str(sigmaVec(2), 3) ')'],...
            ['\gamma=5 (\sigma='  num2str(sigmaVec(3), 3) ')']};

close all
plot(cPlot)
xlim([0, nSim-nDiscard])
legend(plotLabs, 'Location', 'southoutside', 'Orientation', 'horizontal')
formatFig(figSize)
saveas(gcf, [figPath 'p2_partb.png'])

%% Part c
plotIdx  = [4, 1, 5];


cPlot    = []; 
yPlot    = []; 
aPlot    = [];

for j =plotIdx
    cPlot = [cPlot Spec(j).Sim.cPath];
    yPlot = [yPlot Spec(j).Sim.yPath(:)];
    aPlot = [aPlot Spec(j).Sim.aPath(:)];
end

close all
sPlot = (yPlot + (1+Spec(1).r)*aPlot - cPlot) ./((1+Spec(1).r)*aPlot+ yPlot);
mu_s = mean(sPlot);


plotLabs = {['\sigma_\epsilon^2=0.01 (\mu=' num2str(mu_s(1),3) ')'],...
    ['\sigma_\epsilon^2=0.06 (\mu=' num2str(mu_s(2),3) ')'], ...
    ['\sigma_\epsilon^2=0.12 (\mu=' num2str(mu_s(3),3) ')']};
plot(sPlot)
xlim([0, nSim-nDiscard])
legend(plotLabs, 'Location', 'southoutside', 'Orientation', 'horizontal')
formatFig(figSize)
saveas(gcf, [figPath 'p2_partc.png'])


%% Parts d and e
plotIdx  = [1, 6];

cPlot    = []; 
yPlot    = []; 
aPlot    = [];

for j =plotIdx
    cPlot = [cPlot Spec(j).Sim.cPath];
    aPlot = [aPlot Spec(j).Sim.aPath];    
    yPlot = [yPlot Spec(j).Sim.yPath'];
end


% Compute insurance coefficient
dl_c    = log(cPlot(2:end, :)) -     log(cPlot(1:end-1, :));
eps     = log(yPlot(2:end, 1)) - 0.9*log(yPlot(1:end-1, 1));  % Same shocks for both
vcov0   = cov(dl_c(:, 1), eps, 'omitrows');
vcovNBL = cov(dl_c(:, 2), eps, 'omitrows');
psi0    = 1-vcov0(2,1)/vcov0(2,2);
psiNBL    = 1-vcovNBL(2,1)/vcovNBL(2,2);


% Get mean
muVec = mean(cPlot);


% Make plot
plotLabs = {['\phi=0 (\mu=' num2str(muVec(1), 3) ', \psi=' num2str(psi0, 3)  ')'], ...
    ['\phi=NBL (\mu=' num2str(muVec(2),3) ', \psi=' num2str(psiNBL, 3) ')' ]};

close all
plot(cPlot)
xlim([0, nSim-nDiscard])
legend(plotLabs, 'Location', 'southoutside', 'Orientation', 'horizontal')
formatFig(figSize)
saveas(gcf, [figPath 'p2_partd.png'])


%% Appendix: Compare solution methods

compare = 1;

if compare == 1
    j=1;
    
    Spec_j       = Spec(j);
    Spec_j.Na    = 400;  % Make the grid smaller so that the problem is tractable
    Spec_j.gamma = 1;
%    Spec_j.phi   = 'NBL';
    
    SpecEG  = Spec_j;
    SpecVFI = Spec_j;
    
    SpecEG.Res   = solveConsSaving(SpecEG);
    SpecVFI.Res  = solveConsSaving(SpecVFI, 'VFI');
        
    
    % Save policy functions
    aPolEG  = SpecEG.Res.aPol;
    aPolVFI = SpecVFI.Res.aPol;
    aGrid   = SpecEG.Res.aGrid;
     
    cPolEG   = SpecEG.Res.cPol;
    cPolVFI = SpecVFI.Res.cPol;
   
        
    % Compare policy functions
    close all
    subplot(2,1,1)
    plot(aGrid, cPolEG, 'Color', 'b')
    hold on
    plot(aGrid, cPolVFI, 'Color', 'r', 'LineStyle', '--')
    title('c policy')
    
    subplot(2,1,2)
    plot(aGrid, aPolEG, 'Color', 'b')
    hold on
    plot(aGrid, aPolVFI, 'Color', 'r', 'LineStyle', '--')
    title('a policy')
    legend({'EG', 'VFI'}, 'Location', 'best')
    formatFig([6,4])
    saveas(gcf, [figPath 'VFI_EG_comparePolicy.png'])        
    
    
    % Compare time paths
    rng(1)
    ResEG.Sim  = simConsSaving(SpecEG, nSim, nDiscard);
    rng(1)
    ResVFI.Sim = simConsSaving(SpecVFI, nSim, nDiscard);
        
    close all
    subplot(2,1,1)
    plot(ResEG.Sim.aPath, 'Color', 'b')
    title('a path')
    hold on
    plot(ResVFI.Sim.aPath, 'Color', 'r', 'LineStyle', ':')
    xlim([0, nSim-nDiscard])

    subplot(2,1,2)
    plot(ResEG.Sim.cPath, 'Color', 'b')
    hold on
    plot(ResVFI.Sim.cPath, 'Color', 'r', 'LineStyle', ':')
    title('c path')
    xlim([0, nSim-nDiscard])

    formatFig([6,4])

    
    lgd = legend({'EG', 'VFI'}, 'Location', 'southoutside', ...
        'Orientation', 'horizontal', 'Box', 'off');
    lgd.Position= [0.3900    0.0083    0.2500    0.0521];
    
    saveas(gcf, [figPath 'VFI_EG_compare.png'])        
    
end
