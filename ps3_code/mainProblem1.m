% mainProblem1.m  Main file for problem 1.
%% Preliminaries
clear
clc
addpath('functions/')


% Tauchen Parameters
rho     = 0.9;
vareps  = 0.06;
m       = 4;  % grid is 0 +/- sigma^2_w
wBar    = -vareps/(2*(1+rho));
varw    = vareps/(1-rho^2);  % Unconditional variance of w_t
muw     = wBar/(1-rho);  % Unconditional mean


% Simulation/plot parameters
nSim    = 10000;
nBurn   = nSim/10; 
figSize = [6,4];
figPath = 'figures/';

%% Run Tauchen
ResTauchen5    = tauchen(rho, wBar, vareps, 5, m);  % Tauchen, 5 states
ResTauchen10   = tauchen(rho, wBar, vareps, 10, m);  
ResTauchen1000 = tauchen(rho, wBar, vareps, 1000, m);  

% Simulate Tauchen. For consistency, use the same seed
rng(1)
wSim5    = simTauchen(ResTauchen5, nSim);
rng(1)
wSim10   = simTauchen(ResTauchen10, nSim);
rng(1)
wSim1000 = simTauchen(ResTauchen1000, nSim);


% w to income space
ySim5    = exp(wSim5);
ySim10   = exp(wSim10);
ySim1000 = exp(wSim1000);

mean(ySim5)
mean(ySim10)
mean(ySim1000)


%% make plots
close all

% Make plot for process w
figure
subplot(3,1,1)
plot(wSim5(nBurn:end));
title(['5 states (\mu_w=' num2str(mean(wSim5), 3) ', \sigma_w^2=' num2str(var(wSim5), 3), ')'])

subplot(3,1,2)
plot(wSim10(nBurn:end));
title(['10 states (\mu_w=' num2str(mean(wSim10), 3) ', \sigma_w^2=' num2str(var(wSim10), 3), ')'])

subplot(3,1,3)
plot(wSim1000(nBurn:end));
title(['1000 states (\mu_w=' num2str(mean(wSim1000), 3) ', \sigma_w^2=' num2str(var(wSim1000), 3), ')'])
sgtitle(['Analytical unconditional moments (w): \mu_w=' num2str(muw, 3), ', \sigma^2_w=' num2str(varw, 3)])
formatFig(figSize)
saveas(gcf, [figPath, 'p1_wTS.png'])


% Make plot for process y
figure
subplot(3,1,1)
plot(ySim5(nBurn:end));
title(['5 states (\mu_y=' num2str(mean(ySim5), 3) ')'])

subplot(3,1,2)
plot(ySim10(nBurn:end));
title(['10 states (\mu_y=' num2str(mean(ySim10), 3) ')'])

subplot(3,1,3)
plot(ySim1000(nBurn:end));
title(['1000 states (\mu_y=' num2str(mean(ySim1000), 3) ')'])
formatFig(figSize)
saveas(gcf, [figPath, 'p1_yTS.png'])



disp('Unconditional [mean, variance], 5 states')
disp([mean(wSim5) var(wSim5)])


disp('Unconditional [mean, variance], 10 states')
disp([mean(wSim10) var(wSim10)])

disp('Unconditional [mean, variance], analytical')
disp([muw varw])

