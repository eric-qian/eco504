% mainProblem2.m  Main file for problem 2
% TODO: SIMULATION OF LOG CONSUMPTION TENDS TO INFINITY
%% Preliminaries
clear
clc

addpath('functions/')

% Structure giving specifications
Spec = struct;

%% Baseline specification
j=1;
Spec(j).lab = 'Baseline';

% Parameters
Spec(j).phi     = 0;  % Borrowing limit
Spec(j).gamma   = 2;  % CRRA parameter
Spec(j).beta    = 0.95; 
Spec(j).r       = 0.02;
Spec(j).Na      = 10;  
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
Spec(j).gamma = 1;

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

for j = 1
%for j = 1:length(Spec)
    Spec_j      = Spec(j);
    Spec(j).Res = solveConsSaving(Spec_j, 'waka');
end

return
%% Simulate model

nSim = 100000;

for jSpec = 1:length(Spec)
    Spec_j      = Spec(jSpec);
    rng(1)
    Spec(jSpec).Sim = simConsSaving(Spec_j, nSim);
end

%% Part b
plotIdx  = [2, 1, 3];
plotLabs = {'\gamma=1', '\gamma=3', '\gamma=5'};

cPlot    = []; 
aPlot    = []; 

for j =plotIdx
 cPlot = [cPlot Spec(j).Sim.cPath];
  aPlot = [aPlot Spec(j).Sim.aPath];

end


plot(cPlot)
std(cPlot, 0, 1, 'omitnan')
legend(plotLabs, 'Location', 'southoutside', 'Orientation', 'horizontal')


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
std(cPlot, 0, 1, 'omitnan')
legend(plotLabs, 'Location', 'southoutside', 'Orientation', 'horizontal')


%% Part d
plotIdx  = [1, 6];
plotLabs = {'Baseline', 'NBL'};

cPlot    = []; 
yPlot    = []; 

for j =plotIdx
    cPlot = [cPlot Spec(j).Sim.cPath];
    yPlot = [yPlot Spec(j).Sim.yPath'];
end


plot(cPlot)
std(cPlot, 0, 1, 'omitnan')
mean(cPlot, 'omitnan')
legend(plotLabs, 'Location', 'southoutside', 'Orientation', 'horizontal')


%% Part e


dl_c = log(cPlot(2:end, :)) - log(cPlot(1:end-1, :));
eps = log(yPlot(2:end, :)) - 0.9*(yPlot(1:end-1, :));
1- cov(dl_c(:, 1), eps(:, 1), 'omitrows') ./ var(eps, 'omitnan')
1- cov(dl_c(:, 2), eps(:, 2), 'omitrows') ./ var(eps, 'omitnan')



