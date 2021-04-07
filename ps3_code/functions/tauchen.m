function Res = tauchen(rho, wBar, vareps, N, m)
% tauchen()  Tauchen method for an AR(1), with constant
% Process is w_t = wBar + rho w_t-1 + eps_t, 
% Input
%  - rho:    AR coefficient
%  = wBar:   Drift term
%  - vareps: variance of eps_t
%  - N:      Size of grid
%  - m:      Multiple used for constructing the grid


% Unconditionals
varw = vareps/(1-rho^2);  % Unconditional variance of w_t
muw  = wBar/(1-rho);      % Unconditional mean

% Make grid, centered at unconditional mean
wGrid = linspace(muw-(varw^0.5)*m, muw + (varw^0.5)*m, N); 
d     = wGrid(2) - wGrid(1);


% Make stochastic matrix
M  = repmat(wGrid(:)', N, 1) - rho*repmat(wGrid(:), 1, N);  % y_k - rho y_j 
Pi = normcdf(M+d/2, wBar, sqrt(vareps)) -...
     normcdf(M-d/2, wBar, sqrt(vareps));

% Correct boundaries
Pi(:, 1)   =     normcdf(M(:, 1)   + d/2, wBar, sqrt(vareps));
Pi(:, end) = 1 - normcdf(M(:, end) - d/2, wBar, sqrt(vareps));


% Output
Res       = struct;
Res.wGrid = wGrid;
Res.Pi    = Pi;
Res.varw  = varw;

% Save inputs 
Res.in         = struct;
Res.in.rho     = rho;
Res.in.vareps  = vareps;
Res.in.N       = N;
Res.in.m       = m;


end