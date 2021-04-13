function Res = solveConsSaving(Spec, method)
%% Unpack parameters

if ~exist('method')
    method = 'EGM';
end

% Parameters
phi   = Spec.phi;   % Borrowing limit
gamma = Spec.gamma; % CRRA parameter
beta  = Spec.beta;
r     = Spec.r;
epsConv = Spec.epsConv;

% Income process parameters
% w_t = wBar + rho w_t-1 + eps_t, eps_t~N(0, vareps)
vareps = Spec.vareps;  % Variance of epsilon of w process
m      = Spec.m;       % grid is 0 +/- m*sigma^2_w
rho    = Spec.rho;     % AR coefficient of w process
Nw     = Spec.Nw;
Na     = Spec.Na;


%% Make income grid
wBar       = -vareps/(2*(1+rho));
varw       = vareps/(1-rho^2);  % Unconditional variance of w_t
muw        = wBar/(1-rho);  % Unconditional mean

% Make grid using the Tauchen method
ResGrid = tauchen(rho, wBar, vareps, Nw, m);

% Stochastic transition matrix
Pi    = ResGrid.Pi;
yGrid = exp(ResGrid.wGrid);  % Transform to income units


% Upper bound of income process
upper = yGrid(end)/r;

if phi    == 'NBL'  % Natural borrowing limit
    lower = -yGrid(1)/r;
else
    lower = phi;    % Exogenous borrowing limit
end

% Asset grid
aGrid = (exp(linspace(log(1), log(2), Na))-1)*(upper-lower)+lower;
%aGrid = linspace(lower,upper, Na); 


aaGrid = repmat(aGrid(:), 1, Nw);  % Asset grid in matrix
yyGrid = repmat(yGrid(:)', Na, 1);  % Income grid in matrix form



%% EGM

if strcmp(method, 'EGM')
    converged = 0;
    iter      = 1;
    
    cNew = aaGrid*r + yyGrid;  % Consumption policy function
    aNew = nan(size(cNew));    % Savings policy function
    err = NaN;
    
    while converged == 0
        
        % Print to console
        if mod(iter, 25) == 0
            clc
            disp(['Iteration ' num2str(iter), ', err = ' num2str(err) '...'])
        end
        
        
        % Consumption policy function. Na x Ny
        cOld = cNew;
        aOld = aNew;
        

        B      = beta*(1+r) * dcrra(cOld, gamma)* Pi';  % RHS of EE
        cTilde = invdcrra(B, gamma);                    % Solve for consumption                


        aStar  = (cTilde + aaGrid - yyGrid)/(1+r);     % Current asset holdings leading consumer to hold aaGrid tomorrow
                        
        cNew = nan(size(cOld));  % Consumption policy rule
        
        % Interpolate points
        for jy = 1:Nw  % Loop through shocks
            aStar_j     = aStar(:, jy);
            cTilde_j    = cTilde(:, jy);
            cInterp     = griddedInterpolant(aStar_j, cTilde_j, 'linear', 'linear');            
            cNew(:, jy) = cInterp(aGrid);
            
            
            % If borrowing constraint binds a < min(aStar_j)
             jBind          = aGrid < aStar_j(1);
             cNew(jBind,jy) = (1+r)*aGrid(jBind) + yGrid(jy) - aGrid(1);
            
        end
        
        aNew = yyGrid -  cNew + (1+r)*aaGrid;
        
        err = max(abs(aNew(:)-aOld(:))./aOld(:)*100);
        if err < epsConv
            converged = 1;
        end
        
        iter = iter +1;
    end

    
    %% Value function iteration
else
    Nh      = 0;
    aNew    = aaGrid;
    aOld    = nan(size(aNew));
    VNew    = zeros(size(aNew));  % Value function
    VOld    = VNew;
    cNew    = NaN;
    yyyGrid = repmat(yyGrid, 1,1, Na);
    aaaGrid = repmat(aaGrid, 1,1,Na);
    aapGrid = permute(aaaGrid, [3 2 1]);
    
    converged = 0;
    iter      = 1;
    err       = nan;
    
    while converged == 0
        if mod(iter, 25) == 0
            clc
            disp(['Iteration ' num2str(iter), ', err = ' num2str(err) '...'])
        end
        
        aOld = aNew;
        cOld = cNew;
        VOld = VNew;
        
        % Continuation value. 
        Cont      = beta * Pi * VOld';
        Cont      = repmat(Cont, 1, 1, Na);
        Cont      = permute(Cont, [3,1,2]);
        
        
        VVNew     = crra(yyyGrid + (1+r)*aaaGrid - aapGrid, gamma) + Cont;
        [VNew, j] = max(VVNew, [], 3, 'linear');
        aNew      = aapGrid(j);

        for jh = 1:Nh
            VNew = crra(yyGrid + (1+r)*aaGrid - aNew, gamma) + beta * VNew*Pi';
        end
        
        
        err = max(abs((VNew(:)-VOld(:))./VOld(:))*100);
        if err < epsConv
            converged = 1;
        end
        
        iter = iter +1;
                
        cNew = yyGrid -  aNew + (1+r)*aaGrid;        
        
    end


end

%% Prepare output object

Res       = struct;
Res.Spec  = Spec;
Res.Pi    = Pi;
Res.yGrid = yGrid;
Res.aGrid = aGrid;
Res.cPol  = cNew;
Res.aPol  = aNew;


