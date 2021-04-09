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

ResGrid = tauchen(rho, wBar, vareps, Nw, m);

Pi    = ResGrid.Pi;
yGrid = exp(ResGrid.wGrid);

if phi    == 'NBL'
    aGrid = linspace(-yGrid(1)/r, yGrid(end)/r, Na);
else
    aGrid = linspace(phi, yGrid(end)/r, Na);
end




aaGrid = repmat(aGrid(:), 1, Nw);  % Asset grid in matrix
yyGrid = repmat(yGrid(:)', Na, 1);

converged = 0;
iter      = 1;

%% EGM

if strcmp(method, 'EGM')
    
    cNew = aaGrid*r + repmat(yGrid(:)', Na, 1);
    aNew = nan(size(cNew));
    
    
    while converged == 0
        
        % Print to console
        if mod(iter,100) == 0
            clc
            disp(['Iteration ' num2str(iter), ', Delta = ' num2str(max(abs(cNew(:)-cOld(:)))) '...'])
        end
        
        % Consumption policy function. Na x Ny
        cOld = cNew;
        aOld = aNew;
        
        B      = beta*(1+r) * crra(cOld, gamma)* Pi';  % RHS of EE
        cTilde = invCRRA(B, gamma);                    % Solve for consumption
        
        aStar = (cTilde + aaGrid - yyGrid)/(1+r);  % Current asset holdings
        
        
        cNew = nan(size(cOld));
        % Interpolate points
        for jy = 1:Nw
            aStar_j = aStar(:, jy);
            cNew(:, jy) = interp1(aStar_j, cTilde(:, jy), aGrid, 'linear', 'extrap');
            
            
            % If borrowing constraint binds
            jBind          = aGrid < aStar_j(1);
            cNew(jBind,jy) = (1+r)*aGrid(jBind) + yGrid(jy) - aGrid(1);
            
        end
        
        aNew = yyGrid -  cNew + (1+r)*aaGrid;
        
        
        if max(abs(aNew(:)-aOld(:))) < epsConv
            converged = 1;
        end
        
        iter = iter +1;
        
    end
    
    %% Value function iteration
else
    
    aNew    = aaGrid;
    aOld    = nan(size(aNew));
    VNew    = ones(size(aNew));  % Value function
    cNew    = NaN;
    yyyGrid = repmat(yyGrid, 1,1, Na);
    aaaGrid = repmat(aaGrid, 1,1,Na);
    aapGrid = permute(aaaGrid, [3 2 1]);
    
    while converged == 0
        if mod(iter,1) == 0
            clc
            disp(['Iteration ' num2str(iter), ', Delta = ' num2str(max(abs(aNew(:)-aOld(:)))) '...'])
        end
        
        aOld = aNew;
        VOld = VNew;
        
        % Continuation value. 
        Cont = beta * Pi * VOld';
        Cont = repmat(Cont, 1, 1, Na);
        Cont = permute(Cont, [3,1,2]);
        [VNew, j] = max(crra(yyyGrid - (1+r)*aaaGrid - aapGrid, gamma) + Cont, [], 3, 'linear');
        aNew = aaaGrid(j);

        
        
        if max(abs(VNew(:)-VOld(:))) < epsConv
            converged = 1;
        end
        
        iter = iter +1;
        
        
        
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


