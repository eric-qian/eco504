function ResSim = simConsSaving(Spec_j, nSim, nDiscard)

r = Spec_j.r;
nSim = nSim+1;  % Add padding at the end

% Paths
aPath       = nan(nSim, 1);
cPath       = nan(nSim, 1);
yPathIdx    = nan(nSim, 1);
yPathIdx(1) = 1;

% Store policy functions
aPol = Spec_j.Res.aPol;


% Store grids
yGrid = Spec_j.Res.yGrid;
aGrid  = Spec_j.Res.aGrid;


% Initialize
aPath(1) = 0;
aPath(2) = interp1(aGrid, aPol(:, yPathIdx(1)), aPath(1), 'linear', 'extrap');
cPath(1) = yGrid(yPathIdx(1)) - aPath(2) + (1+r)*aPath(1);


% Stochastic matrix
Pi      = Spec_j.Res.Pi;
PiCumul = cumsum(Pi, 2);

%for j=300:nSim-1
for j=2:nSim-1
    % Draw shock
    yPathIdx(j) = find(PiCumul(yPathIdx(j-1), :) >= rand(), 1, 'first');

    % Compute policy. Interpolate.
    aPath(j+1)  = interp1(aGrid, aPol(:, yPathIdx(j)), aPath(j), 'linear', 'extrap');
    cPath(j)    = yGrid(yPathIdx(j)) - aPath(j+1) + (1+r)*aPath(j);  
%     plot(cPath)
%     title(j)
%     pause
%     close all
end

yPath = [yGrid(yPathIdx(1:end-1)) NaN];

% Discard draws
yPath    = yPath(nDiscard:end-1);
aPath    = aPath(nDiscard:end-1);
cPath    = cPath(nDiscard:end-1);
yPathIdx = yPathIdx(nDiscard:end-1);



% Store objects
ResSim       = struct;
ResSim.Spec  = Spec_j;
ResSim.nSim  = nSim-1;
ResSim.nDiscard = nDiscard;
ResSim.yPath = yPath;
ResSim.aPath = aPath;
ResSim.cPath = cPath;
ResSim.yPathIdx  = yPathIdx;


end
