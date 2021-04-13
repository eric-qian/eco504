function ResSim = simConsSaving(Spec_j, nSim, nDiscard)

r    = Spec_j.r;
nSim = nSim+1;  % Add padding at the end

% Paths
aPath       = nan(nSim, 1);
cPath       = nan(nSim, 1);
yPathIdx    = nan(nSim, 1);
yPathIdx(1) = 3;

% Store policy functions
aPol = Spec_j.Res.aPol;
cPol = Spec_j.Res.cPol;

% Store grids
yGrid = Spec_j.Res.yGrid;
aGrid  = Spec_j.Res.aGrid;


% Initialize
aPath(1) = 0;
%aPath(1) = aGrid(round(length(aGrid)/2));
aPath(2) = interp1(aGrid, aPol(:, yPathIdx(1)), aPath(1), 'linear', 'extrap');
cPath(1) = interp1(aGrid, cPol(:, yPathIdx(1)), aPath(1), 'linear', 'extrap');



% Stochastic matrix
Pi      = Spec_j.Res.Pi;
PiCumul = cumsum(Pi, 2);

%for j=300:nSim-1
for j=2:nSim-1
    % Draw shock
    yPathIdx(j) = find(PiCumul(yPathIdx(j-1), :) >= rand(), 1, 'first');

    % Compute policy. Interpolate.
    aPath(j+1)  = interp1(aGrid, aPol(:, yPathIdx(j)), aPath(j), 'linear', 'extrap');
    cPath(j)    = interp1(aGrid, cPol(:, yPathIdx(j)), aPath(j), 'linear', 'extrap');  

    
%     plot(cPath)
%     title(j)
%     pause
%     close all
end

yPath = [yGrid(yPathIdx(1:end-1)) NaN];

% Discard draws
if nDiscard > 0
yPath    = yPath(nDiscard:end-1);
aPath    = aPath(nDiscard:end-1);
cPath    = cPath(nDiscard:end-1);
yPathIdx = yPathIdx(nDiscard:end-1);
end


% Store objects
ResSim          = struct;
ResSim.Spec     = Spec_j;
ResSim.nSim     = nSim-1;
ResSim.nDiscard = nDiscard;
ResSim.yPath    = yPath;
ResSim.aPath    = aPath;
ResSim.cPath    = cPath;
ResSim.yPathIdx = yPathIdx;


end
