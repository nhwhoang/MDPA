doPlot = 1;
dt = 5e-15;
TStop = 3000 * dt;  % stop time (in this case 3000 steps)
InitDist = 0.0;
Method = 'VE'; % VE -- verlot; FD -- Forward Difference (FD is simple but has problem that no conservation of E)

Mass0 = 14 * C.am; % Silicon
Mass1 = 100 * C.am; % Argon
Mass2 = 10 * C.am; % not sure what this is

% Set up parameters for LJ Potnetial distribution
AtomSpacing = 0.5430710e-9;  % also the LJ Potential
LJSigma = AtomSpacing / 2^(1 / 6);
LJEpsilon = 1e-21;

PhiCutoff = 3 * AtomSpacing * 1.1; % cut off (not calculate force) if more than 3 atoms away

T = 30;

% Change type 0 array
AddCircAtomicArray(2*pi, 0, 0, 0, 0, 0, T, 0);

%AddRectAtomicArray(10, 10, 0, 0, 0, 0, 0, T, 0);
% vy0 = -sqrt(0.02*Ep/Mass1);
% AddRectAtomicArray(4,4,0,12*AtomSpacing,0,vy0,0,T,1);

% Fire particle stream at the above array
Ep = 2;
AddParticleStream(5, 0.1, 10, -pi / 2, 1, Ep * C.q_0, 5);

%%% ADD THIS IN FOR PA2
AddParticleStream(6, 0.1, -10, pi / 2, 2, Ep * C.q_0, 6);

% Size and limit of the simulation region
Size = 10*AtomSpacing;
Limits = [-Size +Size -Size +Size]; % square is good
PlDelt = 5 * dt;

PlotFile = 'BlockSt.gif';
PlotPosOnly = 1;
doPlotImage = 0;
PlotSize = [100, 100, 1049, 1049];

ScaleV = .02e-11;
ScaleF = 10;
