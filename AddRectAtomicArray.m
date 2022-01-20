function AddRectAtomicArray(LAtoms, WAtoms, X0, Y0, VX0, VY0, InitDist, Temp, Type)
%% Introduction: Create a rectangular array of atoms
% 1. Inputs:
%   - LAtoms: # atoms in x-direction
%   - WAtoms: # atoms in y-direction
%   - X0,Y0: center of rectangle
%   - VX0,VY0: average (group) velocity of atoms
%   - initDist: initial disturbance of atoms
%   - Temp: temperature of mass
%   - Type: atom (type 1 or 2)
% 2. Outputs (global):
%   - x,y: add new atoms position to this
%   - Vx,Vy: add new atoms velocities to this
%   - nAtoms: updates total
% 3. Functions:
%   - linspace(): uniform distribution of points
%   - rand(): random number from a flat distribution
%   - randn(): random number from a normal distribution


%% Set variables 
global C
global x y AtomSpacing
global nAtoms
global AtomType Vx Vy Mass0 Mass1

%% Set atom types 
if Type == 0
    Mass = Mass0;
else
    Mass = Mass1;
end

%% Intial positions
% Size of the rectangular box of atoms = (# of atom - 1) * atom spacing
L = (LAtoms - 1) * AtomSpacing;  
W = (WAtoms - 1) * AtomSpacing;

numAtoms = LAtoms * WAtoms;

xp(1, :) = linspace(0, L, LAtoms); % generate LAtoms equally spaced between 0 and L
yp(1, :) = linspace(0, W, WAtoms);

x(nAtoms + 1:nAtoms+LAtoms) = xp-L/2;
y(nAtoms + 1:nAtoms+LAtoms) = yp(1)-W/2;

for i = 1:WAtoms-1
    x(nAtoms + i * LAtoms + 1:nAtoms + (i + 1) * LAtoms) = xp - L / 2;
    y(nAtoms + i * LAtoms + 1:nAtoms + (i + 1) * LAtoms) = yp(i + 1) - W / 2;
end

%% Disturb positions: so that the atom is jiggling 
x(nAtoms + 1:nAtoms + numAtoms) = x(nAtoms + 1:nAtoms + numAtoms) + ...
    (rand(1, numAtoms) - 0.5) * AtomSpacing * InitDist + X0;
y(nAtoms + 1:nAtoms + numAtoms) = y(nAtoms + 1:nAtoms + numAtoms) + ...
    (rand(1, numAtoms) - 0.5) * AtomSpacing * InitDist + Y0;

AtomType(nAtoms + 1:nAtoms + numAtoms) = Type;

%% Calculate thermal velocities
if Temp == 0
    Vx(nAtoms + 1:nAtoms + numAtoms) = 0;
    Vy(nAtoms + 1:nAtoms + numAtoms) = 0;
else
    std0 = sqrt(C.kb * Temp / Mass);

    Vx(nAtoms + 1:nAtoms + numAtoms) = std0 * randn(1, numAtoms);
    Vy(nAtoms + 1:nAtoms + numAtoms) = std0 * randn(1, numAtoms);
end

%% Group velocity
Vx(nAtoms + 1:nAtoms + numAtoms) = Vx(nAtoms + 1:nAtoms + numAtoms) - ...
    mean(Vx(nAtoms + 1:nAtoms + numAtoms)) + VX0;
Vy(nAtoms + 1:nAtoms + numAtoms) = Vy(nAtoms + 1:nAtoms + numAtoms) - ...
    mean(Vy(nAtoms + 1:nAtoms + numAtoms)) + VY0;

nAtoms = nAtoms + numAtoms;

end
