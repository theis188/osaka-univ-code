%%   Prepare  Set of Reactions
load([FileName '.mat']);
[ParamInfo, rVnet,Sreg,EnzName,S, V, KVEC, X, SubstrateConc, ParamRange, rV,U, ToRemove] = PrepareLumpedEMVariables(Net);

%% Prepare Solution related functions

% Calculate Dx    
PREDX = S*V;
%Calculate the jacobian matrix
SYMJAC = jacobian(PREDX,X);
%Calculate the flux jacobians
SYMDVDX = jacobian(V,X);
SYMDVDU = jacobian(V,U);

%Calculate the jacobian matrix
SYMJAC = jacobian(PREDX,X);
%Calculate the flux jacobians
SYMDVDX = jacobian(V,X);
SYMDVDU = jacobian(V,U);
%Set up for calculting Vf after sampling the other parameters
SYMVNET = sym('VNET', [length(V) 1]);
SYMK1S = sym('K1S', [size(ParamInfo,1) 1]);
for m=1:size(ParamInfo,1),
    SYMK1S(m) = solve(V(m)-SYMVNET(m), KVEC(ParamInfo(m,1)));
end
% Convert everything into matlab functions for speed
syms t;
DX = matlabFunction(PREDX, 'vars', [{t} {X}, {KVEC} {SubstrateConc} {U}]);
JACOBIAN = matlabFunction(SYMJAC, 'vars', [{t} {X}, {KVEC} {SubstrateConc} {U}]);
K1S = matlabFunction(SYMK1S, 'vars', [{X} {KVEC} {SubstrateConc} {SYMVNET} {U}]);
FLUXES = matlabFunction(V, 'vars', [{X} {KVEC} {SubstrateConc} {U}]);
DVDX = matlabFunction(SYMDVDX, 'vars', [{X}, {KVEC} {SubstrateConc} {U}]);
DVDU = matlabFunction(SYMDVDU, 'vars', [{X}, {KVEC} {SubstrateConc} {U}]);
