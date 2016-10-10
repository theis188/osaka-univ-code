function [ParamInfo, rVnet,Sreg, EnzName,S, V, KVEC, X, SubstrateConc, ParameterRange, rV, U, ToRemove]=PrepareLumpedEMVariables(Net)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The code below sets up all the variables to be used later. There is no
%need to change this
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% input the stoichiometric matrix along with required indices and fluxes

S=Net.S;
Sreg=Net.Sreg;     %%%%%%%% This script will not use regulation
EnzName=Net.EnzName;
S2a = [];
EnzName2 = [];
if isfield(Net, 'S2a')
    S2a = Net.S2a;
    EnzName2 = Net.EnzName2;
end
%% Identify passive transport 
Vout_index = [];
Vin_index = [];
for n=1:length(EnzName)
    if sum(S(:,n)>0) == 0,
        Vout_index = [Vout_index; n];
    end
    if sum(S(:,n)<0) == 0,
        Vin_index = [Vin_index; n];
    end
end
% input the thermodynamic constraints
rVnet=Net.Vref;
Revs = Net.Reversibilities;

%% Create the system of ODEs
TotalParams = FindTotalParams(S, Sreg, Vin_index, Vout_index,Revs,EnzName);
KVEC = sym('KVEC', [TotalParams+length(EnzName2) 1]);
U = sym('U', [length(EnzName)+length(EnzName2) 1]);
X = sym('X', [size(S,1) 1]);

[V, ParamInfo, SubstrateConc, ParameterRange, rV] = SetUpRxnRates(S,Sreg, KVEC, X, Vin_index, Vout_index, U,rVnet,Revs, S2a,EnzName,Net.MetabName);

%% Find conserved moeities and replace them with algebraic expressions
NullMat = null(S','r');
ToRemove = [];
for n=1:size(NullMat,2),
    %Replace every set of conserved moeties with an algebraic
    %expression
    CurrentVec = NullMat(:,n);
    OtherVec = NullMat;
    OtherVec(:,n) = [];
    tempRemov = find(CurrentVec~=0 & all(OtherVec==0,2));
    ToRemove = [ToRemove tempRemov(1)];
    Expression = solve(CurrentVec'*X - sum((CurrentVec)), X(tempRemov(1)));
    V = subs(V, X(tempRemov(1)), Expression);
end
% MetabName(ToRemove) = [];
S(ToRemove,:) = [];
if ~isempty(S2a),
    S2a(ToRemove,:) = [];
end
X(ToRemove) = [];

%Merge stoichiometries
S = [S S2a];
rVnet = [rVnet; zeros(size(S2a,2),1)];
EnzName = [EnzName; EnzName2];


