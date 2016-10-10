%% Model preparation for KAF of Fdh

ModeOpts.ImportParameters = [];
nEnz = 60;
Uini = [ones(nEnz,1); 0];
ModeOpts.Perts = repmat({Uini}, 1,1);

%% Make Perts
%% Make Perts

%Fdh only JCL16F
ModeOpts.Perts{1} = [ones(nEnz,1); 100];
%%%
%Fdh + Frd, Ldh, Adh KO JCL166F
ModeOpts.Perts{2} = [ones(nEnz,2); 0 100];
ModeOpts.Perts{2}(35,:) = .01; % KO succinate 35 25
ModeOpts.Perts{2}(39,:) = .01; % KO succinate 39 26
ModeOpts.Perts{2}(17,:) = .01; % KO lactate   17 9
ModeOpts.Perts{2}(14,:) = .01; % KO ethanol   14 7

%Fdh + Frd, Ldh, Adh, KO + Pta KO JCL299F
ModeOpts.Perts{3} = [ones(nEnz,4); 0 10 20 100];
ModeOpts.Perts{3}(35,:) = .01; % KO succinate
ModeOpts.Perts{3}(39,:) = .01; % KO succinate
ModeOpts.Perts{3}(17,:) = .001; % KO lactate
ModeOpts.Perts{3}(14,:) = .01; % KO ethanol
ModeOpts.Perts{3}(31,:) = [1 0.01 0.001 0.001]; % KO Pta 31 21

%Fdh + Frd, Ldh, Adh, KO + Pta KO + Adh OE JCL299FT (RESTORE PTA @ 10%)
ModeOpts.Perts{4} = [ones(nEnz,4); 0 10 20 100];
ModeOpts.Perts{4}(54,:) = 2; % OE bdh
%ModeOpts.Perts{1}(14,:) = 2; % OE adh
ModeOpts.Perts{4}(35,:) = .01; % KO succinate
ModeOpts.Perts{4}(39,:) = .01; % KO succinate
ModeOpts.Perts{4}(17,:) = .001; % KO lactate
ModeOpts.Perts{4}(14,:) = .01; % KO ethanol
ModeOpts.Perts{4}(31,:) = [1 0.1 0.1 0.1]; % KO Pta 31 21

%Fdh + Frd, Ldh, Adh, KO + Pta KO + Adh OE + Ack KO + EutD KO + AtoB OE SASTIA's PREDICTION STRAIN
ModeOpts.Perts{5} = [ones(nEnz,4); 0 10 20 100];
ModeOpts.Perts{5}(54,:) = 2; % OE bdh
%ModeOpts.Perts{2}(14,:) = 10; % OE adh
ModeOpts.Perts{5}(48,:) = 10; % OE AtoB
ModeOpts.Perts{5}(35,:) = .01; % KO succinate
ModeOpts.Perts{5}(39,:) = .01; % KO succinate
ModeOpts.Perts{5}(17,:) = .001; % KO lactate
ModeOpts.Perts{5}(14,:) = .01; % KO ethanol
ModeOpts.Perts{5}(10,:) = .01; % KO ack
ModeOpts.Perts{5}(31,:) = [1 0.01 0.001 0.001]; % KO Pta 31 21

clear OriginalModeOpts Uini n
tic
% myCluster = parcluster('local');
% myCluster.NumWorkers = 40;
% 
% matlabpool open 39
% pctRunOnAll warning off all
MainAddFlux('ModelNoCOA',500,'SeveralPerturbations',ModeOpts)
matlabpool close
toc
%EmailSender('lafonj@gmail.com','SimulationDone','Finished',[])