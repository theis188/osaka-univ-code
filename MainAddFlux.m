function [SaveFileName] = MainAddFlux(FileName, EnsembleSize, Mode, ModeOpts)
%% Currently Supported Modes:
    % 'SeveralPerturbations',
    % ModeOpts.Perts,ModeOpts.ImportParameters,ModeOpts.Verbose
    % 'SearchRobust',ModeOpts.Perts

    %% Prepare Set of Reactions
MatlabVersion = version;
tmp = strfind(MatlabVersion,'.');
MatlabVersion = MatlabVersion(1:tmp(2)-1);
ModelFile = [FileName '_PrepV' MatlabVersion '.mat'];
if exist(ModelFile, 'file')==2,
    load(ModelFile, '-regexp', '^(?!EnsembleSize$|Mode$|ModeOpts$).');
else    
    PrepareRxns
    save( ModelFile, '-regexp', '^(?!(Mode|ModeOpts|m)$).' )
end

if ~isfield(ModeOpts, 'Verbose')
    ModeOpts.Verbose = 1;
end

if ModeOpts.Verbose,
    fprintf(1,'Preparation is done\n');
end

NoEnzymes = length(EnzName);
Steps = 10;
StepsUp = Steps;
StepsDown = Steps;

switch Mode
    case 'SeveralPerturbations'
        Perturbations = length(ModeOpts.Perts);
        SaveFileName = sprintf([FileName 'Results_' Mode '_%dModels_%s'], EnsembleSize,datestr(now,'yyyymmdd_HHMMSS'));
        if ModeOpts.Verbose
            fprintf(1,'Output file name is: %s.mat\n',SaveFileName);
        end
        ModelResults = repmat([{NaN(size(S,1),Perturbations)} {NaN(size(S,2),Perturbations)} {NaN(size(S,2),Perturbations)} {NaN(Perturbations,1)} {NaN(1,Perturbations)}], EnsembleSize, 1);
 case 'SearchRobust'
        Perturbations = length(ModeOpts.Perts);
        SaveFileName = sprintf([FileName 'Results_' Mode '_%dModels_%s'], EnsembleSize,datestr(now,'yyyymmdd_HHMMSS'));
        if ModeOpts.Verbose
            fprintf(1,'Output file name is: %s.mat\n',SaveFileName);
        end
        ModelResults = repmat([{NaN(size(S,1),Perturbations)} {NaN(size(S,2),Perturbations)} {NaN(size(S,2),Perturbations)} {NaN(Perturbations,1)}], EnsembleSize, 1);
end

Xini = ones(length(X),1);
Uini = (rVnet~=0)+0;
EnsembleKvec = NaN(length(KVEC),EnsembleSize);

if ~isempty(ModeOpts.ImportParameters),
    if strcmp(ModeOpts.ImportParameters(end-3:end),'.mat')
        KvecTemp = load(ModeOpts.ImportParameters);
        EnsembleKvec = KvecTemp.EnsembleKvec;
    elseif strcmp(ModeOpts.ImportParameters(end-3:end),'.csv')
        EnsembleKvec = csvread(ModeOpts.ImportParameters);
    else
        error('Unknown Format')
    end
elseif ~strcmp(Mode,'SearchRobust')
    count=0;
    for Model = 1:EnsembleSize,
        if ModeOpts.Verbose,
            fprintf(1,'SamplingModel: %03d\n',Model);
        end

        Stable = 0;
        
        while ~Stable,
%             Kvec = 10.^((log10(ParamRange(:,1))-log10(ParamRange(:,2))).*rand(length(KVEC), 1)+log10(ParamRange(:,2)));
            Kvec = (((ParamRange(:,1))-(ParamRange(:,2))).*rand(length(KVEC), 1)+(ParamRange(:,2)));
            count=count+1;
            Kvec(ParamInfo(:,1)) = K1S(Xini, Kvec, 1, rVnet, Uini);
            Stable = max(real(eig(JACOBIAN(0, Xini, Kvec, 1,Uini)))) < -1e-6;
        end
        EnsembleKvec(:,Model) = Kvec;
    end
    sprintf('Fraction of stable models = %f', EnsembleSize/count)
%     pause
end

delete(gcp)
parpool(10)

parfor Model = 1:EnsembleSize
    if ModeOpts.Verbose,
        fprintf(1,'AnalyzingModel: %03d\n',Model);
    end
    Kvec = EnsembleKvec(:,Model);
    switch Mode
        case 'SeveralPerturbations'
            [EnzymeConc EnzymeActs EnzymeFluxes Bifs EnzymeTimes] = SeveralPerturbations(ModelFile,Kvec,ModeOpts);
            ModelResults(Model,:) = [{EnzymeConc} {EnzymeActs} {EnzymeFluxes} {Bifs} {EnzymeTimes}];
        case 'SearchRobust'
            [Kvec EnzymeConc EnzymeActs EnzymeFluxes Bifs] = SearchRobust(ModelFile,Kvec,ModeOpts);
            EnsembleKvec(:,Model) = Kvec;
            ModelResults(Model,:) = [{EnzymeConc} {EnzymeActs} {EnzymeFluxes} {Bifs}];    
    end
end

save(sprintf([FileName '_Sampling_%dModels'],EnsembleSize),'EnsembleKvec');
save(SaveFileName,'EnsembleKvec','ModelResults','ModeOpts','ToRemove');

end