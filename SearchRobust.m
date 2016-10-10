function [Kvec EnzymeConc EnzymeActs EnzymeFluxes Bifs] = SearchRobust(Model,Kvec,ModeOpts)

    load(Model)
    Unrobust = 1;
    Xini = ones(length(X),1);
    Uini = (rVnet~=0)+0;
    StartTime = clock;
    Perturbations = length(ModeOpts.Perts);
    ProblemCount = ModeOpts.StartProblemCount;
    nRegularRxns = length(find(Uini>0));
    while Unrobust
        Stable = 0;
        while ~Stable,
%             Kvec = 10.^((log10(ParamRange(:,1))-log10(ParamRange(:,2))).*rand(length(KVEC), 1)+log10(ParamRange(:,2)));
            Kvec = (((ParamRange(:,1))-(ParamRange(:,2))).*rand(length(KVEC), 1)+(ParamRange(:,2)));
            Kvec(ParamInfo(:,1)) = K1S(Xini, Kvec, 1, rVnet,Uini);
            Stable = max(real(eig(JACOBIAN(0, Xini, Kvec, 1, Uini)))) < -1e-6;
        end

        EnzymeConc = NaN(size(S,1),Perturbations);
        EnzymeActs = NaN(size(S,2),Perturbations);
        EnzymeFluxes = NaN(size(S,2),Perturbations);
        Bifs = NaN(Perturbations,1);
        
        BadModel = 0;
        
        ind=1;
        [~,IX] = sort(ProblemCount,'descend');
        while ind<=Perturbations && ~BadModel,
            n = IX(ind);
            U = Uini;
            X = Xini;
            Bif = 0;
            nn=1;
            while nn<=size(ModeOpts.Perts{n},2) && ~Bif
                Uf1 = ModeOpts.Perts{n}(:,nn);
                [uf, conc, Bif] = SimpleODESolverMatlab(1, X, U, Uf1, Kvec,DVDX, DVDU,S );
                U = uf(:,end);
                X = conc(end,:)';
                nn= nn+1;
            end
            EnzymeConc(:,n) = conc(end,:)';
            EnzymeActs(:,n) = uf(:,end);
            EnzymeFluxes(:,n) = FLUXES(conc(end,:)',Kvec,1,uf(:,end));
            Bifs(n) = Bif;
            if max(max(uf(1:nRegularRxns,end)), 1/min(uf(1:nRegularRxns,end)))<2,
                BadModel = 1;
                ProblemCount(n) = ProblemCount(n)+1;
            end
            ind=ind+1;
        end                
        
        NoPerts = size(ModeOpts.Perts,2)/2;
        NonrobustEnzUp = diag(EnzymeActs(1:NoPerts, 1:NoPerts));
        NonrobustEnzDown = diag(EnzymeActs(1:NoPerts, NoPerts+1:end));
        mu = 0;
        sd = .5;
        p = logncdf(NonrobustEnzUp, mu, sd) - logncdf(NonrobustEnzDown, mu, sd);
        Entropy =  -sum(p.*log10(p));
        if Entropy<.3,
            Unrobust = 0;
        end
    end    
    csvFile = ['RobustModels_' Model '.csv'];
    if exist(csvFile, 'file'),
        Kvecs = csvread(csvFile);
        Kvecs = [Kvecs Kvec];
        csvwrite(csvFile,Kvecs)
    else
        csvwrite(csvFile,Kvec)
    end
end