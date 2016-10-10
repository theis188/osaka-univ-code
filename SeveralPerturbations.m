function [EnzymeConc, EnzymeActs, EnzymeFluxes, Bifs, EnzymeTimes] = SeveralPerturbations(Model,Kvec,ModeOpts)

    load(Model)
    Xini = ones(length(X),1);
    Uini = (rVnet~=0)+0;
        
    Perturbations = length(ModeOpts.Perts);
    
    EnzymeConc = NaN(size(S,1),Perturbations);
    EnzymeActs = NaN(size(S,2),Perturbations);
    EnzymeFluxes = NaN(size(S,2),Perturbations);
    EnzymeTimes = NaN(1,Perturbations);
    Bifs = NaN(Perturbations,1);
    
    for n=1:Perturbations,
        U = Uini;
        X = Xini;
        Bif = 0;
        nn=1;
        while nn<=size(ModeOpts.Perts{n},2) && ~Bif
            Uf1 = ModeOpts.Perts{n}(:,nn);
            [uf, conc, Bif, TimeTaken] = SimpleODESolverMatlab(1, X, U, Uf1, Kvec,DVDX, DVDU,S );
            U = uf(:,end);
            X = conc(end,:)';
            nn= nn+1;
        end
        
        EnzymeConc(:,n) = conc(end,:)';
        EnzymeActs(:,n) = uf(:,end);
        EnzymeFluxes(:,n) = FLUXES(conc(end,:)',Kvec,1,uf(:,end));
        EnzymeTimes(n) = TimeTaken;
        Bifs(n) = Bif;
        
    end        
   
    
end