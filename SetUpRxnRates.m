function [ V, ParamInfo, SubstrateConc, ParameterRange, rV] = SetUpRxnRates(S, Sreg, KVEC, X, Vin_index, Vout_index, U,rV,Revs,S2a,EnzName,MetabName)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

    MM_1_1_Params = 4;
    MM_1_Irr_Params = 2;
    MM_2_Irr_Params = 3;
    RandomBiBi_Params = 6;
    OrderedBiUni_Params = 5;
    OrderedUniBi_Params = 5;
%     MassAction_Params = 2;
    OutRxn_Params = 1;
    PTS_Params = 5;

    ParametersUsed = 0;
    V = sym('V', [size(S,2) 1]);
    ParamInfo = zeros(size(S,2), 2);
    SubstrateConc = sym('SubstrateConc');
    ParameterRange = zeros([length(KVEC), 2]);
    for Rxn = 1: size(S,2),

        Subs = abs(sum(S(S(:,Rxn)<0, Rxn)));
        Prods = sum(S(S(:,Rxn)>0, Rxn));
        Rev = Revs(Rxn);
        Regs = Sreg(find(Sreg(:,Rxn)>0), Rxn);
        Inlet = ~isempty(find(Vin_index==Rxn, 1));
        if Inlet,
            Subs = Subs+1;
            if strcmp(EnzName(Rxn),'PTS_r'),
                MetabolitesVector = [X(find(strcmp('PHOSPHO-ENOL-PYRUVATE',MetabName))); X(find(strcmp('PYRUVATE',MetabName))); X(find(strcmp('GLC-6-P',MetabName))); SubstrateConc];
            elseif (Subs == 1 || Subs == 2) && (Prods==1 || Prods==2), 
                MetabolitesVector = [SubstrateConc; X(find(S(:,Rxn)<0)); X(find(S(:,Rxn)== -2)); X(find(S(:,Rxn)>0)); X(find(S(:,Rxn)== 2))];
            else
                MetabolitesVector = [SubstrateConc; X(find(S(:,Rxn)<0));  X(find(S(:,Rxn)>0))];
            end
        else
            if strcmp(EnzName(Rxn),'PTS_r'),
                MetabolitesVector = [X(find(strcmp('PHOSPHO-ENOL-PYRUVATE',MetabName))); X(find(strcmp('PYRUVATE',MetabName))); X(find(strcmp('GLC-6-P',MetabName))); SubstrateConc];
            elseif (Subs == 1 || Subs == 2) && (Prods==1 || Prods==2) && isempty(Regs),
                MetabolitesVector = [X(find(S(:,Rxn)<0)); X(find(S(:,Rxn)== -2)); X(find(S(:,Rxn)>0)); X(find(S(:,Rxn)== 2))];
            else
                MetabolitesVector = [X(find(S(:,Rxn)<0));  X(find(S(:,Rxn)>0))];
            end
        end
        
        if strcmp(EnzName(Rxn),'PTS_r'),
            [V(Rxn), ParameterRange(ParametersUsed+[1:PTS_Params], :)] = PTS_Rate(KVEC(ParametersUsed+[1:PTS_Params]), MetabolitesVector, rV(Rxn), U(Rxn));
            ParamInfo(Rxn,:) = ParametersUsed+[1 PTS_Params];
            ParametersUsed = ParametersUsed + PTS_Params;
        elseif ~isempty(find(Vout_index==Rxn, 1)),
%             if length(MetabolitesVector)==1,
%                 [V(Rxn) ParameterRange(ParametersUsed+[1:MassAction_Params],:)] = MassAction(KVEC(ParametersUsed+[1:MassAction_Params]), MetabolitesVector, length(find(S(:,Rxn)<0)), abs([S(find(S(:,Rxn)<0), Rxn);  S(find(S(:,Rxn)>0), Rxn)]), rV(Rxn));
%                 ParamInfo(Rxn,:) = ParametersUsed+[1 MassAction_Params];
%                 ParametersUsed = ParametersUsed + MassAction_Params;
                [V(Rxn) ParameterRange(ParametersUsed+[1:OutRxn_Params], :)] = OutletRxn(KVEC(ParametersUsed+[1:OutRxn_Params]), MetabolitesVector, U(Rxn));
                ParamInfo(Rxn,:) = ParametersUsed+[1 OutRxn_Params];
                ParametersUsed = ParametersUsed + OutRxn_Params;
        elseif Subs ==1 && Prods == 1 && Rev && isempty(Regs),
            [V(Rxn) ParameterRange(ParametersUsed+[1:MM_1_1_Params], :)] = MichaelisMenten_1To1(KVEC(ParametersUsed+[1:MM_1_1_Params]), MetabolitesVector, rV(Rxn), U(Rxn));
            ParamInfo(Rxn,:) = ParametersUsed+[1 MM_1_1_Params];
            ParametersUsed = ParametersUsed + MM_1_1_Params;
        elseif Subs==2 && Prods==2 && Rev && isempty(Regs),
            [V(Rxn) ParameterRange(ParametersUsed+[1:RandomBiBi_Params],:)] = RandomBiBi(KVEC(ParametersUsed+[1:RandomBiBi_Params]), MetabolitesVector, rV(Rxn), U(Rxn));
            ParamInfo(Rxn,:) = ParametersUsed+[1 RandomBiBi_Params];
            ParametersUsed = ParametersUsed + RandomBiBi_Params;
        elseif Subs==2 && Prods==1 && Rev && isempty(Regs),
            [V(Rxn) ParameterRange(ParametersUsed+[1:OrderedBiUni_Params], :)] = OrderedBiUni(KVEC(ParametersUsed+[1:OrderedBiUni_Params]), MetabolitesVector, rV(Rxn), U(Rxn));
            ParamInfo(Rxn,:) = ParametersUsed+[1 OrderedBiUni_Params];
            ParametersUsed = ParametersUsed + OrderedBiUni_Params;
        elseif Subs==1 && Prods==2 && Rev && isempty(Regs),
            [V(Rxn) ParameterRange(ParametersUsed+[1:OrderedUniBi_Params], :)] = OrderedUniBi(KVEC(ParametersUsed+[1:OrderedUniBi_Params]), MetabolitesVector, rV(Rxn), U(Rxn));
            ParamInfo(Rxn,:) = ParametersUsed+[1 OrderedUniBi_Params];
            ParametersUsed = ParametersUsed + OrderedUniBi_Params;
        elseif ~Rev && Subs==1 && isempty(Regs),
            [V(Rxn) ParameterRange(ParametersUsed+[1:MM_1_Irr_Params], :)] = MM_1_Irr(KVEC(ParametersUsed+[1:MM_1_Irr_Params]), MetabolitesVector, rV(Rxn), U(Rxn));
            ParamInfo(Rxn,:) = ParametersUsed+[1 MM_1_Irr_Params];
            ParametersUsed = ParametersUsed + MM_1_Irr_Params;
        elseif ~Rev && Subs==2 && isempty(Regs),
            [V(Rxn) ParameterRange(ParametersUsed+[1:MM_2_Irr_Params], :)] = MM_2_Irr(KVEC(ParametersUsed+[1:MM_2_Irr_Params]), MetabolitesVector, rV(Rxn), U(Rxn));
            ParamInfo(Rxn,:) = ParametersUsed+[1 MM_2_Irr_Params];
            ParametersUsed = ParametersUsed + MM_2_Irr_Params;
        else
            NoOfParams = 2+length(find(S(:,Rxn)<0))+length(find(S(:,Rxn)>0));
            for mm=Regs'
                if mm<3,    %1 and 2 are non allosteric inhibition and activation
                    NoOfParams = NoOfParams+1;
                else
                    NoOfParams = NoOfParams+2;    %3 and 4 indicate allosteric mechanisms
                end
            end        
            [V(Rxn) ParameterRange(ParametersUsed+[1:NoOfParams], :)] = Liebermeister(Rev,KVEC(ParametersUsed+[1:NoOfParams]), [MetabolitesVector; X(find(Sreg(:,Rxn)>0))], Sreg(find(Sreg(:,Rxn)>0),Rxn), length(find(S(:,Rxn)<0)), abs([S(find(S(:,Rxn)<0), Rxn);  S(find(S(:,Rxn)>0), Rxn)]),U(Rxn),rV(Rxn));
            ParamInfo(Rxn,:) = ParametersUsed + [1 NoOfParams];
            ParametersUsed = ParametersUsed + NoOfParams;
%         else
%             [V(Rxn) ParameterRange(ParametersUsed+[1:MassAction_Params])] = MassAction(KVEC(ParametersUsed+[1:MassAction_Params]), MetabolitesVector, length(find(S(:,Rxn)<0)), abs([S(find(S(:,Rxn)<0), Rxn);  S(find(S(:,Rxn)>0), Rxn)]), rV(Rxn));
%             ParamInfo(Rxn,:) = ParametersUsed+[1 MassAction_Params];
%             ParametersUsed = ParametersUsed + MassAction_Params;
        end
%         VRxn = subs(V(Rxn), SubstrateConc, 1);
    end
    for Rxn=1:size(S2a,2)
        V(size(S,2)+Rxn) = U(size(S,2)+Rxn);
    end
end

function [ V, ParameterRange ] = MichaelisMenten_1To1( k, x, rV, U)
    Vm=k(1);    Keq = k(2);    kms = k(3);    kmp = k(4);
    S = x(1); P = x(2);
    V = U*Vm*(S-P/Keq)/kms/(1+S/kms+P/kmp);
    if rV >0,
        ParameterRange = [0 Inf; 1 10; .1 10; .1 10];
    else
        ParameterRange = [0 Inf; .1 1; .1 10; .1 10];
    end
end

function [ V, ParameterRange ] = RandomBiBi( k, x, rV, U)
    Vm = k(1);    Keq = k(2);    kms1 = k(3);    kms2 = k(4); kmp1 = k(5);   kmp2 = k(6);
    S1 = x(1); S2 = x(2); P1 = x(3); P2=x(4);
    V = U*Vm*(S1*S2-P1*P2/Keq)*(1/kms1/kms2)/(1+S1/kms1+S2/kms2+S1*S2/kms1/kms2+P1/kmp1+P2/kmp2+P1*P2/kmp1/kmp2);
    if rV >0,
        ParameterRange = [0 Inf; 1 10; .1 10; .1 10; .1 10; .1 10];
    else
        ParameterRange = [0 Inf; .1 1; .1 10; .1 10; .1 10; .1 10];
    end    
end

function [ V, ParameterRange ] = OrderedBiUni( k, x, rV, U)
    Vm = k(1);    Keq = k(2);    kms1 = k(3);    kms2 = k(4);   kmp1 = k(5);
    S1 = x(1); S2 = x(2); P1=x(3);
    V = U*Vm*(S1*S2-P1/Keq)*(1/kms1/kms2)/(1+S1/kms1+S2/kms2+S1*S2/kms1/kms2+P1/kmp1);
    if rV >0,
        ParameterRange = [0 Inf; 1 10; .1 10; .1 10; .1 10];
    else
        ParameterRange = [0 Inf; .1 1; .1 10; .1 10; .1 10];
    end    
end

function [ V, ParameterRange ] = OrderedUniBi( k, x, rV, U)
    Vm = k(1);    Keq = k(2);    kms = k(3);    kmp1 = k(4);   kmp2 = k(5);
    S = x(1); P1 = x(2); P2=x(3);
    V = U*Vm*(S-P1*P2/Keq)*(1/kms)/(1+S/kms+P1/kmp1+P2/kmp2+P1*P2/kmp1/kmp2);
    if rV >0,
        ParameterRange = [0 Inf; 1 10; .1 10; .1 10; .1 10];
    else
        ParameterRange = [0 Inf; .1 1; .1 10; .1 10; .1 10];
    end    
end

function [ V, ParameterRange ] = MassAction( k, x, NoSubs, Exponents, rV, U)
    V = U*k(1)*prod(x(1:NoSubs).^Exponents(1:NoSubs)) - k(2)*prod(x(NoSubs+1:end).^Exponents(NoSubs+1:end));
    if rV >0,
        ParameterRange = [0 Inf; .01*rV 1*rV];
    else
        ParameterRange = [0 Inf; rV 2*rV];
    end    
end

function [ V, ParameterRange ] = MM_1_Irr( k, x, ~, U)
    Vm = k(1); Km = k(2);
    S = x(1);
    V = U*Vm/(Km/S + 1);
    ParameterRange = [0 Inf; .1 10];
end

function [ V, ParameterRange ] = MM_2_Irr( k, x, ~, U)
    Vm = k(1); Kma = k(2); Kmb = k(3);
    A = x(1); B = x(2);
    V = U*Vm/(Kma*Kmb/A/B + Kmb/B + Kma/A + 1);
    ParameterRange = [0 Inf; .1 10; .1 10];
end


function [ V, ParameterRange ] = OutletRxn( k, x, U)
    V = U*k(1)*prod(x);
    ParameterRange = [0 Inf];
end

function [ V, ParameterRange ] = PTS_Rate( k, x, rV, U)
    Vm=k(1);    k1 = k(2);    k2 = k(3);    k3 = k(4); ki = k(5);
    PEP = x(1); PYR = x(2); G6P = x(3); Glc = x(4);
    
    V = U*Vm/((k1*PYR/Glc/PEP + k2/Glc + k3*PYR/PEP + 1)*(1+G6P/ki));
    ParameterRange = [0 Inf; .1 10; .1 10; .1 10; .1 10];
end

function [V, ParameterRange] = Liebermeister(Rev,k, x, RegIndices, NoSubs, Exponents, U, rV)
    % k1 = Vf   k2 = keq    k3:2+kNoSubs = kmsub   k2+NoSubs+1:2+NoSubs+NoProds = kmprod     kmNoProds+1:kmNoProds+NoRegs = [kmM PM ki]
    NoRegs = length(RegIndices);
    NoProds = length(x)-NoSubs-NoRegs;
%     T = prod((x(1:NoSubs)./k(3:2+NoSubs)).^Exponents(1:NoSubs))*(1-prod(x(NoSubs+1:NoSubs+NoProds))./prod(x(1:NoSubs))./k(2));
%     D = prod((1+(x(1:NoSubs)./k(3:2+NoSubs))).^Exponents(1:NoSubs))+prod((1+(x(NoSubs+1:NoSubs+NoProds)./k(2+NoSubs+1:2+NoSubs+NoProds))).^Exponents(NoSubs+1:NoSubs+NoProds))-1;
    T = prod((x(1:NoSubs)./k(3:2+NoSubs)))*(1-prod(x(NoSubs+1:NoSubs+NoProds))./prod(x(1:NoSubs))./k(2));
    D = prod((1+(x(1:NoSubs)./k(3:2+NoSubs))))+prod((1+(x(NoSubs+1:NoSubs+NoProds)./k(2+NoSubs+1:2+NoSubs+NoProds))))-1;
    Rreg=1;
    Dreg = 0;
    CurrentK = 2+NoSubs+NoProds+1;
    if rV >0 && Rev,
        ParameterRange(1:2, :) = [0 Inf; 1 10];
    elseif rV>0 && ~Rev
        ParameterRange(1:2, :) = [0 Inf; 100 200];
    else
        ParameterRange(1:2, :) = [0 Inf; .1 1];
    end    
    ParameterRange(3:2+NoSubs+NoProds,:) = repmat([.1 10], NoSubs+NoProds, 1);
    for m=1:length(RegIndices),
        switch RegIndices(m)
            case 1
                %Competitive Inhibition
                Dreg = Dreg + x(NoSubs+NoProds+m)/k(CurrentK);
                ParameterRange(CurrentK, :) = [.1 10];
                CurrentK = CurrentK + 1;
            case 2
                %Non-Allosteric Activation
                Dreg = Dreg + k(CurrentK)/x(NoSubs+NoProds+m);
                ParameterRange(CurrentK, :) = [.1 10];
                CurrentK = CurrentK+1;
            case 3
                %Allosteric Inhibition
                Rreg = Rreg * (k(CurrentK)+(1-k(CurrentK))/(1+x(NoSubs+NoProds+m)/k(CurrentK+1)));
                ParameterRange(CurrentK:CurrentK+1, :) = [0.0000001 0.0000001; .001 10];
                CurrentK = CurrentK+2;
            case 4
                %Allosteric Activation
                Rreg = Rreg * (k(CurrentK)+(1-k(CurrentK))*x(NoSubs+NoProds+m)/k(CurrentK+1)/(1+x(NoSubs+NoProds+m)/k(CurrentK+1)));
                ParameterRange(CurrentK:CurrentK+1, :) = [0.0000001 1; .1 10];
                CurrentK = CurrentK+2;
        end
    end
    
    V = simplify(U*k(1) *Rreg * T/(D+Dreg));

end

% % % % function [V, ParameterRange] = LiebermeisterOutlet(k, x, RegIndices, NoSubs, Exponents, U)
% % % %     % k1 = Vf   k2:1+kNoSubs = kmsub   kmNoSubs+2:kmNoSubs+NoRegs+2 = [kmM PM ki]
% % % % %     T = prod((x(1:NoSubs)./k(2:1+NoSubs)).^Exponents(1:NoSubs));
% % % % %     D = prod((1+(x(1:NoSubs)./k(2:1+NoSubs))).^Exponents(1:NoSubs))-1;
% % % %     T = prod((x(1:NoSubs)./k(2:1+NoSubs)));
% % % %     D = prod((1+(x(1:NoSubs)./k(2:1+NoSubs))))-1;
% % % %     Rreg=1;
% % % %     Dreg = 0;
% % % %     CurrentK = 2+NoSubs;
% % % %     ParameterRange(1,:) = [0 Inf];
% % % %     ParameterRange(2:1+NoSubs,:) = repmat([0 100.00], NoSubs, 1);
% % % %     for m=1:length(RegIndices),
% % % %         switch RegIndices(m)
% % % %             case 1
% % % %                 Dreg = Dreg + x(NoSubs+m)/k(CurrentK);
% % % %                 ParameterRange(CurrentK, :) = [0 100.00];
% % % %                 CurrentK = CurrentK + 1;
% % % %             case 2
% % % %                 Dreg = Dreg + k(CurrentK)/x(NoSubs+m);
% % % %                 ParameterRange(CurrentK, :) = [0 100.00];
% % % %                 CurrentK = CurrentK + 1;
% % % %             case 3
% % % %                 Rreg = Rreg * (k(CurrentK)+(1-k(CurrentK))/(1+x(NoSubs+m)/k(CurrentK+1)));
% % % %                 ParameterRange(CurrentK:CurrentK+1, :) = [.9 1; 0 100.00];
% % % %                 CurrentK = CurrentK + 1;
% % % %             case 4
% % % %                 Rreg = Rreg * (k(CurrentK)+(1-k(CurrentK))*x(NoSubs+m)/k(CurrentK+1)/(1+x(NoSubs+m)/k(CurrentK+1)));
% % % %                 ParameterRange(CurrentK:CurrentK+1, :) = [0 1; 0 100.00];
% % % %                 CurrentK = CurrentK + 1;
% % % %         end
% % % %     end
% % % %     
% % % %     V = simplify(U*k(1) *Rreg * T/(D+Dreg));
% % % % 
% % % % end