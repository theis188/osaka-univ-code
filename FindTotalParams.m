function [ TotalParams ] = FindTotalParams( S, Sreg, Vin_index, Vout_index,Revs,EnzName )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    MM_1_1_Params = 4;
    MM_1_Irr_Params = 2;
    MM_2_Irr_Params = 3;
    RandomBiBi_Params = 6;
    OrderedBiUni_Params = 5;
    OrderedUniBi_Params = 5;
    MassAction_Params = 2;
    OutRxn_Params = 1;
    PTS_Params = 5;

TotalParams = 0;
for m = 1:size(S,2),
    Subs = abs(sum(S(S(:,m)<0, m)));
    Prods = sum(S(S(:,m)>0, m));
    Regs = find(Sreg(:,m)>0);
    Rev = Revs(m);
    Inlet = ~isempty(find(Vin_index==m, 1));
    if Inlet,
        Subs = Subs+1;
    end

    if strcmp(EnzName(m),'PTS_r'),
        TotalParams = TotalParams + PTS_Params;
    elseif ~isempty(find(Vout_index==m, 1)),
%         if length(find(S(:,m)<0))==1,
            TotalParams = TotalParams + OutRxn_Params;
%         elseif length(find(S(:,m)<0))>10,
% %             TotalParams = TotalParams;
%         else
%             TotalParams = TotalParams+ 1+length(find(S(:,m)<0));
%             for mm=Regs'
%                 if mm<3,    %1 and 2 are non allosteric inhibition and activation
%                     TotalParams = TotalParams+1;
%                 else
%                     TotalParams = TotalParams+2;    %3 and 4 indicate allosteric mechanisms
%                 end
%             end        
%         end
    elseif Subs ==1 && Prods == 1 && Rev && isempty(Regs),
        TotalParams = TotalParams + MM_1_1_Params;
    elseif Subs==2 && Prods==2 && Rev && isempty(Regs),
        TotalParams = TotalParams + RandomBiBi_Params;
    elseif Subs==2 && Prods==1 && Rev && isempty(Regs),
        TotalParams = TotalParams + OrderedBiUni_Params;
    elseif Subs==1 && Prods==2 && Rev && isempty(Regs),
        TotalParams = TotalParams + OrderedUniBi_Params;
    elseif ~Rev && Subs==1 && isempty(Regs),
        TotalParams = TotalParams + MM_1_Irr_Params;
    elseif ~Rev && Subs==2 && isempty(Regs),
        TotalParams = TotalParams + MM_2_Irr_Params;
    else
        TotalParams = TotalParams+ 2+length(find(S(:,m)<0))+length(find(S(:,m)>0));
        for mm=Regs'
            if Sreg(mm,m)<3,    %1 and 2 are non allosteric inhibition and activation
                TotalParams = TotalParams+1;
            else
                TotalParams = TotalParams+2;    %3 and 4 indicate allosteric mechanisms
            end
        end        
%     else
%         TotalParams = TotalParams + MassAction_Params;
    end

end
