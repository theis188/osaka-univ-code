clear all

COA299FScreened=[0.336448944896428,0.0279460373438210,0.0339914789447480,0.0343654773553052,0.00161841666820346,0.00355541939728973,0.250114336342905,0.668159982733697]
COA299FStdScreened=[0.199948313704922,0.0167607157367021,0.0374948003973567,0.0203089049412471,0.00761851806209613,0.00152580154785969,0.234083535149447,0.125650578554459]

load('ModelNoCOAResults.mat')

load('ModelCOA.mat')

load('Yields_observed.mat')

metaboliteList = {'Pyruvate_out','Succinate_out','Lactate_out','Formate_out','Acetate_out','Butyrate_out','EtOH_out','Butanol_out'};
ScreenedModels =[];
strain = 3 %1 = JCL16F; 2 = JCL166F; 3 = JCL299F
for i = 1:500;
    if ModelResults{i,3}(61,3)>10;
        ScreenedModels = [ScreenedModels i];
    end
end
if strain ~= 3
    ScreenedModels = 1:500 %%%Uncomment this line for UNSCREENED PREDICTIONS
end
COA = [];
for h = 1:size(ScreenedModels,2);
for i = 1:8
for z = 1:5
    j = ScreenedModels(1,h);
    MetaboliteOut = ModelResults{j,3}(find(strcmp(metaboliteList{i},Net.EnzName)), z);
    GlucCons = ModelResults{j,3}(find(strcmp('PTS_r',Net.EnzName)), z) ...
        + ModelResults{j,3}(find(strcmp('glk',Net.EnzName)), z);
    PredictedYield(z,i,h) = MetaboliteOut / GlucCons;
    COA(h)=ModelResults{j,1}(44,1);
end
end
end

for i = 1:8
    VrefMetabFlux(i) = Net.Vref(find(strcmp(metaboliteList{i},Net.EnzName)));
end
VrefGluc = Net.Vref(find(strcmp('PTS_r',Net.EnzName))) ...
        + Net.Vref(find(strcmp('glk',Net.EnzName)));
VrefYield = VrefMetabFlux / VrefGluc

StdYield = real(std(PredictedYield,[],3));
AvgPredYield = real(mean(PredictedYield,3));
hold on

%% Uncomment this block for figs A & D
 
% x = [1 2 3 4 5 6 7 8] + 0.16;
% bar(permute([MetabYield(strain+1,:) ;AvgPredYield(strain,:)], [2 1]));
% errorbar(x, AvgPredYield(strain,:), StdYield(strain,:), 'k', 'linestyle', 'none');
% set(gca,'XTick',1:8,'XTickLabel',{'PYRUVATE','SUC','LACTATE','FORMATE','ACET','BUTYRATE','ETOH','BUOH'});

%% Uncomment this block for fig B

% xx = [1 2] 
% bar(permute([MetabYield(4,[1 5]); COA299FScreened([1 5]); AvgPredYield(3,[1 5])], [2 1]));
% errorbar([xx xx+0.22], [COA299FScreened([1 5]) AvgPredYield(3,[1 5])], [COA299FStdScreened([1 5]) StdYield(3,[1 5])], 'k', 'linestyle', 'none');
% set(gca,'XTick',1:2,'XTickLabel',{'PYRUVATE','ACET'});

%%

ylim([0 1])
ax = gca
ax.XTickLabelRotation = 60
set(gcf,'color','w');

ylabel('Yield, (mol/mol glc consumed)')
xlabel('Metabolite')
