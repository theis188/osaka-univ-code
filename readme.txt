+ First step to recreate simulations: run the simulation with ModelCOA.mat and ModelNoCOA.mat.
    -Open 'nButanolModeOptsPreparation60Final.m'
    -Change the first argument in the function on line 58 to 'ModelCOA' or 'ModelNoCOA' as desired.
    -Run for each condition
    -The output will be a result file of the form [ModelName]Results_SeveralPerturbations_500Models_yyyymmdd_hhmmss.mat

+ Figure A): 
    -Open Comparison.m
    -On line 14, edit strain variable to take on the desired quantity.
    -On line 6, use the desired result filename (usually ModelCOA)
    -Uncomment lines 50-53 (make sure lines 57-60 are commented)

+ Figure B)
    -Open Comparison.m
    -On line 14, set strain=3
    -On line 6, set result filename to ModelNoCOA results.
    -Uncomment lines 57-60 (make sure lines 50-53 are commented).

+ Figure C)
    -Extract data from desired result filenames (plotted in Excel).
    -In ModelCOA result files, ModelResults{1:500,1}(44,1:5) will be COA-fold change amounts.

+ Figure D): 
    -Open Comparison.m
    -On line 14, edit strain variable to take on the desired quantity.
    -On line 6, use ModelCOA results
    -Uncomment lines 50-53 (make sure lines 57-60 are commented)