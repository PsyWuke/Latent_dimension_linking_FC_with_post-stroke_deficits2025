%% data preparation
rootpath = '/mnt/data8/wuke/project/Corbetta/Index/NodeBasedFC_Gordon333';
stage = 'acute';
out_dir = ['/mnt/data8/wuke/project/Corbetta/PaperResults/part1_PLSC_results' filesep stage];
if strcmp(stage,'acute')
    FC_file = g_ls([rootpath,'/z*_A_*.mat']);
elseif strcmp(stage,'subacute')
    FC_file = g_ls([rootpath,'/z*_C_*.mat']);
elseif strcmp(stage,'chronic')
    FC_file = g_ls([rootpath,'/z*_C2_*.mat']);   
end
X = [];
for i = 1:length(FC_file)
    load(FC_file{i});
    FC_vector = jUpperTriMatToVec(Z,1)';
    X = [X;FC_vector];
end

load('/mnt/data8/wuke/project/Corbetta/Behaviour/corbetta27.2.mat');
domain = fieldnames(corbetta.c128_behavior);
for i=1:length(domain)
    behav_name = ['corbetta.c128_behavior.',domain{i}];
    behav = eval(behav_name);
    if strcmp(domain{i},'attention')
        rawname = fieldnames(behav);
        rawname = rawname([5,6,7,8,13,14,15,16,38,28],1);
        attention = [];
        for j = 1:length(rawname)
            rawdata = eval([behav_name,'.',rawname{j},'.',stage]);
            attention = [attention,rawdata];
        end
        T_atten = array2table(attention,'VariableNames',rawname);
    elseif strcmp(domain{i},'language')
        rawname = fieldnames(behav);
        language = [];
        for j = 1:length(rawname)
            rawdata = eval([behav_name,'.',rawname{j},'.',stage]);
            language = [language,rawdata];
        end
        T_lan = array2table(language,'VariableNames',rawname);
    elseif strcmp(domain{i},'motor')
        rawname = fieldnames(behav);
        rawname = rawname([1:8,17,18,21,28,29,30,33],1);
        motor = [];
        for j = 1:length(rawname)
            rawdata = eval([behav_name,'.',rawname{j},'.chronic']);
            motor = [motor,rawdata];
        end
        T_mot = array2table(motor,'VariableNames',rawname);
    elseif strcmp(domain{i},'memory')
        rawname = fieldnames(behav);
        rawname = rawname([6:7,15:17,19,20],1);
        memory = [];
        for j = 1:length(rawname)
            rawdata = eval([behav_name,'.',rawname{j},'.',stage]);
            memory = [memory,rawdata];
        end
        T_mem = array2table(memory,'VariableNames',rawname);
    end
end
Y = [T_atten,T_lan,T_mot,T_mem];
BehavName = T.Properties.VariableNames;

%% regression
confounds = {'age' 'gender'  'ethnicity' 'smoke' };
regressors=SubInfo{:,confounds};
% Regress out confounds
[X_reg, ~, ~, ~] = glm_regress_matrix(X, regressors, 0, []);
[Y_reg, ~, ~, ~] = glm_regress_matrix(Y, regressors, 0, []);

%% run PLS
disp('Running PLS');
[U,S,V,Lx,Ly,explCovLC,LC_behav_loadings,LC_RSFC_loadings] = ...
    myPLS_analysis(X_reg,Y_reg,1,1);

%% Permutation testing
disp(' Permutation testing over LCs');
nPerms_rest = 1000;
normalization_img = 1; % normalization options for RSFC data
normalization_behav = 1; % normalization options for behavior adta
[pvals_LC,permsamp] = WK_PLS_permut(X_reg,Y_reg,U,S,nPerms_rest,normalization_img,normalization_behav,1000);
% FDR correction over the first 5 LCs
[signif_LC, ~] = FDR(pvals_LC(1:5), 0.05);
disp(['Significant LCs (after FDR correction):' signif_LC]);

%% Bootstrap RSFC & behavior loadings
nBootstraps = 1000; 
signif_LC = 1;
nRois = 333;
out_dir=['/mnt/data8/wuke/project/Corbetta/PaperResults/part1_PLSC_results' filesep stage filesep 'BootStrap_test'];
if ~exist(out_dir)
    mkdir(out_dir);
end
% Re-compute loadings using bootstrap 
[LC_RSFC_loadings_boot,LC_behav_loadings_boot,all_boot_orders] = bootstrap_loadings(X_reg,Y_reg,U,1,nBootstraps,1,1,1000);
save(fullfile(out_dir,[stage '_PLS_bootstrapLoadings_' num2str(nBootstraps) 'bootstraps.mat']),...
    'LC_RSFC_loadings_boot','LC_behav_loadings_boot','X_reg','Y_reg','U','signif_LC','nBootstraps','all_boot_orders');
% Compute confidence intervals, z-scores, p-values of loadings
[std_behav_boot,zscore_behav_boot,pvals_behav_boot,std_RSFC_boot,zscore_RSFC_boot,pvals_RSFC_boot] = bootstrap_stats...
    (behav_loadings,LC_behav_loadings_boot,RSFC_loadings,LC_RSFC_loadings_boot,nRois,1,out_dir,stage);
save(fullfile(out_dir,[stage '_PLS_bootstrapResults_' num2str(nBootstraps) 'bootstraps.mat']),...
    'std_behav_boot','zscore_behav_boot','pvals_behav_boot','std_RSFC_boot','zscore_RSFC_boot','pvals_RSFC_boot',...
    'behav_loadings','LC_behav_loadings_boot','RSFC_loadings','LC_RSFC_loadings_boot','nRois','signif_LC','out_dir');
