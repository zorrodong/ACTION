clear
addpath(genpath('../../modules'));

%% Read RNASeq dataset (subset of the Mouse sci-ATAC-seq Atlas)
input_path = '../../input/datasets/MouseBrain/RNA';
fname = fullfile(input_path, 'preprocessed.mat');
if(~exist(fname, 'file'))
    RNASeq.expression = mmread(fullfile(input_path, 'expression.mm'));
    RNASeq.gene_names = readtable(fullfile(input_path, 'gene_names.txt'), 'Delimiter', '\t', 'FileType', 'text');
    RNASeq.sample_annotations = readtable(fullfile(input_path, 'sample_annotations.txt'), 'Delimiter', '\t', 'FileType', 'text');
    save(fname, 'RNASeq');
else
    load(fname);
end
Labels = RNASeq.sample_annotations.Labels;
UL = unique(Labels);
[~, labels] = ismember(Labels, UL);

%% Reduce gene expression profile
PCA_dim = 30;
method_no = 1;
tic; [S_r, ~, explained_var] = reduceGeneExpression(RNASeq.expression, PCA_dim, method_no); toc

plot(explained_var); xlabel('# components'); ylabel('Explained variance');

S = S_r'*S_r;
plotKernel_HeatMap(S, Labels);
drawnow();
%% Run ACTION
[C_trace, H_trace] = runACTION(S_r, 2, 20, 8);        

ACTION_ARI = zeros(20, 1);
ACTION_NMI = zeros(20, 1);

for i = 2:20
    H = H_trace{i};

    [~, predicted_celltyps] = max(H);
    ACTION_ARI(i) = adjustedrand(labels, predicted_celltyps);
    ACTION_NMI(i) = nmi(labels, predicted_celltyps);
    fprintf('ACTION:: k = %d, NMI = %.2f, ARI = %.2f\n', i, ACTION_NMI(i), ACTION_ARI(i));
end

figure; plot(2:20, ACTION_ARI(2:20)); xlabel('# archs'); ylabel('ARI');
figure; plot(2:20, ACTION_NMI(2:20)); xlabel('# archs'); ylabel('NMI');

k = numel(UL);
H = H_trace{k};
[~, predicted_celltyps] = max(H);

%% Build ACTIONet
    build