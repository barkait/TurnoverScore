clear all;close all;clc;

%% IMPORT GTEX
GTEX_PATH = "X:/Common/Joint_projects/Bilomics/datasets/GTEX/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.txt";
T = readtable(GTEX_PATH);
clear g
g.gene_name=T.Description;
g.tissue=T.Properties.VariableNames(3:end);
for i=1:length(g.tissue)
    g.tissue{i}(findstr(g.tissue{i},'_'))='-';
end
g.mat=table2array(T(:,3:end));


% identify duplicate genes
[uniqueA, ~, idx] = unique(g.gene_name);
counts = histcounts(idx, 1:max(idx)+1);
duplicates = uniqueA(counts > 1);
% stack them
tmp = zeros(length(duplicates), length(g.tissue));
for ii = 1:length(duplicates)
    tmp(ii,:) = sum(g.mat(find(strcmpi(g.gene_name, duplicates(ii))),:));
end
% delete old duplicates
for ii = 1:length(duplicates)
    dup_idx = find(strcmpi(g.gene_name, duplicates(ii)));
    g.mat(dup_idx,:) = [];
    g.gene_name(dup_idx) = [];
end
% attaching the sum of the duplicate values at the end of the table
g.gene_name = [g.gene_name; duplicates];
g.mat = [g.mat; tmp];

% create mat norm 
g.mat_norm=g.mat./sum(g.mat);
% create table
g.df = array2table(g.mat_norm,'RowNames',g.gene_name,'VariableNames',g.tissue);

%% Load Jabri
load("X:\Common\Joint_projects\Bilomics\datasets\Jabri_celiac\jabri_struct.mat")
ind_j = find(j.celiac_status == 0);
ind_j = intersect(ind_j, find(j.age > 30));
j.df = array2table(j.mat_norm(:,ind_j));
j.df.Properties.VariableNames = j.sample_names(ind_j);
j.df.Properties.RowNames = j.gene_name;
j.mean_df = median(j.df,2);
j.mean_df.Properties.VariableNames = "duodenum";
j.mean_df.gene_name = j.gene_name;

%writetable(j.df,"jabri_controls_over_30.csv",'WriteRowNames', true);

%% create signature matrix of UGI organs - GTEX+JABRI
my_gtex_tissues = ["Stomach", "Esophagus-Mucosa"];
%g.df(:,my_gtex_tissues)
tissue_bulk_gtex = g.df(:,my_gtex_tissues);
tissue_bulk_gtex.gene_name = tissue_bulk_gtex.Properties.RowNames;
tissue_bulk_gtex = innerjoin(tissue_bulk_gtex, j.mean_df, 'Keys',{'gene_name','gene_name'});
tissue_bulk_gtex.Properties.RowNames = tissue_bulk_gtex.gene_name;
tissue_bulk_gtex.gene_name = [];

tissue_bulk_gtex.Properties.VariableNames = lower(tissue_bulk_gtex.Properties.VariableNames);
tissue_bulk_gtex.Properties.VariableNames(2) = "esophagus";
tissue_bulk_gtex = tissue_bulk_gtex(:, ["esophagus","stomach","duodenum"]);
writetable(tissue_bulk_gtex,"ugi_sig_mat_gtex_jabri.csv",'WriteRowNames', true);
