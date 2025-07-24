close all hidden;clear all;clc %#ok<*CLALL>
% load data
%PROJ_BASE = "X:\Common\Joint_projects\Bilomics\Analysis_202503\bilomics_dup\01-shed_preop_umi\";
PROJ_BASE = pwd()+"/01-shed_preop_umi_with_dups/";
obs_path = strcat(PROJ_BASE, "obs.csv");
var_path = strcat(PROJ_BASE, "var.csv");
seq_data_path = strcat(PROJ_BASE, "seq_data.csv");

% create struct
clear t
t.seq_data = readtable(PROJ_BASE + "seq_data.csv");
t.bc = t.seq_data.Properties.VariableNames;
t.gene_name = t.seq_data.gene_name;
t.seq_data = table2array(t.seq_data(:,2:end));
t.mat_norm = t.seq_data./sum(t.seq_data);

% load some fields from OBS.csv
obs = readtable(PROJ_BASE + "obs.csv",'Format','auto');
t.idx = obs.Var1;
field_to_keep = ["sample_name","notes","drain","cancer","time","batch","age","bmi","sex","diagnosis","surgeryType","meds","pmh","psh","id"];
for ii = 1:length(field_to_keep)
    t.(field_to_keep(ii)) = obs.(field_to_keep(ii));
end

% validate data proecssing (in that specific case - its already processed)
t.sample_name_print = strrep(t.sample_name,'_','-');
t = process_structure(t);

%% identify duplicated IDs
dup_ids = [];
for ii = 1:length(t.id)
    my_id = t.id(ii);
    if sum(t.id == my_id) > 1
        dup_ids = [dup_ids, my_id];
    end
end
dup_ids = unique(dup_ids);

t.identify_id = t.id;
for ii = 1:length(dup_ids)
    my_id = dup_ids(ii);
    my_id_idx = find(t.identify_id == my_id);
    t.identify_id(my_id_idx) = t.identify_id(my_id_idx);
end
%% Clustergram
EXP_THRESH = 1e-4;
SZ = 25;
PN= 1e-6;
ind_samples=1:length(t.cancer);%find(strcmpi(t.cancer,'Control'));

ind_exp=find(median(t.mat_norm,2)>EXP_THRESH);

mat = t.mat_norm(ind_exp,:);
mat=log10(mat+EXP_THRESH);
%Z=(mat-mean(mat,2))./std(mat,[],2);
Z = zscore(mat')';


x = cellstr(num2str(t.identify_id));%= t.sample_name;%
y = repmat({'black'}, length(x), 1);  % All columns green
y(find(t.identify_id == 146))={'red'};
y(find(t.identify_id == 125))={'green'};
y(find(t.identify_id == 120))={'cyan'};
col_color_struct = struct('Labels', x, 'Colors', y);


pdist_type = 'euclidean';
cgo = clustergram(Z, 'Rowlabels', t.gene_name(ind_exp), ...
    'columnlabels', t.identify_id, ...
    'ColumnPDist', pdist_type, 'RowPDist', pdist_type, ...
    'Colormap', redbluecmap, ...
    'ColumnLabelsColor', col_color_struct, ...
    'LabelsWithMarkers', true);