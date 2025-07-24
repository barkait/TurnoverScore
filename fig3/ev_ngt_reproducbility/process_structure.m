function t=process_structure(t,only_mito_flag)

% This function receives the mcSCRBseq UMI table, removes mitochondrial
% genes and non-coding genes

if nargin<2
    only_mito_flag=0;
end
%% Sort the exp. values in each sample and show the fraction taken by the top gene
[y,ord]=sort(t.mat_norm);
figure;
bar(y(end,:));
% add ticks from 1 to 30 to the current axes
set(gca,'xtick',1:length(t.sample_name));
% add the names of the samples to the x axis
set(gca,'XTickLabel',t.sample_name_print)
% rotate sample names 45 degrees
set(gca,'XTickLabelRotation',45);
ylabel('Fraction of total UMI');
title("MAX fraction of top gene")
%set(gcf,'position',[1000         785        1432         553]);
%set(gcf,'position',[869         173        1432         644]);

% print on screen the most highly abundant gene in each sample
t.gene_name(ord(end,:))

%% demonstrate that mito genes are problematic
ind_mito=find(startsWith(t.gene_name,'MT-'));
ind_mito=union(ind_mito,find(strcmpi(t.gene_name,'MTRNR2L12')));
ind_mito=union(ind_mito,find(strcmpi(t.gene_name,'MTRNR2L1')));
ind_mito=union(ind_mito,find(strcmpi(t.gene_name,'MTRNR2L8')));
ind_mito=union(ind_mito,find(strcmpi(t.gene_name,'AURKAIP1')));

sprintf("Mitochondrial genes count: %d",length(ind_mito))


s_mito=sum(t.mat_norm(ind_mito,:));
if only_mito_flag == 1
    figure;
    bar(s_mito);
    set(gca,'xtick',1:length(t.sample_name));
    set(gca,'XTickLabel',t.sample_name);
    set(gca,'XTickLabelRotation',45);
    title('sum mito genes');
end

%% filter ribosomal genes
if only_mito_flag == 0
    ind_ribo = union(find(startsWith(t.gene_name, "RPL")), find(startsWith(t.gene_name, "RPS")));
    if length(ind_ribo) > 0
        sprintf("Ribosomal genes count: %d",length(ind_ribo))
        s_ribo=sum(t.mat_norm(ind_ribo,:));
        figure;
        bar(1:length(s_ribo),[s_ribo;s_mito]);
        set(gca,'xtick',1:length(t.sample_name));
        set(gca,'XTickLabel',t.sample_name_print);
        set(gca,'XTickLabelRotation',45);
        title('Sum of ribosomal & mitochondrial genes');
        legend({"Ribosomal","Mitochondrial"})
        fprintf("Filtering out mitochondrial and ribosomal genes...\n")
    end
else
    fprintf("Filtering out mitochondrial genes...\n")
end

%% filter out mitochondrial genes
indin=setdiff(1:length(t.gene_name),ind_mito);

if only_mito_flag == 0
    indin = intersect(indin, setdiff(1:length(t.gene_name),ind_ribo));
end

t.gene_name_orig=t.gene_name;
t.gene_name=t.gene_name(indin);
t.seq_data=t.seq_data(indin,:);
t.mat_norm=t.seq_data./sum(t.seq_data);
[y,ord]=sort(t.mat_norm);
% figure;
% bar(y(end,:));
% set(gca,'xtick',1:length(t.sample_name));
% set(gca,'XTickLabel',t.sample_name);
% set(gca,'XTickLabelRotation',45);
% title('most highly expressed gene');


%% Extract only protein coding genes and lincRNAs
%T=readtable('X:\Common\useful_datasets\Mouse_GRcm38_91_ensemblBioMart_ref.csv');

T=readtable('X:\Common\useful_datasets\Human_GRch38_91_ensemblBioMart_ref.csv');
typ=table2cell(T(:,3));
gn=table2cell(T(:,2));
gns=gn(find(strcmpi(typ,'protein_coding')));

[~,indin] = intersect(t.gene_name, gns);
% indin=[];
% for i=1:length(t.gene_name)
%     indd=find(strcmpi(gns,t.gene_name{i}));
%     if ~isempty(indd)
%         indin=[indin i];
%     end
% end
t.gene_name=t.gene_name(indin);
t.seq_data=t.seq_data(indin,:);
t.mat_norm=t.seq_data./sum(t.seq_data);

%% examine the sum of remaining reads
SUM_UMI_THRESH=10000;
figure;bar(log10(sum(t.seq_data)));
set(gca,'xtick',1:length(t.sample_name));
set(gca,'XTickLabel',t.sample_name_print);
set(gca,'XTickLabelRotation',45);
ylabel('Log_1_0(Sum of UMIs)');
line(xlim,[log10(SUM_UMI_THRESH) log10(SUM_UMI_THRESH)],'color','k','linestyle','--','linewidth',2)
title("Sum of UMIs after filtration")
%set(gcf,'position',[1000         785        1432         553]);

%% Allon Klein Normalization
KLEIN_THRESH=0.1;
ind_norm_genes=find(max(t.mat_norm,[],2)<KLEIN_THRESH);
t.mat_norm=t.mat_norm./sum(t.mat_norm(ind_norm_genes,:));
