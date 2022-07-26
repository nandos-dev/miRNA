%% MicroRNA analysis 
% Template by [Fernando Ramirez & Thinh Nguygen]
% Code by [Fernando Ramirez & Thinh Nguygen]

%% Use the unfiltered Excel data Remove anay miRNA that is detected in 3 or fewer samples
CT_values = readmatrix('CRPS_unfiltered.xlsx');
miRNAs = readvars('CRPS_unfiltered.xlsx');
CT_values(:,1) = []; %extracting the all rows of first column 
vars_miRNAs = readcell('CRPS_unfiltered.xlsx');
%separating out control and patients 

control = find(contains(vars_miRNAs(1,:),'control')); 
patient = find(contains(vars_miRNAs(1,:), 'patient')); 

elimrows = true(size(CT_values,1),1); %removes any row taht does now contain value 
for row = 1:size(CT_values,1) %iterate over the size of the matrix 
    if size(unique(CT_values(row,:)),2) <=4 %argument of 2 
        elimrows(row) = false; %elimination of the rows that return 0 as false
    end
end
CT_values = CT_values(elimrows,:);
miRNAs = miRNAs(elimrows,:); 


%% Replace any undetected value (Inf) with the average of the expression of
% that miRNA in the rest of the samples

temparr = []; tempave = []; TF = []; %prealocate double to store output 

for row = 1:size(CT_values,1) 
    temparr = CT_values(row,:);
    TF = isinf(temparr);
    temparr = temparr(~TF);
    tempave = mean(temparr);
    for col = 1:size(CT_values,2)
        if isinf(CT_values(row,col))
            CT_values(row,col) = tempave; 
        end
    end
end


%% Print the first 5 row (genes) and first 6 columns (samples) of your data
% of your data 

T_names = table(miRNAs); 
T_data = array2table(CT_values);
T = [T_names T_data];
fprintf('First 5 genes and 6 samples of the filtered data\n')
disp(T(1:5,1:6))


% You should end up with datda that is the same as the filtered Excel data. 
% If you get different filtered data, and cannot find 
% out why you are getting different results, abandon your data, use the 
% filtered Excel data. 
%%%%% Filtered Data Matched! %%%%%


%% For each sample, find the CT0 values by averaging the CT values of RNU44,
% RNU48, and MammU6. Subtract the CT values of that sample
% with the CT0 of that sample. 
expr = 'RNU44|RNU48|MammU6'; indexrow = zeros(size(miRNAs));
%prealocate the indexrow (miRNAs) 

%iterative over size of miRNAS using regex, extract these values 
for row = 1:size(miRNAs)
    temp = regexp(miRNAs{row},expr, 'once');
    if isempty(temp)
        indexrow(row) = 0; %does not match the expression (expr) 
    else
        indexrow(row) = row; %does match, output row index 
    end
end
temp = find(indexrow == 0); indexrow(temp,:) = [];
%indicies of the expression values to calculate average with. 

CT0 = mean(CT_values(indexrow,:),1); %get CT0 values, across rows `1`
deltaCT_norm = CT_values;

for row = 1:size(deltaCT_norm,1)
    deltaCT_norm(row,:) = deltaCT_norm(row,:) - CT0; 
end
%% The result of the previous step will give you normalized deltaCT values. 
%%%%% Filtered Data Matched! %%%%%
%% For each miRNA, find the average of all healthy deltaCT values and 
% separately, find the average of the patient deltaCT values 
deltaCT_ctrl = zeros(size(deltaCT_norm,1),1); % preallocation for healthy values
deltaCT_pat = zeros(size(deltaCT_norm,1),1); % preallocation for patients values

for row = 1:size(deltaCT_norm,1)
    deltaCT_ctrl(row) = mean(deltaCT_norm(row,control-1)); 
    deltaCT_pat(row) = mean(deltaCT_norm(row,patient-1));
end
%subtracting one from control and patient 


%% For each miRNA, subtract the average healthy deltaCT values from the 
% average patient deltaCT values. This gives you the deltadelta CT. 

deltadeltaCT = deltaCT_pat - deltaCT_ctrl;
%size of the matrices match 

%% Calculate 2^deltadeltaCT, which are the fold change values. 

FC = 2.^(-deltadeltaCT);

%% Replace andy fold change value x that is less than 1, with its negative inverse (-1/x)
for row = 1:size(FC,1) %size across rows 
    if FC(row)<1
        FC(row) = (-1/FC(row));
    end
end
%few number of genes -- 49 genes that are  511 that are upregulated. 
%% Print the names and fold changes of the top-10 most changing miRNAs 
%(either up or down regulaltion, i.e., order by decending
%absolute fold change)
[~,FCindex] = sort(abs(FC),'descend'); %absolute value of the fold changes
% to print 
miRNAnames = table(miRNAs(FCindex)); 
miRNAnames.Properties.VariableNames{1} = 'miRNAs'; %rename column
T_FC = array2table(FC(FCindex)); %printing display to formatted table 
T_FC.Properties.VariableNames{1} = 'Fold Changes Values'; % rename column
TFC_fold = [miRNAnames, T_FC];
TFC_fold.Properties.Description = 'Top 10 Most Changing miRNAs';
fprintf('Fold Changes of Top 10 Most Changing miRNAs\n')
disp(TFC_fold(1:10,:))


%% Using the deltaCT values, find the significantly different miRNAs between 
%controls and patients. 


[pvals] = mattest(deltaCT_norm(:,control-1),deltaCT_norm(:,patient-1)); %shuffle based test 


%% Print the names and p-values of the top-10 most significantly different 
% miRNAs (ordered by pvalue)

[~,pindx] = sort(pvals,'descend');
T_namesp = table(miRNAs(pindx)); 
T_namesp.Properties.VariableNames{1} = 'miRNAs'; %rename column
T_pvals = array2table(pvals(pindx));
T_pvals.Properties.VariableNames{1} = 'P-Values'; % rename column
TFC_pvalue = [T_namesp T_pvals];
TFC_pvalue.Properties.Description = 'Top 10 Most Significantly Different miRNAs';
fprintf('p-values Top 10 Most Significantly Different miRNAs\n')
disp(TFC_pvalue(1:10,:))

%Detemine the similarities across tables 

%newnames = matlab.lang.makeUniqueStrings([TFC_fold.Properties.VariableNames, TFC_pvalue.Properties.VariableNames]);
%TFC_fold.Properties.VariableNames = newnames(1:numel(TFC_fold.Properties.VariableNames));
%TFC_pvalue.Properties.VariableNames = newnames(numel(TFC_fold.Properties.VariableNames)+1:end);
%T3 = [TFC_fold, TFC_pvalue];

%% The print-out of gene request above are limited to keep your output small. 
% you will need to define you own p value and or fold change thresholds
% to be used for selecting significantly different miRNAs for the steps
% below. Use all differentially expressed genes for the steps below, not
% just the ones you printed above. 

Foldchange = array2table(abs(FC));
Pvalues = array2table(pvals);
Alldata = [T_namesp, Foldchange, Pvalues]
Alldata.Properties.VariableNames = {'miRNA','foldchange','pvalues'};

%miRNAs_of_interest = Alldata(:,'pvalues')<=0.01

%I = signif_dpvals(:,'p-values')<=0.01 & abs(signif_dpvals(:,'negfc'))>=1.5;
%miRNAs_of_interest = (Alldata.foldchange >= 1.00 & Alldata.pvalues <=0.05)
%Alldata(miRNAs_of_interest,:)


%Adding fold change information to the dpvals object. 
%signif_pvals = pvals(pvals(:,1) <= 0.01,:); 
%signif_pvals = Foldchange(Foldchange(:,1) >= 1,:); 
%I = pvals(pvals(:,1) <= 0.01,:) & abs(FC(:,FC))>=1.5; 

%signifance of 137 out of 560 selected miRNAs 
%signif_pvals = [pvals bioma.data.DataMatrix(FC,pvalues, 'ColNames',{'fc'},{'pvals'})]
%I = abs(signif_pvals(:,'fc'))>=1.5;
%signif_pvals = [pvals bioma.data.DataMatrix(FC, 'ColNames',{'pvals'})]

%dsigfc = signif_pvals(I,:);
%dsigfc = dsigfc.sortrows('fc');
%fprintf('Found %d genes with pvalue<=0.01 and FC>=1.5. Showing top 10:\n',size(dsigfc,1));
%disp(dsigfc)


%% Find which mRNAs are the predicted targets of the significant miRNAs from
% the CRPS study using TargetScan. 


[targets] = bmes_targetscandb_mir2target(); 
[targets] = bmes_targetscandb_mir2target('hsa-miR-30a',0.9);
[targets] = bmes_targetscandb_mir2target('hsa-miR-03a',0.9);

% generefseqid is equal to the the targets 
% take the list of the genes and perform the gene-set enrichment 
% take the miRNA and take the genes that are predicted to target. 

%score threshold of .90, 
%function creates an SQL database --> construct the query to get to the query, ask 
%databased, to produce the 
%go to worksheet called the mir2target, column generefseq ID (column), Score (column)

%delete these comments
%background -- target scan website, purpose of automation 
%download code and analysis tooling 
%novel -- using 
%sequence similarities, miRNA (highly similar) -- e.g. mir-20-5p and -20-3p 
%mir20, mir30 (sequences), may be apart of same family (unfrequent) 
%minimizing the size of the database, characters in the list (in the mir20 familiy) 
%sacan created miRNA to Target Gene can be obtained from single target 
%sql lite database only contains one table 
%thus, finding gene targets of any of the hsa-miR-30a, and 30ap, or -50a
%improvements to check to see if the first hsa-miR-30a-3p already exist no need to check the 
%other gene sequences. 
%d.query('SELECT generefseqid FROM mir2target WHERE mirna IN ("hsa-miR-30a")')
%install the database extensions for obtain the col1 as the table headings 
%target scan is a prediction method. 
%


%% Perform enrichment of the targets using the DAVID webservice. List the 
% top 3 most significantly enriched pathways and top 3 most significantly 
% enriched Gene Ontology terms, along with their p-values and the number 
% of genes from your gene list found in these pathways and terms. 







%% Gene set enrichment with tens of genes will not get any enriched annotations 
%and the annotations from many thousands of genes will not be meaningful.
%For targetscanb, you can additionally adjust the confidence score
%threshold to control how many genes you get. You can also adjust pvalue or
%fold change thresholds above to control the number of most-different
%micrornas and this to limit the resulting gene set. 
























