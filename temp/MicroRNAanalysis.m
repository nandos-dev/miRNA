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


%% The print-out of gene request above are limited to keep your output small. 
% you will need to define you own p value and or fold change thresholds
% to be used for selecting significantly different miRNAs for the steps
% below. Use all differentially expressed genes for the steps below, not
% just the ones you printed above. 

Foldchange = array2table(abs(FC));
Pvalues = array2table(pvals);
Alldata = [T_namesp, Foldchange, Pvalues]
Alldata.Properties.VariableNames{1} = 'miRNAs'
Alldata.Properties.VariableNames{2} = 'foldchange'
Alldata.Properties.VariableNames{3} = 'pvalues'

diff_miRNA = Alldata.foldchange >= 1.5 & Alldata.pvalues <=0.01
signf_miRNA = Alldata(diff_miRNA,:)

%% Find which mRNAs are the predicted targets of the significant miRNAs from
% the CRPS study using TargetScan. 

unique_targets = [];
targets = bmes_targetscandb_mir2target('hsa-miR-30d',0.9); % score threshold
fprintf('Number of gene from gene list found in these pathways \n')
length(targets)

%unique_targets.append(target,unique_targets);

%for item = 1:size(Alldata.('miRNAs'))
    %target = bmes_targetscandb_mir2target(item,0.9)
    %[~,idx] = unique(target.target,'rows','stable')
%end
%U = unique(size(Alldata(:,1)))
%for k = 1:Alldata.('miRNAs')
    %X = U(k) == size(Alldata(:,1))
   
%end

%% Perform enrichment of the targets using the DAVID webservice. List the 
% top 3 most significantly enriched pathways and top 3 most significantly 
% enriched Gene Ontology terms, along with their p-values and the number 
% of genes from your gene list found in these pathways and terms. 


%%%%%%%%%%%%%%% NOTES %%%%%%%%%%%%%%%%%%%%%%
%need to run each miRNA in forward loop that contained in the variable `signf_miRNA`
%to obtain the list of genes that are differentially expressed. 
%next take a unique list of these genes, and query them in david. 
%example of output genes ---- as follows -----

   % {'NM_001198820'}
   % {'NM_001168237'}
   % {'NM_005159'   }
    
% This will become the master list and perform the rich set enrichment
% genes, 
% Next and finally, search this query in david to produce the set of genes
% that are enriched. 
% in david you need to select GENBANK_ACCESSION and Gene List to accurately
% select the list from homo sapiens. 

figure(1)
imshow(imread('enriched_pathways.png'))
figure(2)
imshow(imread('Gene_Ontology_terms.png'))

%% Gene set enrichment with tens of genes will not get any enriched annotations 
%and the annotations from many thousands of genes will not be meaningful.
%For targetscanb, you can additionally adjust the confidence score
%threshold to control how many genes you get. You can also adjust pvalue or
%fold change thresholds above to control the number of most-different
%micrornas and this to limit the resulting gene set. 
























