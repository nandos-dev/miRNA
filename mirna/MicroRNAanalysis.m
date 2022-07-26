%% MicroRNA analysis 
% Template by [Fernando Ramirez & Thinh Nguyen]
% Code by [Fernando Ramirez & Thinh Nguyen]

%% Use the unfiltered Excel data
% matrix is a much more optimizing method to import data
data_miRNAs = readmatrix('CRPS_unfiltered.xlsx');
miRNAs = readvars('CRPS_unfiltered.xlsx');
data_miRNAs(:,1) = [];
vars_miRNAs = readcell('CRPS_unfiltered.xlsx');

control = find(contains(vars_miRNAs(1,:),'control'));
patients = find(contains(vars_miRNAs(1,:),'patient'));

       
%% Remove anay miRNA that is detected in 3 or fewer samples
elimrows = true(size(data_miRNAs,1),1);
for row = 1:size(data_miRNAs,1)
    if size(unique(data_miRNAs(row,:)),2) <= 4
        elimrows(row) = false;
    end
end
data_miRNAs = data_miRNAs(elimrows,:);
miRNAs = miRNAs(elimrows,:);

%% Replace any undetected value (Inf) with the average of the expression of
% that miRNA in the rest of the samples
temparr = []; tempave = []; TF = [];
for row = 1:size(data_miRNAs,1)
    temparr = data_miRNAs(row,:);
    TF = isinf(temparr);
    temparr = temparr(~TF);
    tempave = mean(temparr);
    for col = 1:size(data_miRNAs,2)
        if isinf(data_miRNAs(row,col))
            data_miRNAs(row,col) = tempave;
        end
    end
end

%% Print the first 5 row (genes) and first 6 columns (samples) of your data
% of your data
T_names = table(miRNAs);
T_data = array2table(data_miRNAs);
T = [T_names T_data];
fprintf('First 5 genes and 6 samples of the filtered data\n')
disp(T(1:5,1:6))

%% You should end up with datda that is the same as the filtered Excel data. 
% If you get different filtered data, and cannot find 
% out why you are getting different results, abandon your data, use the 
% filtered Excel data.
%%%%% Filtered Data Matched! %%%%%

%% For each sample, find the CT0 values by averaging the CT values of RNU44,
% RNU48, and MammU6. Subtract the CT values of that sample
% with the CT0 of that sample.
% Finding RNU44 RNU48 & MammU6 data rows
expr = 'RNU44|RNU48|MammU6'; indexrow = zeros(size(miRNAs));
for row = 1:size(miRNAs)
    temp = regexp(miRNAs{row},expr, 'once');
    if isempty(temp)
        indexrow(row) = 0;
    else
        indexrow(row) = row;
    end
end
temp = find(indexrow == 0); indexrow(temp,:) = [];
CT0 = mean(data_miRNAs(indexrow,:),1); % get CT0 values
deltaCT_norm = data_miRNAs;
for row = 1:size(deltaCT_norm,1)
    deltaCT_norm(row,:) = deltaCT_norm(row,:) - CT0;
end

%% The result of the previous step will give you normalized deltaCT values. 
%%%%% Normalized deltaCT values obtained! %%%%%
%% For each miRNA, find the average of all healthy deltaCT values and 
% separately, find the average of the patient deltaCT values
deltaCT_ctrl = zeros(size(deltaCT_norm,1),1); % preallocation for healthy values
deltaCT_pat = zeros(size(deltaCT_norm,1),1); % preallocation for patients values
for row = 1:size(deltaCT_norm,1)
    deltaCT_ctrl(row) = mean(deltaCT_norm(row,control-1));
    deltaCT_pat(row) = mean(deltaCT_norm(row,patients-1));
end

%% For each miRNA, subtract the average healthy deltaCT values from the 
% average patient deltaCT values. This gives you the deltadelta CT. 
deltadeltaCT = deltaCT_pat - deltaCT_ctrl;

%% Calculate 2^deltadeltaCT, which are the fold change values. 
FC = 2.^(-deltadeltaCT);
%% Replace andy fold change value x that is less than 1, with its negative inverse (-1/x)
for row = 1:size(FC,1)
    if FC(row) < 1
        FC(row) = (-1/FC(row));
    end
end
%% Print the names and fold changes of the top-10 most changing miRNAs 
%(either up or down regulaltion, i.e., order by decending
%absolute fold change)
[~,FCindex] = sort(abs(FC),'descend');
T_namesFC = table(miRNAs(FCindex)); 
T_namesFC.Properties.VariableNames{1} = 'miRNAs'; %rename column
T_FC = array2table(FC(FCindex));
T_FC.Properties.VariableNames{1} = 'Fold Changes Values'; % rename column
TFC = [T_namesFC T_FC];
TFC.Properties.Description = 'Top 10 Most Changing miRNAs';
fprintf('Top 10 Most Changing miRNAs\n')
disp(TFC(1:10,:))

%% Using the deltaCT values, find the significantly different miRNAs between 
%controls and patients. 
pvals = mattest(deltaCT_norm(:,control-1),deltaCT_norm(:,patients-1));
%% Print the names and p-values of the top-10 most significantly different 
% miRNAs (ordered by pvalue)
[~,pindex] = sort(pvals,'descend');
T_namesp = table(miRNAs(pindex)); 
T_namesp.Properties.VariableNames{1} = 'miRNAs'; %rename column
T_pvals = array2table(pvals(pindex));
T_pvals.Properties.VariableNames{1} = 'P-Values'; % rename column
TFC = [T_namesp T_pvals];
TFC.Properties.Description = 'Top 10 Most Significantly Different miRNAs';
fprintf('Top 10 Most Significantly Different miRNAs\n')
disp(TFC(1:10,:))

%% The print-out of gene request above are limited to keep your output small. 
% you will need to define you own p value and or fold change thresholds
% to be used for selecting significantly different miRNAs for the steps
% below. Use all differentially expressed genes for the steps below, not
% just the ones you printed above. 

newdata = cat(2,miRNAs,num2cell(FC));
newdata = cat(2,newdata,num2cell(pvals));
TF_miRNA = false(size(newdata,1),1);
for row = 1:size(newdata,1)
    if (abs(newdata{row,2}) >= 1.5) && (newdata{row,3} <= 0.01)
        TF_miRNA(row) = true;
    end
end

signf_miRNA = newdata(TF_miRNA,:);

%% Find which mRNAs are the predicted targets of the significant miRNAs from
% the CRPS study using TargetScan. Sequelquery 
targettemp = [];newlist = {};
for row = 1:size(signf_miRNA,1)
    tempstruct = bmes_targetscandb_mir2target(newdata{row,1},0.9);
    if isempty(tempstruct) ~= true
        for row = 1:size(tempstruct)
            newlist = [newlist; tempstruct{row}];
        end
    end
end
newlist = unique(string(newlist));
fprintf('Number of genes from gene list found in pathways is %d',length(newlist))

%% Perform enrichment of the targets using the DAVID webservice. List the 
% top 3 most significantly enriched pathways and top 3 most significantly 
% enriched Gene Ontology terms, along with their p-values and the number 
% of genes from your gene list found in these pathways and terms. 

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





















