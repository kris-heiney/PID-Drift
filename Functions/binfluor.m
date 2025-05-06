% Function to 'state' or discretize the fluorescence data

function [fluor_binned, bins] = binfluor(fluor, numStates, binmethod, iscut)

% Check if fluor is number of samples x number of cells. If not then
% transpose.
if size(fluor,2) > size(fluor,1)
    fluor = fluor';
end
numCells = size(fluor,2);
numSamples = size(fluor,1);

if ~exist('iscut','var')
    iscut=0;
end
if iscut==1
    fluor_cut = CutOutliers(fluor);
else
    fluor_cut = fluor;
end

% Define bins for states
fluor_binned = nan(size(fluor));
bins=nan(numCells,numStates+1);
bins(:,1)=min(fluor_cut);
bins(:,end)=max(fluor_cut)';
for cel=1:numCells
    bins(cel,2)=prctile(abs(fluor_cut(fluor_cut(:,cel)<0,cel)),99);
end
fluor_mean=mean(fluor,'omitnan');
active_cellind=find(~isnan(fluor_mean));
switch binmethod
    case 'uniform'
        if numStates>2
            for cel=1:length(active_cellind)
                celind=active_cellind(cel);
                bins(celind,2:end) = linspace(bins(celind,2), bins(celind,end),numStates);
                fluor_binned(:,celind) = discretize(fluor_cut(:,celind), bins(celind,:));
            end
        end
    case 'percentile'
        for cel = 1:length(active_cellind)
            celind=active_cellind(cel);        
            fluor_cel_pos = fluor_cut(fluor_cut(:,celind)>bins(celind,2),celind);
            bins(celind,1)=min(fluor_cut(:,celind),[],'omitnan');
            if numStates==3
                bins(celind,2)=prctile(fluor_cel_pos,80);
                bins(celind,3)=prctile(fluor_cel_pos,90);
                bins(celind,4)=max(fluor_cel_pos,[],'omitnan');
            elseif numStates==4
                bins(celind,2)=prctile(fluor_cel_pos,72);
                bins(celind,3)=prctile(fluor_cel_pos,84);
                bins(celind,4)=prctile(fluor_cel_pos,96);
                bins(celind,5)=max(fluor_cel_pos,[],'omitnan');
            elseif numStates==5
                bins(celind,2)=prctile(fluor_cel_pos,65);
                bins(celind,3)=prctile(fluor_cel_pos,76);
                bins(celind,4)=prctile(fluor_cel_pos,87);
                bins(celind,5)=prctile(fluor_cel_pos,98);
                bins(celind,6)=max(fluor_cel_pos,[],'omitnan');
            end
            fluor_binned(:,celind) = discretize(fluor(:,celind), bins(celind,:));
        end
end
end