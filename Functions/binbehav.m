% Function to 'state' the data (state bins)
% cutoffmethod: 'percent' or 'abs'
%     'percent' cuts the upper and lower percentages of the data given in
%     x_cutoffs
%     'abs' uses x_cutoffs as absolute thresholds
% binmethod: 'uniform' or 'activity'
%     'uniform' discretizes x into numStates bins of uniform width
%     'activity' discretizes x so that the amount of activity (in fluor) is
%     uniform across the bins (use only for stimulus data)
% numStates can be 0 to allow for states to be 0:1:max(spikecount)


function [x_states, bins] = binbehav(x, numStates, binmethod, isHeading, fluor)

% Check if x a vector (behaviour) or matrix (neuron response)
if size(x,1)<size(x,2)
    x=x';
end
if exist("fluor","var")
    if size(fluor,1)<size(fluor,2)
        fluor=fluor';
    end
end

% If heading data, cut everything beyond 90 degrees
if isHeading == 1
    x(x>90) = 90;
    x(x<-90) = -90;
end

% Allow for number of bins to not be predetermined; use if data is spike
% counts to allow for counts to be integers up to max spike count
if numStates == 0
    numStates = ceil(max(x));
end

% Define bins for states
Xmax = max(x,[],'omitnan');
Xmin = min(x,[],'omitnan');

% Find bin edges
switch binmethod
    % Uniform bins - This amlifies movements in the mid-range, where the
    % cell is silent
    case 'uniform'
        x_states = nan(size(x));
        binSize = (Xmax - Xmin)/numStates;
        if ~isnan(binSize)
            bins = Xmin:binSize:Xmax; 
            x_states = discretize(x, bins);
        else
            disp('Problem occured in binbehav during uniform binning') %this should not happen as 
        end
    % Activity-dependent bins
    case 'activity'
        [x_sort, sortinds] = sort(x, 'ascend');
        % Sum activity over cells and sort by stimulus
        popActivity = sum(abs(fluor), 2, 'omitnan');
        popActivitySort = popActivity(sortinds);
        % Cumulative sum of stimulus-sorted activity
        cPopActivity = cumsum(popActivitySort);
        % Find where cumulative activity exceeds percentage of total
        % cumulative activity
        bins = zeros(numStates-1,1);
        for b = 1:(numStates-1)
            ind = find(cPopActivity > cPopActivity(end)*b/numStates, 1);
            bins(b) = x_sort(ind);
        end
        bins = [Xmin; bins; Xmax];
        % Discretize variables into state bins
        x_states = discretize(x, bins);        
end

end