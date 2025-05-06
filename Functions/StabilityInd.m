function [meanTuningStability, meanDayTuningStability, tuningStability, deltaSession, tuningCurves] = StabilityInd(file,trialPick,numBhvStates,numPerm,doBinning,isHeading)
%% This function calculates drift metrics
rng('shuffle')

% Load file
load(file,'Data_mtx');

% Parameters
numSessions = max(Data_mtx(:,end));
numCells = size(Data_mtx,2)-4;

% Daily drift
[daily_drift, tuningCurves] = daily_TCdiff(file, trialPick, numPerm, numBhvStates, doBinning, isHeading);

%% Find tuning curve change between pairs of days
deltaSession = nan(numCells, nchoosek(numSessions, 2));
tuningStability = nan(numCells, nchoosek(numSessions, 2), 5);
meanTuningStability = nan(numCells, 5);
meanDayTuningStability = nan(numCells, 5);
% tuningStability is a n x ds x 4 cell array where
% n is the number of neurons
% ds is the number of session pairs
% All cases have compensation for within-day variability
% (:,:,1) is absolute difference
% (:,:,2) is absolute difference after z-scoring
% (:,:,3) is cosine similiarity
% (:,:,4) is cosine similarity after z-scoring

for cel = 1:numCells
    k = 1;
    for s0 = 1:numSessions - 1
        for s1 = s0+1:numSessions
            % Tuning curves
            tuning0 = tuningCurves(:,s0,cel);
            tuning1 = tuningCurves(:,s1,cel);
            keepbins = ~isnan(tuning0) & ~isnan(tuning1);
            tuning0 = tuning0(keepbins);
            tuning1 = tuning1(keepbins);
            % 1: Absolute difference
            tuningStability(cel,k,1) = -sum(abs(tuning1 - tuning0))/numBhvStates ...
                + (0.5*(mean(daily_drift(s0,cel,:,1))+mean(daily_drift(s1,cel,:,1))));
            % 3: Cosine similarity
            tuning0_norm = tuning0/norm(tuning0);
            tuning1_norm = tuning1/norm(tuning1);
            tuningStability(cel,k,3) = dot(tuning0_norm, tuning1_norm)/(0.5*(mean(daily_drift(s0,cel,:,3))+mean(daily_drift(s1,cel,:,3))));
            % Z-score tuning curves for each (cell, session) pair
            tun0_z = (tuning0 - mean(tuning0))/std(tuning0);
            tun1_z = (tuning1 - mean(tuning1))/std(tuning1);
            % 2: Absolute difference of z-scored
            tuningStability(cel,k,2) = -sum(abs(tun1_z - tun0_z))/numBhvStates ...
                + (0.5*(mean(daily_drift(s0,cel,:,2))+mean(daily_drift(s1,cel,:,2))));
            % 4: Cosine similarity
            tun0_norm = tun0_z/norm(tun0_z);
            tun1_norm = tun1_z/norm(tun1_z);
            tuningStability(cel,k,4) = dot(tun0_norm, tun1_norm)/(0.5*(mean(daily_drift(s0,cel,:,4))+mean(daily_drift(s1,cel,:,4))));
            % 5: Pearson correlation
            CCbs = corr(tun0_z, tun1_z);
            CCws = 0.5*(mean(daily_drift(s0,cel,:,5), 'omitnan') + mean(daily_drift(s1,cel,:,5), 'omitnan'));
            %CCws = mean(daily_drift(s0,cel,:,5), 'omitnan');
            if CCws > CCbs && CCws > 0 && CCbs > 0
                tuningStability(cel,k,5) = 1 - (CCws - CCbs)/(CCws + CCbs);
            elseif CCbs < 0
                tuningStability(cel,k,5) = 0;
            end
            
            % Change in days
            deltaSession(cel,k) = s1 - s0;
            k = k+1;
        end
    end
    
    % Average over session pairs for all cells
    ts = squeeze(tuningStability(cel,:,:));
    meanTuningStability(cel,:) = mean(ts, 1, 'omitnan');
    meanDayTuningStability(cel,:) = mean(ts.*deltaSession(cel,:)', 1, 'omitnan');%sum up the d-day differences and weigh by how many samples we have for each day difference
end

end