
function [daily_drift, tuningcurves] = daily_TCdiff(file, trialPick, numPerm, numBhvStates, doBinning, isHeading)
%% This function is similar to TCdiff0_MJ
load(file,'Data_mtx');
rng('shuffle')

% Trial type
inds = Data_mtx(:,end-1)==trialPick;
if doBinning == 0
    bhv_states = Data_mtx(inds,end-3);
else
    bhv_states = binbehav(Data_mtx(inds,end-3), numBhvStates, 'uniform', isHeading);
end

% Parameters
numSessions = max(Data_mtx(:,end));
numCells = size(Data_mtx,2)-4;
if ~exist('numPerm','var') 
    numPerm = 30;
end

% Initialize variables
daily_drift=nan(numSessions,numCells,numPerm,5);
tuningcurves=nan(numBhvStates,numSessions,numCells);

% Get tuning curves per session per cell
tic
for s = 1:numSessions
    ind_s = Data_mtx(inds,end)==s;
    bhv_states_s = bhv_states(ind_s);
    ind_s = Data_mtx(:,end-1)==trialPick & Data_mtx(:,end)==s;
    fluor_s = Data_mtx(ind_s,1:end-4);
    trials=Data_mtx(ind_s,end-2);
    for cel=1:numCells
        tuningcurves(:,s,cel)=TuningCurve(fluor_s(:,cel), bhv_states_s, numBhvStates);
        for i=1:numPerm
            M=randsample(length(trials),floor(length(trials)/2),false);
            tc1=TuningCurve(fluor_s(M,cel), bhv_states_s(M), numBhvStates);
            tc2=TuningCurve(fluor_s(setdiff(1:length(trials),M),cel), bhv_states_s(setdiff(1:length(trials),M)), numBhvStates);
            tc1_z=zscore(tc1);
            tc2_z=zscore(tc2);
            keepbins = ~isnan(tc1) & ~isnan(tc2);
            tc1 = tc1(keepbins);
            tc2 = tc2(keepbins);
            tc1_z = tc1_z(keepbins);
            tc2_z = tc2_z(keepbins);

            daily_drift(s,cel,i,1) = sum(abs(tc1-tc2))/numBhvStates;
            daily_drift(s,cel,i,2) = sum(abs(tc1_z-tc2_z))/numBhvStates;
            daily_drift(s,cel,i,3) = dot(tc1/norm(tc1), tc2/norm(tc2));
            daily_drift(s,cel,i,4) = dot(tc1_z/norm(tc1_z), tc2_z/norm(tc2_z));
            daily_drift(s,cel,i,5) = corr(tc1_z,tc2_z);
        end
    end
    disp(s)
    toc
end

end

function z = zscore(x)

z = (x - mean(x,'omitmissing'))/std(x,'omitmissing');
end