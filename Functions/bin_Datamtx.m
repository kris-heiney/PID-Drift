function bin_Datamtx(file, numBhvStates, savePath, isHeading)

load(file,'Data_mtx');

% Parameters
numSessions=max(Data_mtx(:,end));
numFluorStates = 3;
Data_mtx_binned = nan(size(Data_mtx));
Data_mtx_binned(:,end-1:end) = Data_mtx(:,end-1:end);

% Bin the position and the fluorescence data
fluor = Data_mtx(:,1:end-4);
behavior = Data_mtx(:,end-3);
fluor_states_all = binfluor(fluor, numFluorStates, 'uniform', 1); 
bhv_states_all = binbehav(behavior, numBhvStates, 'activity', isHeading, fluor);

for trialPick = 2:3
    for s = 1:numSessions
        ind_s = (Data_mtx(:,end-1)==trialPick) & Data_mtx(:,end)==s;
        bhv_states = bhv_states_all(ind_s);
        fluor_states = fluor_states_all(ind_s,:);
        Data_mtx_binned(ind_s,1:end-4) = fluor_states;
        Data_mtx_binned(ind_s,end-3) = bhv_states;
    end
end
saveFilename = strcat(savePath, file(1:end-4), '_binned.csv');
writematrix(Data_mtx_binned, saveFilename)

end

