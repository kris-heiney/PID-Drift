% Tuning curves are calulated by position bin and by averaging the activity
% of the mouse in the position bins

function tuningCurve = TuningCurve(fluor, behav, numBehavstates)
if size(fluor,1) < size(fluor,2)
    fluor = fluor';
end
if size(behav,1) < size(behav,2)
    behav = behav';
end

% Get tuning curve vectorized version
inds_mtx = bsxfun(@eq,behav, 1:numBehavstates);
N = sum(inds_mtx);
xatY_mtx = inds_mtx.*fluor;
tuningCurve = (sum(xatY_mtx,[],'omitnan')./N)';

end