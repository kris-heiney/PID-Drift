function Xcut = CutOutliers(X) 
%% big outliers might be artifacts and might bias the results
% X is fluorescence data, where sessions are merged: cell x time
% Xcut is a filtered X where the outliers are cut at a certain level
if size(X,1)<size(X,2)
    X=X';
end

Xcut=X;
Xmean=mean(X,'omitnan');
X0mean=(X'-Xmean')';
Xmedian=mean(X0mean,'omitnan');
active_cellind=find(~isnan(Xmean));
for i=1:length(active_cellind)
    ind_i=active_cellind(i);
    min_vec=mink(X(~isnan(X(:,ind_i)),ind_i),5);
    max_vec=maxk(X(~isnan(X(:,ind_i)),ind_i),5);
    bounds=1.5*[Xmedian(ind_i)+prctile(X0mean(:,ind_i),5),Xmedian(ind_i)+prctile(X0mean(:,ind_i),99)];% asymmetric because positive jumps in the data are informative about spikes but it is hard to interpret the negative jumps so we are stricter about cutting those
    bounds(1)=max(min_vec(end),bounds(1));
    bounds(2)=min(max_vec(end),bounds(2));
    Xcut(X(:,ind_i)>bounds(2),ind_i)=bounds(2);
    Xcut(X(:,ind_i)<bounds(1),ind_i)=bounds(1);
end
end