function outdata = Switch_MergeExperiments(indata,pairedsessions,experimentID,mergesides,sideindex)
    outdata = indata;
    if mergesides
        dimsize = ndims(indata);
        outdata = permute(outdata,[1 sideindex setdiff([2:dimsize],sideindex)]);
        outdata(experimentID==2,1,:,:,:,:,:,:,:,:,:,:) = mean(outdata(experimentID==2,1:2,:,:,:,:,:,:,:,:,:,:),2);
        outdata = squeeze(outdata(:,1,:,:,:,:,:,:,:,:,:,:));
    end

    %average repeated measures
    for ipair = 1:size(pairedsessions,2)
        outdata(pairedsessions(1,ipair),:,:,:,:,:,:,:,:,:,:) = mean(outdata(pairedsessions(1:2,ipair),:,:,:,:,:,:,:,:,:,:),1);
    end
    outdata(pairedsessions(2,:),:,:,:,:,:,:,:,:,:,:) = [];
end