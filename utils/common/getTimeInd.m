function time_index = getTimeInd(time_vector,limits)
%function time_index = getTimeInd(time_vector,limits)
%time_vector: 1xn vector of time points
%limits: 1x2 vector of time range limits
%time_index: 1xn logical vector with ones between limits(1) and limits(2)
    time_index = time_vector>=limits(1) & time_vector<=limits(2);
end
