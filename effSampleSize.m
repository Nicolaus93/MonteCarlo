function ess = effSampleSize(w)
    [N, ~] = size(w);
    diff = bsxfun(@rdivide,N*w,sum(w))-1;
    CV = sqrt(sum(diff.^2)/N);
    ess = N./(1+CV.^2);