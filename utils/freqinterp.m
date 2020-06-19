function [F,H]=freq_interp(F,H,opts)

if strcmp(optget(opts,'freqinterp',''),'lin')
    num_samples = optget(opts,'num_samples',length(F));
    F1=linspace(min(F),max(F),num_samples);
    H1=permute(H,[3,1,2]);
    m = size(H1,2);
    H=[];
    for c=1:m
        for d=1:m
            H(:,c,d)=interp1(F,H1(:,c,d),F1,'spline');
        end
    end
    F=F1;
    H=permute(H,[2,3,1]);
end
