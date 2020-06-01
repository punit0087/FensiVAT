function [smp,rp,m] = MMRS(x, cp, ns )
[n,p]=size(x);

[m,rp]=MM(x,cp); %%MM sampling 

[~,i]=min(rp,[],2);
smp=[];
for t=1:cp,
    t
    s = find(i==t); 
    nt = ceil(ns*length(s)/n); %
    ind = ceil(rand(nt,1)*length(s));  %%Random Sampling
    smp=[smp; s(ind)];
end;
smp=unique(smp);
end