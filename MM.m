function [m,Rp] = MM(X,cp)

[n,p]=size(X);

m=ones(cp,1);
m(1)=ceil(rand(1)*size(X,1)); %%radonmly choose first point
d=distance2(X(m(1),:),X)';
Rp(:,1)=d;
for t=2:cp,
    t;
    d=min(d,Rp(:,t-1));
    [~,m(t)]=max(d);
    Rp(:,t)=distance2(X(m(t),:),X)';
end;