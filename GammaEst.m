function  bm=GammaEst(X0,r,w,alpha)

n=500; r0=-log(1-r);
b=(linspace(r0,5,n))';
sumB=zeros(length(b),1);
for j=1:w-1
    t=w-j; 
    SA=gammainc(b.*t,alpha).*X0(j);
    SB=X0(j)-(1-r)^t*X0(j);
    sumB=sumB+SA-SB;
end
m=find(sumB>0,1);
%[~,m]=min(abs(sumB));
if m<n
    bm=b(m);
else
    disp(m)
    error('no local optimal point')
end
