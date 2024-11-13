function [gamma,dI,b0,statvals,OLSvals]=RecoveryEstimate(M,w0,a,d)
% Recovery parameter estimation
% M - Data matrix,  w0 - window sizes, 
% a - Gamma dist parameter, d - discretize (1 = yes)

% Number of data points
L=length(M(:,1));

% M = Deaths	N	S	I	R
D=M(:,1); I=M(:,4); R=M(:,5); 

% Forward differencing
dR=R(2:end)-R(1:end-1);
%Backward differencing
%dR=[dR(2:end);0];

dI=zeros(L-1,1);


% Comulative cases (CC) and incidents (Inc)
CC=I+R+D; Inc=CC(2:end)-CC(1:end-1);

I=I(1:end-1);
p=L-1; q=length(w0);

% Discreitization of a
if d==1
    n=41; 
    alpha=linspace(1,a,n);  
    gamma = zeros(p,n+1,q);
    b0=zeros(p,n+1,q);
else
    alpha=a;
    gamma = zeros(p,3,q);
    b0=zeros(p,3,q);
end

% OLS results:
% Uncentered R squared
R2y1=zeros(p,q);
% Standard errors
SEy1=zeros(p,q);
% MSE of residuals
MSEy1=zeros(p,q);
% Mean of residuals
Meany1=zeros(p,q);

% Kolmogorov-Smirnov test results
% Null hypothesis test is that the data comes from a standard normal distribution
KSy1=zeros(p,q); 
KSy1p=zeros(p,q);

% Shapiro-Wilk parametric hypothesis test results
% Null hypothesis test is that the data comes from a normal distribution
SWy1=zeros(p,q); SWy1p=zeros(p,q);

for j0=1:q

w=w0(j0);

for i=1:p-w+1
    ep = w+i-1; % estimation point
    
    % dR
    X = I(i:w+i-1).^0.5;
    y = dR(i:w+i-1)./X; 
    Rmdl = fitlm(X,y,'Intercept',false);

    % Regression quality
    SEy1(ep,j0)   = Rmdl.Coefficients.SE;
    MSEy1(ep,j0)  = Rmdl.RMSE;
    ResR=table2array(Rmdl.Residuals(:,1));
    z1= ResR(2:end)-ResR(1:end-1); 
    z1 = z1(~isnan(z1));
    Meany1(ep,j0) = mean(z1);
    R2y1(ep,j0) =1-(ResR'*ResR)/sum((y).^2); % Uncentered

    % Kolmogorov-Smirnov and Shapiro-Wilk tests
    if std(z1) ==0
        z1=0;
    else
        z1=z1./std(z1);
    end
    
    if length(z1)>2
        [KSy1(ep,j0),KSy1p(ep,j0)]=kstest(z1);
    else
        KSy1(ep,j0)=0;
        KSy1p(ep,j0)=1;
    end

    z1=z1-mean(z1);

    if length(z1)>2
        [SWy1(ep,j0),SWy1p(ep,j0), ~] = swtest(z1, 0.05);
    else
        SWy1(ep,j0)=0;
        SWy1p(ep,j0)=1;
    end
    
    r = Rmdl.Coefficients.Estimate;
    gamma(ep,1,j0)=r;

    % Estimation correction - exponential distribution
    r0=-log(1-r); % k=1
    gamma(ep,2,j0)=r0;

    % Estimation correction - Gamma(k,b) distribution correction
    X0=Inc(i:w+i-1); X0(1)=I(i);
    if d==0
        b0=GammaEst(X0,r,w,alpha);
        gamma(ep,3,j0)=b0/alpha;
    else
        for k=2:n
            bm=GammaEst(X0,r,w,alpha(k));
            gamma(ep,k+1,j0)=bm/alpha(k);
            b0(ep,k+1,j0)=bm;
        end
    end
    dI(ep,j0) = 1/w*(D(w+i-1)-D(i))/(CC(w+i-1)-CC(i));
end
end

field1 = 'R2_R';  value1 = R2y1;
field2 = 'SE_R';  value2 = SEy1;
field3 = 'MSE_R';  value3 = MSEy1;
field4 = 'Mean_R';  value4 = Meany1;

OLSvals=struct(field1,value1,field2,value2,field3,value3,field4,value4);

field1 = 'KS_R';  value1 = KSy1;
field2 = 'KSp_R';  value2 = KSy1p;
field3 = 'SW_R';  value3 = SWy1;
field4 = 'SWp_R';  value4 = SWy1p;

statvals=struct(field1,value1,field2,value2,field3,value3,field4,value4);



