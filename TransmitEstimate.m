function [beta,statvals,OLSvals]=TransmitEstimate(M,w0)
% Transmission parameter estimation
% M - Data matrix,  w0 - window sizes, 


% Number of data points
L=length(M(:,1));

% M = Deaths	N	S	I	R	
N=M(:,2); S=M(:,3); I=M(:,4);

% Forward differencing
dS=S(2:end)-S(1:end-1); 
%Backward differencing
%dS=[dS(2:end);0];

S=S(1:end-1); I=I(1:end-1); N=N(1:end-1);
p=L-1; q=length(w0);

beta = zeros(p,q);

% OLS results:
% Uncentered R squared
R2y0=zeros(p,q); 
% Standard errors
SEy0=zeros(p,q);
% MSE of residuals
MSEy0=zeros(p,q);
% Mean of residuals
Meany0=zeros(p,q);

% Kolmogorov-Smirnov test results
% Null hypothesis test is that the data comes from a standard normal distribution
KSy0=zeros(p,q);
KSy0p=zeros(p,q);

% Shapiro-Wilk parametric hypothesis test results
% Null hypothesis test is that the data comes from a normal distribution
SWy0=zeros(p,q); SWy0p=zeros(p,q);

for j0=1:q

w=w0(j0);

for i=1:p-w+1
    ep = w+i-1; % estimation point
    
    % dS
    X = (S(i:w+i-1).*I(i:w+i-1)./N(i:w+i-1)).^0.5;
    y = dS(i:w+i-1)./X; 
    Smdl = fitlm(X,y,'Intercept',false);
    beta(ep,j0) = -Smdl.Coefficients.Estimate;

    % Regression quality
    SEy0(ep,j0)   = Smdl.Coefficients.SE;
    MSEy0(ep,j0)  = Smdl.RMSE;
    ResS=table2array(Smdl.Residuals(:,1));
    z0= ResS(2:end)-ResS(1:end-1); 
    z0 = z0(~isnan(z0));
    Meany0(ep,j0) = mean(z0);
    R2y0(ep,j0) =1-(ResS'*ResS)/sum((y).^2); % Uncentered

    % Kolmogorov-Smirnov and Shapiro-Wilk tests
    if std(z0) ==0
        z0=0;
    else
        z0=z0./std(z0);
    end
    
    if length(z0)>2
        [KSy0(ep,j0),KSy0p(ep,j0)]=kstest(z0);
    else
        KSy0(ep,j0)=0;
        KSy0p(ep,j0)=1;
    end    
    %[KSy0(ep,j0),KSy0p(ep,j0)]=kstest(z0);

    z0=z0-mean(z0);

    if length(z0)>2
        [SWy0(ep,j0),SWy0p(ep,j0), ~] = swtest(z0, 0.05);
    else
        SWy0(ep,j0)=0;
        SWy0p(ep,j0)=1;
    end

end
end

field1 = 'R2_S';  value1 = R2y0;
field2 = 'SE_S';  value2 = SEy0;
field3 = 'MSE_S';  value3 = MSEy0;
field4 = 'Mean_S';  value4 = Meany0;

OLSvals=struct(field1,value1,field2,value2,field3,value3,field4,value4);

field1 = 'KS_S';  value1 = KSy0;
field2 = 'KSp_S';  value2 = KSy0p;
field3 = 'SW_S';  value3 = SWy0;
field4 = 'SWp_S';  value4 = SWy0p;


statvals=struct(field1,value1,field2,value2,field3,value3,field4,value4);