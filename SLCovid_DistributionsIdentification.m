clear; clc;

A=readtable('SLCovidData.xlsx');

% Datetime
t = table2array(A(:,1)); 
% M = Deaths	N	S	I	R
M = table2array(A(:,2:end));
I=M(:,4); p =length(t); t0=t(1:end-1);

f1 = figure;

a0=[3,6]; w0=[10,30]; 
for i=1:2
    w=w0(i); a=a0(i); C=linspace(1,a,41);
    [beta,stats_S,OLS_S]=TransmitEstimate(M,w);
    [gamma,dI,b0,stats_R,OLS_R]=RecoveryEstimate(M,w,a,1);

    Rt=beta./(dI+gamma);
    
    set(0, 'CurrentFigure', f1)
    subplot(2,2,2*i-1)
    plot(C,Rt(58,2:end)); hold on;
    plot(C,Rt(59,2:end)); hold on
    plot(C,Rt(60,2:end)); hold on;
    plot(C,ones(length(C),1)); hold off    
    legend("May 7","May 8", "May 9")

    subplot(2,2,2*i)
    plot(C,Rt(239,2:end)); hold on;
    plot(C,Rt(240,2:end)); hold on
    plot(C,Rt(241,2:end)); hold on;
    plot(C,ones(length(C),1)); hold off    
    legend("Nov 4","Nov 5", "Nov 6")

end