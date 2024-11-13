clear; clc;

A=readtable('SLCovidData.xlsx');

A=A(209:372,:);

% Datetime
t = table2array(A(:,1)); 
% M = Deaths	N	S	I	R
M = table2array(A(:,2:end));
I=M(:,4); p =length(t); t0=t(1:end-1);

% Window and Gamma(a,b)
w0=[7,10,30];  a0=[1,1.2,4];
for i=1:3
    w=w0(i); a=a0(i);
    [beta,stats_S,OLS_S]=TransmitEstimate(M,w);
    [gamma,dI,b0,stats_R,OLS_R]=RecoveryEstimate(M,w,a,0);
    Rt=beta./(dI+gamma);

    disp('    w          b0'); disp([w b0]);

    fprintf('Local peak: %s\n',t(31));
    t1=t0(29+find(Rt(30:end,1)<1,1));
    fprintf('The detected peak for Geometric distribution: %s\n',t1-1);
    t2=t0(29+find(Rt(30:end,end)<1,1));
    fprintf('The detected peak for Gamma distribution: %s\n\n',t2-1);
    
    fprintf('Rt = %1.2f at: %s\n',Rt(65,3),t(65))
    fprintf('Rt = %1.2f at: %s\n\n',Rt(68,3),t(68))


    fprintf('Local peak: %s\n',t(72));
    t1=t0(69+find(Rt(70:end,1)<1,1));
    fprintf('The detected peak for Geometric distribution: %s\n',t1-1);
    t2=t0(69+find(Rt(70:end,end)<1,1));
    fprintf('The detected peak for Gamma distribution: %s\n\n',t2-1);

    fprintf('Local peak: %s\n',t(114));
    t1=t0(112+find(Rt(113:end,1)<1,1));
    fprintf('The detected peak for Geometric distribution: %s\n',t1-1);
    t2=t0(112+find(Rt(113:end,end)<1,1));
    fprintf('The detected peak for Gamma distribution: %s\n\n',t2-1);

    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

end

% Choose Window and Gamma(a,b) for plot
%w0=[10,30];  a0=[1.2,4];
w0=[7,30];  a0=[1,4];

f1 = figure;
f2 = figure;

for i=1:2
    w=w0(i); a=a0(i);
    [beta,stats_S,OLS_S]=TransmitEstimate(M,w);
    [gamma,dI,b0,stats_R,OLS_R]=RecoveryEstimate(M,w,a,0);
    Rt=beta./(dI+gamma);

    set(0, 'CurrentFigure', f1)
    subplot(3,2,i)
    plot(t0,Rt(:,end,1),'k',t,ones(p,1),'r'); hold on
    plot(t0,Rt(:,1,1),'b'); hold off
    if i==1
        ylim([0 6]);
    else
        ylim([0 4]);
    end
    if i==2
        dl = datetime('01-Nov-2020');
    else
        dl = datetime('01-OCt-2020');
    end
    dr = datetime('01-Apr-2021');
    xlim([dl dr])

    subplot(3,2,i+2)
    plot(t0,beta,'k');
    if i==2
        dl = datetime('01-Nov-2020');
    else
        dl = datetime('01-OCt-2020');
    end
    dr = datetime('01-Apr-2021');
    xlim([dl dr])


    subplot(3,2,i+4)
    plot(t0,gamma(:,1,1),'b'); hold on
    plot(t0,gamma(:,end,1),'k'); 
    if i==2
        dl = datetime('01-Nov-2020');
    else
        dl = datetime('01-OCt-2020');
    end
    dr = datetime('01-Apr-2021');
    xlim([dl dr])


    set(0, 'CurrentFigure', f2)
    subplot(2,4,i)
    plot(t0,OLS_S.R2_S)
    if i==2
        dl = datetime('01-Nov-2020');
    else
        dl = datetime('01-OCt-2020');
    end
    dr = datetime('01-Apr-2021');
    xlim([dl dr])

    subplot(2,4,i+2)
    plot(t0,OLS_R.R2_R)
    if i==2
        dl = datetime('01-Nov-2020');
    else
        dl = datetime('01-OCt-2020');
    end
    dr = datetime('01-Apr-2021');
    xlim([dl dr])

    subplot(2,4,i+4)
    if i==1
        plot(t0,stats_S.SW_S)
        dl = datetime('01-Oct-2020');
    else
        plot(t0,stats_S.KS_S)
        dl = datetime('01-Nov-2020');
    end
    dr = datetime('01-Apr-2021');
    xlim([dl dr])

    subplot(2,4,i+6)
    if i==1
        plot(t0,stats_R.SW_R)
        dl = datetime('01-Oct-2020');
    else
        plot(t0,stats_R.KS_R)
        dl = datetime('01-Nov-2020');
    end
    dr = datetime('01-Apr-2021');
    xlim([dl dr])

end
