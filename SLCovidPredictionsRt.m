clear; clc;
% Prediction of Rt

A=readtable('SLCovidData.xlsx');

% Datetime
t = table2array(A(:,1)); 
% M = Deaths	N	S	I	R
M = table2array(A(:,2:end));

% prediction window w, a moving window length w0 
w=15; Rt_pred=zeros(w,2);

% Moving window lenghts are w0=7,10,30
% w0=10; a=1.4; % First Wave
% w0=10; a=1.2; % Second Wave
% w0=30; a=4;    
% w0=7; a=1.0; 

% Choose the appropriate parameter set to run the predictions.
k0=1;
if k0==1
    t0=80; % First Wave
    z=[7,10,30];  a0=[1,1.4,4];
else
    t0=273; % Second Wave
    z=[7,10,30];  a0=[1,1.2,4];
end


for i=1:3
    w0=z(i); a=a0(i);
    disp('    w        a'); disp([w0 a]);

    M0=M(1:t0,:);
    disp(t(t0))
    [Rt_pred(:,1),~,~,~]=Rtpredict(M0,w,w0,a,1); % Geometric
    [Rt_pred(:,2),~,~,~]=Rtpredict(M0,w,w0,a,0); % Gamma
    % disp('  Geometric    Gamma')
    % disp(Rt_pred)

    y=find(Rt_pred(:,1)>1,1,'last');
    if y==w
        t1='none';
    else
        t1=t(t0+y);
    end      
    fprintf('The detected peak for Geometric distribution: %s\n',t1);

    y=find(Rt_pred(:,2)>1,1,'last');
    if y==w
        t2='none';
    else
        t2=t(t0+y);
    end   
    fprintf('The detected peak for Gamma distribution : %s\n\n',t2);
    
    M0=M(1:t0+3,:);
    disp(t(t0+3))
    [Rt_pred(:,1),~,~,~]=Rtpredict(M0,w,w0,a,1); % Geometric
    [Rt_pred(:,2),~,~,~]=Rtpredict(M0,w,w0,a,0); % Gamma
    % disp('  Geometric    Gamma')
    % disp(Rt_pred)

    y=find(Rt_pred(:,1)>1,1,'last');
    if y==w
        t1='none';
    else
        t1=t(t0+3+y);
    end      
    fprintf('The detected peak for Geometric distribution : %s\n',t1);

    y=find(Rt_pred(:,2)>1,1,'last');
    if y==w
        t2='none';
    else
        t2=t(t0+3+y);
    end   
    fprintf('The detected peak for Gamma distribution : %s\n\n',t2);

    disp('%%%%%%%%%%%%%%%%')
end