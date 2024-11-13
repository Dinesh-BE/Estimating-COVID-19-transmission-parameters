clear; clc;
% Prediction of Cases

A=readtable('SLCovidData.xlsx');

% Datetime
t = table2array(A(:,1)); 
% M = Deaths	N	S	I	R
M = table2array(A(:,2:end));

% prediction window w, a moving window length w0 
w=30; Inc=zeros(w,2); New_C=zeros(w,2); Prev=zeros(w,2);


% Moving window lenghts are w0=7,10,30
% Choose the appropriate parameter set to run the predictions.
% w0=10; a=1.4; % First Wave
% w0=10; a=1.2; % Second Wave
% w0=30; a=4;    
% w0=7; a=1.0; 

% Choose the data 
t0=80; % First Wave
% t0=273; % Second Wave
fprintf('Prediction date is %s\n',t(t0+7));
fprintf('Prevalence cases at peak %d\n',M(t0+7,4));
ActInc=M(t0+7,4)+M(t0+7,5)+M(t0+7,1)...
    -(M(t0+6,4)+M(t0+6,5)+M(t0+6,1));
fprintf('Incidence cases at peak  %d\n\n',ActInc);

w0=30; a=4;

M0=M(1:t0,:);
fprintf('Data Availability %s\n',t(t0));
disp(t(t0))
[~,Prev(:,1),New_C(:,1),Inc(:,1)]=Rtpredict(M0,w,w0,a,1); % Geometric
[~,Prev(:,2),New_C(:,2),Inc(:,2)]=Rtpredict(M0,w,w0,a,0); % Gamma
disp('  Total        Prevelence    Incidence' )
disp('Geome  Gamma  Geome  Gamma  Geome  Gamma')
disp(round([New_C(7,:) Prev(7,:) Inc(7,:)]))

ActInc=M(t0+2:t0+w+1,4)+M(t0+2:t0+w+1,5)+M(t0+2:t0+w+1,1)...
    -(M(t0+1:t0+w,4)+M(t0+1:t0+w,5)+M(t0+1:t0+w,1));
%plot(t(t0+1:t0+w),ActInc,t(t0+1:t0+w),Inc)
disp('%%%%%%%%%%%%%%%%')
M0=M(1:t0+3,:);
fprintf('Data Availability %s\n',t(t0+3));
[~,Prev(:,1),New_C(:,1),Inc(:,1)]=Rtpredict(M0,w,w0,a,1); % Geometric
[~,Prev(:,2),New_C(:,2),Inc(:,2)]=Rtpredict(M0,w,w0,a,0); % Gamma
disp('  Total        Prevelence    Incidence' )
disp('Geome  Gamma  Geome  Gamma  Geome  Gamma')
disp(round([New_C(4,:) Prev(4,:) Inc(4,:)]))
disp('%%%%%%%%%%%%%%%%')
M0=M(1:t0+5,:);
fprintf('Data Availability %s\n',t(t0+5));
[~,Prev(:,1),New_C(:,1),Inc(:,1)]=Rtpredict(M0,w,w0,a,1); % Geometric
[~,Prev(:,2),New_C(:,2),Inc(:,2)]=Rtpredict(M0,w,w0,a,0); % Gamma
disp('  Total        Prevelence    Incidence' )
disp('Geome  Gamma  Geome  Gamma  Geome  Gamma')
disp(round([New_C(2,:) Prev(2,:) Inc(2,:)]))