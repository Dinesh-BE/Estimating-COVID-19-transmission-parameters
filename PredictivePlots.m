clear; clc;
% Prediction of Cases

A=readtable('SLCovidData.xlsx');
warning('off','MATLAB:nearlySingularMatrix')
warning('off','MATLAB:rankDeficientMatrix')


% Datetime
t = table2array(A(:,1)); 
% M = Deaths	N	S	I	R
M = table2array(A(:,2:end));

% prediction window w, a moving window length w0 
w=15; Inc=zeros(w,2); New_C=zeros(w,2); Prev=zeros(w,2);
w0=30; a=4;

% Choose the data 
t0=80; % First Wave
%t0=273; % Second Wave

PrevVAR3= zeros(w,3);
IncVAR3 = zeros(w,3);
k=1;


for j=[0,3,5]
    A0=A(1:t0+j,:);
    y=[A0.S A0.I A0.R];
    
    
    Mdl = vecm(3,2,3);
    EstMdl = estimate(Mdl,y);
    [Y,~] = forecast(EstMdl,w+1,y);
    PrevVAR3(:,k)=Y(2:end,2);
    IncVAR3(:,k)=Y(2:end,2)+Y(2:end,3) -Y(1:end-1,2)-Y(1:end-1,3);
    k=k+1;
end

PrevVAR2= zeros(w,3);
IncVAR2 = zeros(w,3);
k=1;

for j=[0,3,5]
    A0=A(1:t0+j,:);
    y=[A0.S A0.I A0.R];
    
    
    Mdl = vecm(3,2,2);
    EstMdl = estimate(Mdl,y);
    [Y,~] = forecast(EstMdl,w+1,y);
    PrevVAR2(:,k)=Y(2:end,2);
    IncVAR2(:,k)=Y(2:end,2)+Y(2:end,3) -Y(1:end-1,2)-Y(1:end-1,3);
    k=k+1;
end


M0=M(1:t0,:);
fprintf('Plot1 - Data Availability %s\n',t(t0));
[~,Prev(:,1),New_C(:,1),Inc(:,1)]=Rtpredict(M0,w,w0,a,1); % Geometric
[~,Prev(:,2),New_C(:,2),Inc(:,2)]=Rtpredict(M0,w,w0,a,0); % Gamma
ActInc=M(t0+1:t0+w,4)+M(t0+1:t0+w,5)+M(t0+1:t0+w,1)...
    -(M(t0:t0+w-1,4)+M(t0:t0+w-1,5)+M(t0:t0+w-1,1));
ActInc=ActInc.*(ActInc>0);
Inc=Inc.*(Inc>0);
s=t0+1:t0+w;

subplot(2,3,1)
plot(t(s),ActInc,t(s),Inc,...
    t(s),IncVAR3(:,1), t(s),IncVAR2(:,1))
legend('Actual', 'Geometric','Gamma','VAR(3,2)','VAR(2,2)')
subplot(2,3,4)
plot(t(s),cumsum(ActInc),t(s),cumsum(Inc),...
    t(s),cumsum(IncVAR3(:,1)), t(s),cumsum(IncVAR2(:,1)))
%legend('Actual', 'Geometric','Gamma','VAR(3,2)','VAR(2,2)')
disp('%%%%%%%%%%%%%%%%')
M0=M(1:t0+3,:);
fprintf('Plot2 - Data Availability %s\n',t(t0+3));
[~,Prev(:,1),New_C(:,1),Inc(:,1)]=Rtpredict(M0,w,w0,a,1); % Geometric
[~,Prev(:,2),New_C(:,2),Inc(:,2)]=Rtpredict(M0,w,w0,a,0); % Gamma
ActInc=M(t0+4:t0+w+3,4)+M(t0+4:t0+w+3,5)+M(t0+4:t0+w+3,1)...
    -(M(t0+3:t0+w+2,4)+M(t0+3:t0+w+2,5)+M(t0+3:t0+w+2,1));
ActInc=ActInc.*(ActInc>0);
Inc=Inc.*(Inc>0);
s=t0+4:t0+w+3;

subplot(2,3,2)
plot(t(s),ActInc,t(s),Inc,...
    t(s),IncVAR3(:,2), t(s),IncVAR2(:,2))
%legend('Actual', 'Geometric','Gamma','VAR(3,2)','VAR(2,2)')

subplot(2,3,5)
plot(t(s),cumsum(ActInc),t(s),cumsum(Inc),...
    t(s),cumsum(IncVAR3(:,2)), t(s),cumsum(IncVAR2(:,2)))

disp('%%%%%%%%%%%%%%%%')
M0=M(1:t0+5,:);
fprintf('Plot3 - Data Availability %s\n',t(t0+5));
[~,Prev(:,1),New_C(:,1),Inc(:,1)]=Rtpredict(M0,w,w0,a,1); % Geometric
[~,Prev(:,2),New_C(:,2),Inc(:,2)]=Rtpredict(M0,w,w0,a,0); % Gamma
ActInc=M(t0+6:t0+w+5,4)+M(t0+6:t0+w+5,5)+M(t0+6:t0+w+5,1)...
    -(M(t0+5:t0+w+4,4)+M(t0+5:t0+w+4,5)+M(t0+5:t0+w+4,1));
ActInc=ActInc.*(ActInc>0);
Inc=Inc.*(Inc>0);
s=t0+6:t0+w+5;
subplot(2,3,3)
plot(t(s),ActInc,t(s),Inc,...
    t(s),IncVAR3(:,3), t(s),IncVAR2(:,3))
%legend('Actual', 'Geometric','Gamma','VAR(3,2)','VAR(2,2)')
subplot(2,3,6)
plot(t(s),cumsum(ActInc),t(s),cumsum(Inc),...
    t(s),cumsum(IncVAR3(:,3)), t(s),cumsum(IncVAR2(:,3)))


