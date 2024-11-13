function [Rt,I0,New_C,Inc]=Rtpredict(M0,w,w0,a,Gm)

% Inputs:
% w - Prediction Window
% w0 - moving window of length for parameter estimation
% a - Shape parameter for the Gamma distribution
% If Gm = 0, Gamma distribution
% If Gm =1, Geomeric distribution
% If Gm =2, Exp distribution

% Outputs
% Rt - effective reproduction number
% I0 - Prevalence cases
% New_C - Commulative new casses
% Inc - New incidence


% Number of data points
L=length(M0(:,1));
Rt=zeros(w,1); S0=zeros(w,1); I0=zeros(w,1); R0=zeros(w,1);
CC=zeros(w,1); New_C=zeros(w,1); Inc=zeros(w,1);

M=M0(end-w0:end,:);

[beta,~,~]=TransmitEstimate(M,w0);  
[gamma,dI,~,~,~]=RecoveryEstimate(M,w0,a,0);
% Rt_p=beta./(dI+gamma); % [Geometric, Exponential, Gamma]
% disp([Rt_p(end,1) Rt_p(end,3)]) % uncomment to see Rt

for j=1:w
    D=M(end,1);
    N=M(end,2);
    S=M(end,3); 
    I=M(end,4);
    R=M(end,5);
    if j==1
        CC_p=I+R+D;
        Inc(j)=CC_p-M(end-1,4)-M(end-1,5)-M(end-1,1);
    end

    S=round(S-beta(end)*S*I/N);
    if Gm==1
        R=round(R+gamma(end,1)*I); % Geometric
    elseif Gm==2
        R=round(R+gamma(end,2)*I); % Exp
    else
        R=round(R+gamma(end,end)*I); % Gamma
    end
    N=round(N-dI(end)*I);
    I=N-R-S;
    if I<0
        I=0;
    end
    % [N;S;I;R]

    D=D+dI(end)*I;
    M=[M(2:end,:);D N S I R];
    M=[M(2:end,:);D N S I R];
    if w0<30
        [gamma,~,~,~,~]=RecoveryEstimate(M,w0,a,0);
    else
        [gamma,dI,~,~,~]=RecoveryEstimate(M,w0,a,0);
    end
    [beta,~,~]=TransmitEstimate(M,w0);
    S0(j)=S; I0(j)=I; R0(j)=R;


    if Gm==1
        Rt(j)=beta(end)/(dI(end)+gamma(end,1)); % Geometric
    elseif Gm==2
        Rt(j)=beta(end)/(dI(end)+gamma(end,2)); % Exp
    else
        Rt(j)=beta(end)/(dI(end)+gamma(end,end)); % Gamma
    end
    CC(j)=I+R+D;
    New_C(j)=CC(j)-CC_p;
    if j>1
        Inc(j)=CC(j)-CC(j-1);
    end
end



