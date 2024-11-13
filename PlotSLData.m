clear; clc;

A=readtable('SLCovidData.xlsx');

A=A(1:383,:);

% Datetime
t = table2array(A(:,1)); 
% M = Deaths	N	S	I	R
M = table2array(A(:,2:end));
I=M(:,4); D=M(:,1); R=M(:,5);
Inc=I(2:end)+R(2:end)+D(2:end)-(I(1:end-1)+R(1:end-1)+D(1:end-1));
Inc=[0;Inc];
p =length(t);

figure;
subplot(3,1,1)
plot(t,I)
subplot(3,1,2)
plot(t,Inc)
subplot(3,1,3)
plot(t,cumsum(Inc))