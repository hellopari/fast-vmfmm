close all; clear; clc;tic
gene=load('data/Sporulation.mat');  %si =-0.1260
x=gene.data;
data=bsxfun(@rdivide,x,sqrt(sum(power(x,2),2)));

opts = mkopts_avdp;  

 result=vdpgm(data',opts);
 toc
 K=result.K-1
 si=result.si


% v =1:1:10;
% C = nchoosek(v,2);
% 
% for i=1:length(C)
%      result(i)=vdpgm(data',opts,C(i,:)); 
% end
% 
% s=[result(:).si];
% k=[result(:).K];
% a=s+k;
% [x,y] = max(a,[],2)
% b=1;

