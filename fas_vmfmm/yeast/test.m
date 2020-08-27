close all; clear; clc;
tic
gene=load('data/yeast.mat'); % si=  0.5824
x=gene.data;
data=bsxfun(@rdivide,x,sqrt(sum(power(x,2),2)));

opts = mkopts_avdp;  
 result=vdpgm(data',opts); 
 toc
 K=result.K-1