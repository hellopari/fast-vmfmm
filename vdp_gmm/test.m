close all; clear; 



  %gene=load('data/Human_Fibroblasts.mat'); 
  gene=load('data/yeast.mat');
  %gene=load('data/Sporulation.mat');
data=gene.data;

% GMModel = fitgmdist(data,4);
% T1 = cluster(GMModel,data);
%  s=silhouette(data, T1);
%   si=mean(s)

 opts =mkopts_bj(6);
result=vdpgm(data',opts);









