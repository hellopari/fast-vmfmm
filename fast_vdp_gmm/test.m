close all; clear; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 代码修改:
% 2:数据归一化
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 gene=load('data/Human_Fibroblasts.mat'); %si=0.4102  k=3  hp_prior.xi0 = 0.05;
 % gene=load('data/yeast.mat'); %si=0.5116  k= 5  hp_prior.xi0 = 0.04;
 %gene=load('data/Sporulation.mat'); %si=0.5442  k=5  hp_prior.xi0 = 0.1;
 x=gene.data;
 data=bsxfun(@rdivide,x,sqrt(sum(power(x,2),2)));

 opts = mkopts_avdp;
result=vdpgm(data',opts);
t2=cputime; 








