close all; clear; 



  %gene=load('data/Human_Fibroblasts.mat'); %si=0.0876  k=3
 gene=load('data/yeast.mat'); %si=0.6234  k=3
   %gene=load('data/Sporulation.mat'); %si=0.3919  k=4
x=gene.data;
 data=bsxfun(@rdivide,x,sqrt(sum(power(x,2),2)));

 opts = mkopts_avdp;
% opts =mkopts_bj(6);
result=vdpgm(data',opts);
t2=cputime; 








