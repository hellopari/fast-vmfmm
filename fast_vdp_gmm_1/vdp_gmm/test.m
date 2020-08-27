close all; clear; 



  gene=load('data/Human_Fibroblasts.mat'); %si=0.4507 0.1947 0.2404  0.3695 0.1735 k=3
  %gene=load('data/yeast.mat'); %si= 0.3709  k=5
  %gene=load('data/Sporulation.mat'); %si=0.3439  k=5
data=gene.data;

 %opts = mkopts_avdp;
 opts =mkopts_bj(5);
result=vdpgm(data',opts);
t2=cputime; 








