close all; clear; clc;
t1=cputime;
 data=load('experiment_data/4c3456.mat'); %0.8284


% 有标签数据
x=data(1).RandVMF(:,1:end-1);
label=data(1).RandVMF(:,end);
% 数据点的分布
% figure;
% scatter3(x(:,1), x(:,2), x(:,3),3, 'r', 'filled'), axis([-1 1 -1 1 -1 1]);
opts = mkopts_avdp;   
 result=vdpgm(x',label,opts); 
 t2=cputime; 
% 时间
t=t2-t1







