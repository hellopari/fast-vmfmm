close all; clear; clc;
tic
 data=load('experiment_data/4c3456.mat'); %0.8284


% �б�ǩ����
x=data(1).RandVMF(:,1:end-1);
label=data(1).RandVMF(:,end);
% ���ݵ�ķֲ�
% figure;
% scatter3(x(:,1), x(:,2), x(:,3),3, 'r', 'filled'), axis([-1 1 -1 1 -1 1]);
opts = mkopts_avdp;   
 result=vdpgm(x',label,opts); 

% ʱ��
toc







