close all; clear; tic
data=load('experiment_data/2c58.mat');  %0.8351

% �б�ǩ����
x=data(1).RandVMF(:,1:end-1);
label=data(1).RandVMF(:,end);
% ���ݵ�ķֲ�
% figure;
% scatter3(x(:,1), x(:,2), x(:,3),3, 'r', 'filled'), axis([-1 1 -1 1 -1 1]);

opts = mkopts_avdp;  
 result=vdpgm(x',label,opts); 
toc






