close all; 
 gene=load('data/yeast.mat'); 
%gene=load('data/Human_Fibroblasts.mat'); 
%gene=load('data/Sporulation.mat');
x=gene.data;
Y_pred =kmeans(x,6);
save('kmeans_Sporulation.mat','Y_pred');
 s=silhouette(x, Y_pred);
 si=mean(s)
