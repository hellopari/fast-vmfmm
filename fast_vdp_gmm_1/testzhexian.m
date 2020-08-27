close all; clear; clc;
gene= importdata('data/yeast.mat');
x=gene.data;

clusters=load('Y_pred_yeast_gmm.mat');

pred=clusters.Y_pred;
t=length(unique(pred))
figure 
for c = 1:length(unique(pred))
    subplot(3,3,c);
    plot(x((pred == c),:)');
    xlabel('Time');
    ylabel('Expression Level');
   title(c);
    axis tight
   
end

suptitle('Hierarchical Clustering of Profiles');