close all; clear; clc;
gene= importdata('data/Human_Fibroblasts.mat');
x=gene;

clusters=load('Y_pred.mat');

pred=clusters.Y_pred;

figure 
for c = 1:length(unique(pred))-1
    subplot(3,3,c);
    plot(x((pred == c),:)');
    xlabel('Time');
    ylabel('Expression Level');
   title(c);
    axis tight
   
end

suptitle('Hierarchical Clustering of Profiles');