close all; clear; clc;
gene= importdata('data/Sporulation.mat');
x=gene;

clusters=load('kmeans_Sporulation.mat');

pred=clusters.Y_pred;
k=length(unique(pred));

figure 
for c = 1:k
    subplot(3,3,c);
    plot(x((pred == c),:)');
    xlabel('Time');
    ylabel('Expression Level');
   title(c);
    axis tight
   
end

suptitle('Hierarchical Clustering of Profiles');