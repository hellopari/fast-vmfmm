function plotClustering( data,z,titleStr )

figure;
colors = {'k','b','r','g'};
zVals = unique(z)';
plotCount=0;
for k = zVals
    plotCount = plotCount+1;
    inds = (z == k);
     if plotCount>length(colors)
        colors = [colors rand(3,1)];
    end
    scatter3(data(inds,1),data(inds,2),data(inds,3),2,colors{plotCount});hold on;
end
hold off;
title(titleStr);

end

