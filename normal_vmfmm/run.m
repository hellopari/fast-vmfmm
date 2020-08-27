% Example code to run the Variational Bayesian Clustering toolbox
%
% By: Shreyas Seshadri (shreyas.sesahdri@aalto.fi), Ulpu Remes and Okko Rasaen, Last update:19.10.2016
% (C) MIT license
% For license terms and references, see README.txt

close all; clear;clc;

tic
%% replace1
% a=1;
a='yeast'  
gene=load('data/yeast.mat'); 

% a='Human_Fibroblasts'
% gene=load('data/Human_Fibroblasts.mat');
% a='Sporulation'
% gene=load('data/Sporulation.mat');
data=gene.data;
 dataOrig=data;
% figure;
% scatter3(data(:,1), data(:,2), data(:,3),3, 'r', 'filled'), axis([-1 1 -1 1 -1 1]);
%% replace2
weightopt = 1; %1=DD 2=DP, 3=PYP 
if weightopt==1
    'DD'
else
    'DP'
end
modelopt = 4; %1='GMM-Full', 2='GMM-Diag', 3='GMM-Fixed' or 4='VMM'

% Weight type (DD, DP and PYP)
if weightopt == 1
    op.Pi_Type = 'DD'; % 'DP', 'DD' or 'PYP'(include prior.g)
    prior.alpha = 1;
elseif weightopt == 2
    op.Pi_Type = 'DP'; % 'DP', 'DD' or 'PYP'(include prior.g)
    prior.alpha = 1;
elseif weightopt == 3
    op.Pi_Type = 'PYP'; % 'DP', 'DD' or 'PYP'(include prior.g)    
    prior.g=0.5;
    prior.alpha = 1;
end

% model type (GMM (cov - full, diagonal or Fixed) or VMM)
if modelopt == 1
    op.model_Type = 'GMM-Full'; %  'GMM-Full', 'GMM-Diag', 'GMM-Fixed' or 'VMM'(check appropriate priors required)
    prior.m0 = mean(data,1);
    prior.W0 = inv(cov(data)); %prior cov of the gaussian used to model the means
    prior.beta0 = 1;
    prior.v0 = D+2;
elseif modelopt == 2
    op.model_Type = 'GMM-Diag'; %  'GMM-Full', 'GMM-Diag', 'GMM-Fixed' or 'VMM'(check appropriate priors required)
    prior.m0 = mean(data,1);
    prior.b0 = diag(inv(cov(data)));
    prior.a0 = 1;
    prior.beta0 = 1;
elseif modelopt == 3
    op.model_Type = 'GMM-Fixed'; %  'GMM-Full', 'GMM-Diag', 'GMM-Fixed' or 'VMM'(check appropriate priors required)
    prior.m0 = mean(data,1);
    prior.Pres0 = (1/0.05)*eye(D);%inv(cov(data)); %prior cov of the gaussian used to model the means
    op.PresMain = (1/1e-3)*eye(D);%inv(cov(data)); % known covariance of the GMMs
elseif modelopt == 4
    op.model_Type = 'VMM'; %  'GMM-Full', 'GMM-Diag', 'GMM-Fixed' or 'VMM'(check appropriate priors required)
    % IMPORTANT- UNIT NORM DATA FOR VMM
    data=bsxfun(@rdivide,data,sqrt(sum(power(data,2),2))); % normalise to unit length     
    prior.mu=sum(data)/norm(sum(data)); % assume data points close to mean
     prior.beta=0.05; % do not trust mean parameter too much
    % concentration parameter prior:
    if strcmp(a,'yeast')
        prior.a=0.1; 
        prior.b=0.05; 
    elseif strcmp(a,'Human_Fibroblasts')
        prior.a=0.4; 
        prior.b=0.6; 
    elseif strcmp(a,'Sporulation')
        prior.a=0.7; 
        prior.b=0.7; 
    else
        prior.a=1; % assume unconcentrated data
        prior.b=0.01; 
    end
end


% general VB priors
op.init_Type = 'random'; % initilization of clusters -  'random' or 'self'
op.repeats = 1;% how many repaeats
%replace3
op.K =6;% truncation limit
op.stopCrit = 'freeEnergy'; %'number' of runs or 'freeEnergy'
op.freethresh = 1e-6;
op.max_num_iter =700; % max number of runs
%replace4
op.reorder = 0; % roderder (usualy u sed for NP(DP and PYP) methods)
seed =1;

result = VB_mixModel(data,prior,op,seed);

foundClusters = result.z(1); %are the found cluters

y1=foundClusters{1};
k=length(unique(y1))
 s=silhouette(data, foundClusters{1});
 si=mean(s)
 toc

 


% if D==2
%     plotClustering(dataOrig,foundClusters{1},'Found Clusters');
% end
% if D==3
%     plotClustering3d(dataOrig,foundClusters{1},'Found Clusters');
% end
% [acp,asp,acc]=acp_asp_acc(foundClusters{1},z)

% posterior probablities 
postProbs = updateR(data,result.model.post,op,result.model.extra);
