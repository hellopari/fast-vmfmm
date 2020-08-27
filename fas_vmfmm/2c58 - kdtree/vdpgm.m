function [results] = vdpgm(given_data,label,opts)
start_time = clock;
if nargin == 1
  opts = struct();
end
if issparse(given_data)
  given_data = full(given_data);
end
if ~ isfield(opts, 'algorithm')
  opts.algorithm = 'vdp';
end
if ~ isfield(opts, 'collapsed_means')
  opts.collapsed_means = 0;
end
if ~ isfield(opts, 'do_sort')
  opts.do_sort = '0';
end
if ~ isfield(opts, 'get_r')
  opts.get_r = 0;
end
if ~ isfield(opts, 'weight')
  opts.weight = 1;
end
if ~ isfield(opts, 'get_log_likelihood')
  opts.get_log_likelihood = 0;
end
if ~ isfield(opts, 'use_kd_tree')
  opts.use_kd_tree = 1;
end
if ~ isfield(opts, 'threshold')
  opts.threshold = 1.0e-5;
end
if ~ isfield(opts, 'sis')
  opts.sis = 0;
end
if ~ isfield(opts, 'initial_depth')
  opts.initial_depth = 3;
end
if ~ isfield(opts, 'initial_K')
  opts.initial_K = 1;
end
if ~ isfield(opts, 'ite')
  opts.ite = inf;
end
if ~ isfield(opts, 'do_split')
  opts.do_split = 0;
end
if ~ isfield(opts, 'do_merge')
  opts.do_merge = 0;
end
if ~ isfield(opts, 'do_greedy')
  opts.do_greedy = 1;
end
if ~ isfield(opts, 'max_target_ratio')
  opts.max_target_ratio = 0.5;
end
if ~ isfield(opts, 'init_of_split')
  % 'pc', 'rnd', 'rnd_close' or 'close_f'
  opts.init_of_split = 'pc';
end
if ~ isfield(opts, 'recursive_expanding_depth')
  opts.recursive_expanding_depth = 2;
end
if ~ isfield(opts, 'recursive_expanding_threshold')
  opts.recursive_expanding_threshold = 1.0e-1;
end
if ~ isfield(opts, 'recursive_expanding_frequency')
  opts.recursive_expanding_frequency = 3;
end
if isfield(opts, 'seed')
  rand('state', opts.seed);
else
  seed = rand('state');
  results.seed = seed;
end

data.given_data = given_data;
if opts.use_kd_tree
  partitions = init_kdtree_partitions(given_data, opts.initial_depth);
  data.kdtree = partitions;
 
end

% the hyperparameters of priors
hp_prior = mk_hp_prior(data, opts);


if strcmp(opts.algorithm,'vdp')
    assert(min(sum(power(data.given_data',2),2))>1-exp(-20)); % check unit norm
    assert(max(sum(power(data.given_data',2),2))<1+exp(-20));
    opts.max_kk = round(size(data.given_data',2)/2); 
    while isfinite(besseli(size(data.given_data',2)/2,opts.max_kk+10))
        opts.max_kk = opts.max_kk+10;
    end
    extra_V.kk = [30,30];
    extra_V.ln_kk = psi(hp_prior.a)-log(hp_prior.b);
    extra_V.kk_approx=extra_V.kk-(hp_prior.a>1)*(1./hp_prior.b);

end

if opts.sis > 0
  r = sequential_importance_sampling(data, hp_prior, opts);
elseif isfield(opts, 'r')
  r = opts.r;
else
  r = rand_r(data, opts.initial_K, opts);
end

[hp_posterior,extra_V] = mk_hp_posterior(data, r, hp_prior, opts,extra_V);
if opts.do_greedy
  [free_energy, hp_posterior, data,new_extra_V] = greedy(data, hp_posterior, hp_prior, opts,extra_V,label);
else
  [free_energy, hp_posterior, data] = split_merge(data, hp_posterior, hp_prior, opts);
end
kk=hp_posterior.a./hp_posterior.b
results.algorithm = opts.algorithm;
results.elapsed_time = etime(clock, start_time);
results.free_energy = free_energy;
results.hp_prior = hp_prior;
results.hp_posterior = hp_posterior;
results.K = length(hp_posterior.beta);
results.opts = opts;


Y = label;
cluster=1:size(hp_posterior.mu,2);
t1=mk_map_log_p_of_x_given_c(given_data, cluster, hp_posterior, opts,new_extra_V);
[~,Y_pred] = max(t1',[],2);

[FMeasure,Accuracy] = Fmeasure(Y', Y_pred')
 Nmi=nmi(Y', Y_pred')

if opts.get_r
  results.r = mk_r(data, hp_posterior, hp_prior, opts);
  if opts.use_kd_tree
    results.r = mk_non_kdtree_r(data, results.r);
  end
end
if opts.get_log_likelihood
  results.log_likelihood = mk_log_likelihood(data, hp_posterior, hp_prior, opts);
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [free_energy, hp_posterior, data,extra_V] = greedy(data, hp_posterior, hp_prior, opts,extra_V,label);
free_energy= mk_free_energy(data, hp_posterior, hp_prior, opts,extra_V);
disp_status(free_energy, hp_posterior, opts);
while 1
  disp('finding the best one....')
  [new_free_energy, new_hp_posterior, new_data, c,new_extra_V] = find_best_splitting(data, ...
                                                    hp_posterior, ...
                                                    hp_prior, opts,extra_V);
  if c == -1
    break
  end

  disp(['finding the best one.... done.  component ' num2str(c) ' was split.'])

  disp_status(new_free_energy, new_hp_posterior, opts);

  [new_free_energy, new_hp_posterior, new_data,~,new_extra_V] = update_posterior2(new_data, ...
                                                    new_hp_posterior,new_extra_V, ...
                                                    hp_prior, opts, 1, opts.ite);



  if free_energy_improved(free_energy, new_free_energy, 0, opts) == 0
     break
  end

  
      
  free_energy = new_free_energy;
  hp_posterior = new_hp_posterior;
  extra_V=new_extra_V;
  data = new_data;

end
disp_status(free_energy, hp_posterior, opts);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [free_energy, hp_posterior, data, c,extra_V] = find_best_splitting(data, ...
                                                  hp_posterior, ...
                                                  hp_prior, opts,extra_V);
c_max = 10;
K = size(hp_posterior.mu, 2);
candidates = find(hp_posterior.Nc>2);
if isempty(candidates)
  free_energy = 0;
  c = -1;
  return
end
r = mk_r(data, hp_posterior, hp_prior, opts,extra_V);
new_free_energy = ones(1,max(candidates))*inf;
%%%%%%%%%%%%%%%%%%%%%
fc = mk_E_log_q_p_eta(data, hp_posterior, hp_prior, opts,extra_V); % 1*K
S_tk = mk_S_tk(data, hp_posterior, hp_prior, opts,extra_V); % N*K

%%%%%%%%%%%%%%%%%%%%%
for c = candidates(1:min(c_max, length(candidates)))
  disp(['splitting ' num2str(c) '...'])
  [new_data(c), new_r,new_extra_V, info] = split(c, data, r, hp_posterior, hp_prior, opts,extra_V, 1);
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  new_c = info.new_c;
  if new_c==0
       continue
  end
  relating_n = find(sum(new_r(:,[c new_c]),2) > 0.5);
  if isempty(relating_n)
    continue
  end
  new_K = size(new_r, 2);
  
  sub_r = new_r(relating_n, [c new_c new_K]);
  sub_extra_V.kk_approx=new_extra_V.kk_approx([c new_c new_K]);
  sub_extra_V.kk=new_extra_V.kk([c new_c new_K]);
  sub_extra_V.ln_kk=new_extra_V.ln_kk([c new_c new_K]);

  if opts.use_kd_tree
    sub_data.kdtree = new_data(c).kdtree(relating_n);
  else
    sub_data.given_data = new_data(c).given_data(:,relating_n);
  end
  [sub_hp_posteior,sub_extra_V] = mk_hp_posterior(sub_data, sub_r, hp_prior, opts,sub_extra_V);
  [sub_f, sub_hp_posteior, dummy, sub_r,sub_extra_V] = update_posterior2(sub_data, ...
                                                    sub_hp_posteior,sub_extra_V, ...
                                                    hp_prior, opts, 0, 10, 1);
  if size(sub_r,2) < 3
    continue
  end

  if length(find(sum(sub_r,1)<1.0e-10)) > 1
    continue
  end
  new_S_tk = S_tk;
  if opts.use_kd_tree
    updated_data.kdtree = new_data(c).kdtree(info.updated_I);
    new_S_tk(info.updated_I,:) = mk_S_tk(updated_data, hp_posterior, hp_prior, opts,extra_V);
  end
  sub_S_tk = mk_S_tk(new_data(c), sub_hp_posteior, hp_prior, opts,sub_extra_V);
  insert_indices = [c new_c new_K:(new_K+size(sub_r,2)-3)];
  new_S_tk(:,insert_indices) = sub_S_tk;
  new_fc = fc;
  new_fc(insert_indices) = mk_E_log_q_p_eta(sub_data, sub_hp_posteior, hp_prior, opts,sub_extra_V);
  new_free_energy(c) = mk_free_energy(new_data(c), sub_hp_posteior, hp_prior, opts,sub_extra_V, new_fc, new_S_tk);
  new_r(relating_n,:) = 0;
  new_r(relating_n,insert_indices) = sub_r;
  new_extra_V.kk(insert_indices) = sub_extra_V.kk;
  new_extra_V.kk_approx(insert_indices) = sub_extra_V.kk_approx;
  new_extra_V.ln_kk(insert_indices) = sub_extra_V.ln_kk;
  new_r_cell{c} = new_r;
  sub_extra_cell{c} = new_extra_V;
end
[free_energy, c] = min(new_free_energy);
if isinf(free_energy)
  c = -1;
  return
end
data = new_data(c);
[hp_posterior,extra_V] = mk_hp_posterior(data, new_r_cell{c}, hp_prior, opts,sub_extra_cell{c});




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function log_p_of_x_given_c = mk_map_log_p_of_x_given_c(data, clusters, hp_posterior, opts,extra_V);
% clusters: e.g [1:K]
[D,N] = size(data);
K = length(clusters);
log_p_of_x_given_c = zeros(K,N);
for i = 1:K
    c = clusters(i);
    nu=size(hp_posterior.mu,1)/2-1;
    distances = data'*hp_posterior.mu(:,c);
    t3=bsxfun(@times,distances, extra_V.kk_approx(:,c));
    t1=log(besseli(nu, extra_V.kk_approx(:,c))+exp(-700));
    log_p_of_x_given_c(i,:) =nu*log( extra_V.kk_approx(:,c))-(nu+1)*log(2*pi)-t1+t3;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [new_data, new_r,new_extra_V, info] = split(c, data, r, ...
                                               hp_posterior, hp_prior, opts,extra_V, update_kdtree)
% r: N*K
if nargin < 8
  update_kdtree = 1;
end
new_data = data;

[dummy, indices] = max(r,[],2);
target_partitions_all = find(indices==c)';
if ~ isempty(target_partitions_all)
    m = [data.kdtree(target_partitions_all).mean];
    log_p_of_m_given_c = mk_map_log_p_of_x_given_c(m, c, hp_posterior, opts,extra_V); % 1*|target_partitions|
    [dummy, I] = sort(log_p_of_m_given_c,2, 'descend');
    target_partitions = target_partitions_all(I(1:ceil(length(I)*opts.max_target_ratio)));
    [new_data, r, updated_I] = expand_all_nodes(new_data, r, hp_posterior, ...
                                                     hp_prior, target_partitions, opts,extra_V);
    info.updated_I = updated_I;
else
    info.updated_I = [];
end

arg1_data = [new_data.kdtree(:).mean];
dir = divide_by_principal_component(arg1_data, ...
                                hp_posterior.B{c}, ...
                                   hp_posterior.mu(:,c));
I = find(dir>=0);

r_c1 = zeros(size(r,1),1);
r_c2 = r(:,c);
r_c1(I) = r(I,c);
r_c2(I) = 0;
N=size(r,1);
new_r = zeros(size(r,1), size(r,2)+1);
new_r(:,[1:end-2 end]) = r;
new_r(:,c) = r_c1;
new_c = size(new_r, 2) - 1;
new_r(:,new_c) = r_c2;
info.new_c = new_c;

new_extra_V.kk_approx=ones(1,size(r,2)+1)*100;
new_extra_V.kk_approx(:,[1:end-2 end]) = extra_V.kk_approx;
new_extra_V.kk_approx(:,c)=extra_V.kk_approx(c)*(length(I)/N);
new_extra_V.kk_approx(:,new_c)=extra_V.kk_approx(c)*((N-length(I))/N);

new_extra_V.kk=ones(1,size(r,2)+1)*100;
new_extra_V.kk(:,[1:end-2 end]) = extra_V.kk;
new_extra_V.kk(:,c)=extra_V.kk(c)*(length(I)/N);
new_extra_V.kk(:,new_c)=extra_V.kk(c)*((N-length(I))/N);

new_extra_V.ln_kk=ones(1,size(r,2)+1)*4;
new_extra_V.ln_kk(:,[1:end-2 end]) = extra_V.ln_kk;
new_extra_V.ln_kk(:,c)=extra_V.ln_kk(c)*(length(I)/N);
new_extra_V.ln_kk(:,new_c)=extra_V.ln_kk(c)*((N-length(I))/N);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [new_data, new_r, updated_I] = expand_all_nodes(data, ...
                                                  r, hp_posterior, ...
                                                  hp_prior, ...
                                                  target_partitions, opts,extra_V);
if size(target_partitions, 1) ~= 1
  error('target_partitions must be a row vector.')
end
new_data = data;
for a=target_partitions
  if isempty(data.kdtree(a).children)
    children = mk_children_of_partition(data.kdtree(a));
    data.kdtree(a).children = children;
  else
    children = data.kdtree(a).children;
  end
  if length(children) == 1
    continue
  end
  new_data.kdtree(a) = children(1);
  new_data.kdtree(end+1) = children(2);
end % while a <= length(data.kdtree)
updated_I = [target_partitions length(data.kdtree)+1:length(new_data.kdtree)];
sub_data.given_data = data.given_data;
sub_data.kdtree = new_data.kdtree(updated_I);
sub_r = mk_r(sub_data, hp_posterior, hp_prior, opts,extra_V);
new_r = zeros(length(new_data.kdtree),size(r,2));
new_r(1:size(r,1),:) = r;
new_r(updated_I,:) = sub_r;
disp(['### building kd-tree done; #partition ' num2str(length(data.kdtree)) ...
      ' -> ' num2str(length(new_data.kdtree))])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [data, r] = expand_recursively_until_convergence(data, ...
                                                  r, hp_posterior, ...
                                                  hp_prior, opts,extra_V);
kdtree_size = length(data.kdtree);
for a=1:length(data.kdtree)
  [data, r] = expand_recursively_until_convergence2(data, a, ...
                                                    r, hp_posterior, ...
                                                    hp_prior,opts,extra_V, ...
                                                    opts.recursive_expanding_depth);
end
[prob, best_c] = max(r, [], 2);

disp(['### building kd-tree done; #partition ' num2str(kdtree_size) ...
      ' -> ' num2str(length(data.kdtree))])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [data, r] = expand_recursively_until_convergence2(data, a, ...
                                                  r, hp_posterior, ...
                                                  hp_prior,opts,extra_V,depth);
% r : N*K
if depth == 0
  return
end
if isempty(data.kdtree(a).children)
  children = mk_children_of_partition(data.kdtree(a));
  data.kdtree(a).children = children;
else
  children = data.kdtree(a).children;
end
if length(children) == 1
  return
end
sub_data.given_data = data.given_data;
sub_data.kdtree = children;
sub_r = mk_r(sub_data, hp_posterior, hp_prior,opts,extra_V); % 2*K
diff = sub_r - repmat(r(a,:), 2, 1);
if sum(sum(diff.*diff, 1), 2)/prod(size(sub_r)) < opts.recursive_expanding_threshold
  return
end
b = length(data.kdtree)+1;
data.kdtree(a) = children(1);
data.kdtree(b) = children(2);
r([a b],:) = sub_r;
[data, r] =  expand_recursively_until_convergence2(data, a, ...
                                                  r, hp_posterior, ...
                                                  hp_prior, opts,extra_V, depth-1);
[data, r] =  expand_recursively_until_convergence2(data, b, ...
                                                  r, hp_posterior, ...
                                                  hp_prior, opts,extra_V,depth-1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [free_energy, hp_posterior, data, r,extra_V] = update_posterior2(data, hp_posterior,extra_V,hp_prior, opts, upkdtree, ite, do_sort);
% update r: N*K
disp(['### updating posterior ...'])
free_energy = inf;
if nargin == 5
  upkdtree = 0;
end
if nargin < 7
  ite = inf;
end
if nargin < 8
  do_sort = 1;
end
i = 0;
last_Nc = 0;
start_sort = 0;
while 1
  i = i+1;
  [new_free_energy, S_tk] = mk_free_energy(data, hp_posterior, hp_prior, opts,extra_V);
   disp_status(new_free_energy, hp_posterior, opts);
  if (~isinf(ite) && i>=ite) || ...
        (isinf(ite) && free_energy_improved(free_energy, new_free_energy, 0, opts) == 0)
    free_energy = new_free_energy;
    if do_sort && opts.do_sort && ~ start_sort
      start_sort = 1;
    else
      break
    end
  end
  last_Nc = hp_posterior.Nc;
  free_energy = new_free_energy;
  [r, data] = mk_r(data, hp_posterior, hp_prior, opts,extra_V, S_tk);
  freq = opts.recursive_expanding_frequency;
  if opts.use_kd_tree & upkdtree & (freq==1 | mod(i,freq)==1)
    [data, r] = expand_recursively_until_convergence(data, ...
                                                      r, ...
                                                      hp_posterior, ...
                                                      hp_prior,opts,extra_V);
  end
  % check if the last component is small enough
  
  if isequal(opts.algorithm, 'vdp') & sum(r(:,end)) >1
    r(:,end+1) = 0;
    extra_V.kk_approx(:,end+1) = 1;
    extra_V.kk(:,end+1)=1;
    extra_V.ln_kk(:,end+1)=4; 
  end
  
  if start_sort
    [r,~,extra_V] = sort_r(data, r, opts,extra_V);
  end
 
  if isequal(opts.algorithm, 'vdp') & sum(r(:,end-1)) <1
      r(:,end-1) = [];
    extra_V.kk_approx(:,end-1) = [];
    extra_V.kk(:,end-1) = [];
    extra_V.ln_kk(:,end-1) = [];
  end

  [hp_posterior,extra_V] = mk_hp_posterior(data, r, hp_prior, opts,extra_V);
%   t1=sum(hp_posterior.a)
%   t2=sum(hp_posterior.b)
%   t3=hp_posterior.a./hp_posterior.b
% t4=sum(hp_posterior.a./hp_posterior.b)
end
% disp_status(free_energy, hp_posterior, opts);
disp(['### updating posterior ... done.'])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [r, I,extra_V] = sort_r(data, r, opts,extra_V);
disp('sorting...')
if opts.use_kd_tree
  Nc = [data.kdtree(:).N]*r; % 1*K
else
  Nc = sum(r, 1); % 1*K
end
[dummy,I] = sort(Nc(1:end-1), 2, 'descend');
I(end+1) = length(Nc);


r = r(:,I);
extra_V.kk = extra_V.kk(I);
extra_V.ln_kk = extra_V.ln_kk(I);
extra_V.kk_approx = extra_V.kk_approx(I);
disp('sorting... done.')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function disp_status(free_energy, hp_posterior, opts);
if isequal(opts.algorithm, 'vdp')
  Nc = hp_posterior.true_Nc;
else
  Nc = hp_posterior.Nc;
end
disp(['F=' num2str(free_energy) ...
      ';    Nc=[' num2str(Nc, ' %0.5g ') ...
      '];'])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bool = free_energy_improved(free_energy, new_free_energy, warn_when_increasing, opts);
diff = new_free_energy - free_energy;
if abs(diff/free_energy) < opts.threshold
  bool = 0;
elseif diff > 0
  if warn_when_increasing
    if abs(diff/free_energy) > 1.0e-3
      error(['the free energy increased.  the diff is ' num2str(new_free_energy-free_energy)])
    else
      warning(['the free energy increased.  the diff is ' num2str(new_free_energy-free_energy)])
    end
  end
  bool = 0;
elseif diff == 0
  bool = 0;
else
  bool = 1;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function partitions = init_kdtree_partitions(given_data, depth);
[D,N] = size(given_data);
root = mk_kdtree_partition(given_data, [1:N]);
partitions = expand_recursively(root, depth);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function partitions = expand_recursively(partition, depth);
children = mk_children_of_partition(partition);
if depth == 1 | length(children) == 1
  partitions = children;
else
  partitions1 = expand_recursively(children(1), depth-1);
  partitions2 = expand_recursively(children(2), depth-1);
  partitions = [partitions1 partitions2];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function children = mk_children_of_partition(partition);
if partition.N == 1
  children = partition;
  return
end
dir = divide_by_principal_component(partition.given_data(:,partition.indices), ...
                                    partition.sigma, ...
                                    partition.mean);
positive_I = find(dir>=0);
negative_I = find(dir<0);
if length(positive_I) == 0 | length(negative_I) == 0
  children = partition;
  return
end
children(1) = mk_kdtree_partition(partition.given_data, partition.indices(positive_I));
children(2) = mk_kdtree_partition(partition.given_data, partition.indices(negative_I));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%�������ɷ�
function direction = divide_by_principal_component(data, covariance, mean);
N = size(data, 2);
if size(data,1) <= 16
  [V,D] = eig(covariance);
  [eig_val, principal_component_i] = max(diag(D));
  principal_component = V(:,principal_component_i);
else
   %principal_component��ʾ��x������ͬ�ĵ�λ������eig_val��ʾx�Ķ�����
  [principal_component,eig_val] = power_method(covariance);
end
direction = sum((data - repmat(mean, 1, N)).*repmat(principal_component, 1, N), 1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function partition = mk_kdtree_partition(given_data, indices)
partition.N = length(indices);
if partition.N == 0
  error('indices must have at least one index.')
end
data_a = given_data(:,indices);
partition.given_data = given_data;
partition.indices = indices;
partition.sum_x = sum(data_a, 2); % D*1 �������
mean = partition.sum_x / partition.N; % D*1
partition.mean = mean; 
partition.sum_xx = data_a*data_a'; % D*D
%partition.sigma = partition.sum_xx / partition.N - mean*mean';
partition.sigma = corr(data_a');
partition.children = [];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function r = rand_r(data, K, opts);
% r: N*K
if opts.use_kd_tree
  N = length(data.kdtree);
else
  N = size(data.given_data, 2);
end
if isequal(opts.algorithm, 'vdp')
  r = zeros(N, K+1);    %4*2
else
  r = zeros(N, K);
end
r(:,1:K) = rand(N, K);
r = normalize(r, 2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [hp_posterior,extra_V]= mk_hp_posterior(data, r, hp_prior, opts,extra_V)
if opts.use_kd_tree
  N = length(data.kdtree);
  D = size(data.kdtree(1).sum_x, 1);
  Na = [data.kdtree(:).N];
  if isequal(opts.algorithm, 'vdp')
    true_Nc = Na*r; % 1*K
    r(:,end) = 0;
  end
  Nc = Na*r; % 1*K
  sum_x = [data.kdtree(:).sum_x] * r;
else
  [D,N] = size(data.given_data);
  if isequal(opts.algorithm, 'vdp')
    true_Nc = sum(r, 1); % 1*K
    r(:,end) = 0;
  end
  Nc = sum(r, 1); % 1*K
  sum_x = data.given_data * r;
end
threshold_for_N = 1.0e-200;
I = find(Nc>threshold_for_N);
inv_Nc = zeros(1,K);
inv_Nc(I) = 1./Nc(I);
Nc = Na*r; % 1*K

for k=1:K
    for T=1:N
        Z_tk=r(T,k);
        Nt = [data.kdtree(T).N];
        term1(:,T)=Nt * Z_tk;
        term2(:,T)=[data.kdtree(T).mean]*Nt * Z_tk;%D*1
    end
    hp_posterior.mu(:,k) = hp_prior.mu'*hp_prior.beta+sum(term2,2);%D*K
    hp_posterior.beta(k) = sqrt(sum(power(hp_posterior.mu(:,k),2)));
    hp_posterior.mu(:,k) = bsxfun(@rdivide,hp_posterior.mu(:,k), hp_posterior.beta(k));

    nu=D/2-1;
    bes_kk(k) = d_besseli(nu,extra_V.kk_approx(k));
    bes_bo(k) = d_besseli(nu,hp_prior.beta.*extra_V.kk_approx(k));
    bes_bk(k) = d_besseli_lower(nu,hp_posterior.beta(k).*extra_V.kk_approx(k));

    hp_posterior.a(k)= hp_prior.a + sum(term1)*nu + hp_posterior.beta(k).*extra_V.kk_approx(k).*bes_bk(k);
    hp_posterior.b(k) = hp_prior.b + sum(term1)*bes_kk(k)+hp_prior.beta.*bes_bo(k); 

    assert(min(hp_posterior.a)>0);
    assert(min(hp_posterior.b)>0);

    
    extra_V.kk_approx(k)=extra_V.kk(k);
        
    extra_V.kk_approx(k) = min(extra_V.kk_approx(k),opts.max_kk); % must choose a point where besseli is finite
    extra_V.kk(k) =hp_posterior.a(k)./hp_posterior.b(k);
    extra_V.ln_kk(k)= psi(hp_posterior.a(k))-log(hp_posterior.b(k));
end

hp_posterior.gamma = zeros(2,K);
hp_posterior.gamma(1,:) = 1 + true_Nc;
hp_posterior.gamma(2,:) = hp_prior.alpha + sum(true_Nc) - cumsum(true_Nc,2);

hp_posterior.Nc = Nc; 
hp_posterior.true_Nc = true_Nc;
hp_posterior.r = r; 
% t4=hp_posterior.beta.*extra_V.kk_approx.*bes_bk
% vmf�ķ���
t1 = reshape([data.kdtree(:).sum_xx],D,D,N);
sum_x = repmat(reshape([data.kdtree(:).sum_x],D,1,N),[1,D,1]);
Na_1 = repmat(reshape([data.kdtree(:).N],1,1,N),[D,D,1]);


for c=1:K
    t2 = sum_x.*repmat(hp_posterior.mu(:,c)',[D,1,N]);
    S1 = (t1 - 2*t2 )./Na_1+repmat(hp_posterior.mu(:,c)*hp_posterior.mu(:,c)',[1,1,N]);
% S1= [data.kdtree(:).sum_x]-hp_posterior.mu(:,c);
   S=S1/(Nc(:,c));
     hp_posterior.B{c} = sum(S,3);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fc = mk_E_log_q_p_eta(data, hp_posterior, hp_prior, opts,extra_V);
% returns E[ log q(eta)/p(eta) ]_q
% fc : 1 by K
D = size(hp_posterior.mu, 1);
K = size(hp_posterior.r, 2);

nu=D/2-1;
m_o=repmat(hp_prior.mu,K,1)';
E_m1 =nu*log(hp_posterior.beta./hp_prior.beta)-hp_posterior.beta+ hp_prior.beta+1000;
E_m2 =hp_posterior.beta.*extra_V.kk-hp_prior.beta.*extra_V.kk.*sum(m_o.*hp_posterior.mu);
t1=hp_posterior.beta.*extra_V.kk;    %fc����ǣ�ű���

E_m_H_m = E_m1+E_m2;
E_k_H_k1= hp_posterior.a.*log(hp_posterior.b)-gammaln(hp_posterior.a)+(hp_posterior.a-1).*extra_V.ln_kk-hp_posterior.b.*extra_V.kk;
E_k_H_k2 = hp_prior.a.*log(hp_prior.b)-gammaln(hp_prior.a)+(hp_prior.a-1).*extra_V.ln_kk-hp_prior.b.*extra_V.kk;
E_k_H_k = E_k_H_k1-E_k_H_k2;
fc = E_m_H_m+E_k_H_k;

debug
if length(find(E_m_H_m<-1.0e-8,1)) > 0
  E_m_H_m
  error('E_log_q_p_mean is negative.')
end
if length(find(E_k_H_k<-1.0e-8,1)) > 0
  E_k_H_k
  error('E_log_q_p_mean is negative.')
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [free_energy,S_tk] = mk_free_energy(data, hp_posterior, ...
                                                  hp_prior, opts,extra_V, ...
                                                  fc, S_tk);
if nargin == 5
  fc = mk_E_log_q_p_eta(data, hp_posterior, hp_prior, opts,extra_V); % 1*K
  S_tk = mk_S_tk(data, hp_posterior, hp_prior, opts,extra_V); % N*K
end

E_log_p_of_V = ...
  gammaln(sum(hp_posterior.gamma, 1)) ...
  - gammaln(1+hp_prior.alpha) ...
  - sum(gammaln(hp_posterior.gamma), 1) ...
  + gammaln(hp_prior.alpha) ...
  + ((hp_posterior.gamma(1,:)-1) ...
     .*(psi(hp_posterior.gamma(1,:))-psi(sum(hp_posterior.gamma,1)))) ...
  + ((hp_posterior.gamma(2,:)-hp_prior.alpha) ...
     .*(psi(hp_posterior.gamma(2,:))-psi(sum(hp_posterior.gamma,1))));

extra_term = sum(E_log_p_of_V);
t1=[data.kdtree(:).N]*log_sum_exp(S_tk, 2);
t2=sum(fc);
free_energy = extra_term +sum(fc)- [data.kdtree(:).N]*log_sum_exp(S_tk, 2);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function S_tk= mk_S_tk(data, hp_posterior, hp_prior, opts,extra_V);
% S_tk: N*K
% q(z_n=c|x_n) = lambda_n_c / sum_c lambda_n_c


N = length(data.kdtree);
K = size(hp_posterior.a, 2);
E_log_p_of_x= zeros(N,K);
S_tk = zeros(N,K);

for k=1:K
    for T=1:N
       E_log_p_of_x(T,k) = approximatebound([data.kdtree(T).mean],hp_posterior.mu(:,k),extra_V.kk(k),extra_V.ln_kk(k),extra_V.kk_approx(k),hp_posterior.Nc(k));
    end
    E_log_p_of_z_given_other_z_c = ...
        psi(hp_posterior.gamma(1,k)) ...
        - psi(sum(hp_posterior.gamma(:,k),1)) ...
        + sum(psi(hp_posterior.gamma(2,[1:k-1])) - psi(sum(hp_posterior.gamma(:,[1:k-1]),1)), 2);
    S_tk(:,k) = E_log_p_of_x(:,k) + E_log_p_of_z_given_other_z_c;
end
S_tk(:,end) = S_tk(:,end) - log(1- exp(psi(hp_prior.alpha) - psi(1+hp_prior.alpha)));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [r, data, S_tk] = mk_r(data, hp_posterior, hp_prior, opts, extra_V,S_tk);
% r: N*K
if nargin == 5
  S_tk = mk_S_tk(data, hp_posterior, hp_prior, opts,extra_V);
end
r = exp(normalizeln(S_tk, 1));
if opts.weight ~= 1
  r = r .* repmat(opts.weight, 1, size(r,2));
else
  r = r;
end

function M =normalizeln(M ,dimension);
M = lpt2lcpt(M, dimension);

function lcpt=lpt2lcpt(lpt,dimension);

the_other_dimension=-(dimension-1.5)+1.5;
lpt=permute(lpt,[dimension,the_other_dimension]);
% now we can calculate as if dimension=1.

log_sum_exp_lpt = log_sum_exp(lpt,2); % Mx1
lcpt = lpt - repmat(log_sum_exp_lpt,1,size(lpt,2));

lcpt=permute(lcpt,[dimension,the_other_dimension]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function hp_prior = mk_hp_prior(data, opts)

if ~isfield(opts, 'mu')
    hp_prior.mu=sum(data.given_data')/norm(sum(data.given_data'));
    %prior.mu=sum(data)/norm(sum(data));
else 
  hp_prior.mu= 0.01;
end

if isfield(opts, 'alpha')
  hp_prior.alpha = opts.alpha;
else
  hp_prior.alpha = 1;
end

hp_prior.beta=1; 
hp_prior.a=1; % assume unconcentrated data
hp_prior.b=0.01; % but do not trust that too much







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function m=normalize(m,dim)
dims = ones(1, ndims(m));
dims(dim) = size(m, dim);
m = m ./ repmat(sum(m, dim), dims);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Y] = detln( X )
% Y = logdet( X )
% return a log determinant of X
[d err] = chol(X);
if err
  error('error in Choleski disposition for detln');
end
Y = sum(log(diag(d))) *2;

function MIhat = nmi(A, B)
%NMI Normalized mutual information
% A, B: 1*N;
if length(A) ~= length(B)
    error('length( A ) must == length( B)');
end
N = length(A);
A_id = unique(A);
K_A = length(A_id);
B_id = unique(B);
K_B = length(B_id);
% Mutual information
A_occur = double (repmat( A, K_A, 1) == repmat( A_id', 1, N ));
B_occur = double (repmat( B, K_B, 1) == repmat( B_id', 1, N ));
AB_occur = A_occur * B_occur';
P_A= sum(A_occur') / N;
P_B = sum(B_occur') / N;
P_AB = AB_occur / N;
MImatrix = P_AB .* log(P_AB ./(P_A' * P_B)+eps);
MI = sum(MImatrix(:));
% Entropies
H_A = -sum(P_A .* log(P_A + eps),2);
H_B= -sum(P_B .* log(P_B + eps),2);
%Normalized Mutual information
MIhat = MI / sqrt(H_A*H_B);

function [FMeasure,Accuracy] = Fmeasure(P,C)
% PΪ�˹���Ǵ�
% CΪ�����㷨������
N = length(C);% ��������
p = unique(P);
c = unique(C);
P_size = length(p);% �˹���ǵĴصĸ���
C_size = length(c);% �㷨����Ĵصĸ���
% Pid,Rid���������ݣ���i�з������ݴ�����������ڵ�i����
Pid = double(ones(P_size,1)*P == p'*ones(1,N) );
Cid = double(ones(C_size,1)*C == c'*ones(1,N) );
CP = Cid*Pid';%P��C�Ľ���,C*P
Pj = sum(CP,1);% ��������P��C�������еĸ���
Ci = sum(CP,2);% ��������C��P�������еĸ���
 
precision = CP./( Ci*ones(1,P_size) );
recall = CP./( ones(C_size,1)*Pj );
F = 2*precision.*recall./(precision+recall);
% �õ�һ���ܵ�Fֵ
FMeasure = sum( (Pj./sum(Pj)).*max(F) );
Accuracy = sum(max(CP,[],2))/N;


% Local Variables: ***
% mode: matlab ***
% End: ***
