3c157      N=100：500：700   k=10：20：30     mu:001;010;100
3c149      N=100：400：900  	
3c258     N=200：500：800   k=10：5：15
4c2583     N=200：500：800 300  k=10：20：30: 40   mu:001;010;100; [0,8,-20]
4c2583dak     N=200：500：800 300  k=100：60：150: 40   mu:001;010;100; [0,8,-20]
2c30500      N=30：500 	k=100：60 mu0,0,1;0,1,0
4c2583lable   N=200：500：800 300   k=70：80：100: 200   mu:001;010;100; [0,8,-20]
3c258lable   N=200：500：800   k=70：80：100   mu:001;010;100
3c123label N=100：200：300     k=25 25 25   mu:[0,-10,1] [-10,1,0] [0,8,-20];
3c433label N=400：300：300    k=15 15 15   mu:[0,-20,8] [-20,1,4] [1,-9,-20];
3c433label_lessk  N=400：300：300    k=10 10 10   mu:[0,-20,8] [-20,1,4] [1,-9,-20];
4c4563label  N=400：500：600 300   k=20 22 25 30 mu:001;010;100; [0,8,-20]

Variational Dirichlet Process Gaussian Mixture Models
written by Kenichi Kurihara
distributed under the modified BSD license.



ALGORITHMS (see Refrences, the end of this document)

1. accelerated variational Dirichlet process Gaussian mixture model
2. collapsed variational stick-breaking Dirichlet process Gaussian mixture model
3. variational Gaussian mixture model with a collapsed Dirichlet prior.
4. variational Dirichlet process Gaussian mixture model by Blei and Jordan.


USAGE

>> result = vdpgm(X, opts);

The first argument is data. Each data point is a column vector of X.

The second argument opts is the option of this program which
determines an algorithm and hyperparameters. You can set opts as you
want, or basic option generators are also available.

>> opts = mkopts_avdp;     % for the algorithm 1
>> opts = mkopts_csb(10);  % for the algorithm 2 with T=10 
>> opts = mkopts_cdp(10);  % for the algorithm 3 with K=10 
>> opts = mkopts_bj(10);   % for the algorithm 4 with T=10 
  
Although opts accepts many options, some options are exclusive.

The output result is a structure containing parameters for posteriors.
Maybe, the most useful result is result.q_of_z which is the posterior
probability of assignments. q_of_z is a N by K (or T) matrix
s.t. sum_c q_of_z(i,c) = 1 for any c. q_of_z is available only when
opts.get_q_of_z is set to 1.

Other useful stats:
分量c的协方差的期望值
- The expected value of the covariance of component c,
>> result.hp_posterior.B{c} / results.hp_posterior.eta(c)
分量c的质心的期望值
- The expected value of the centroid of component c,
>> result.hp_posterior.m(:,c)

One may want to know the number of discovered clusters.
If opts.algorithms is 'vdp', it is 
>> results.K - 1

Otherwise, results.K is initialized by opts, and does not change.  In
these cases, the number of clusters is K s.t.
result.hp_posterior.Nc(K+1) is close enough to zero.



REFERENCES

    * Kenichi Kurihara, Max Welling and Yee Whye Teh,
      Collapsed Variational Dirichlet Process Mixture Models,
      the Twentieth International Joint Conference on Artificial Intelligence (IJCAI 2007). 

    * Kenichi Kurihara, Max Welling and Nikos Vlassis,
      Accelerated Variational Dirichlet Mixture Models,
      Advances in Neural Information Processing Systems 19 (NIPS 2006). 

    * David M. Blei and Michael I. Jordan,
      Variational Inference for Dirichlet Process Mixtures,
      Bayesian Analysis, Vol.1, No.1, 2005.
