# MarginalLikelihoodEstimator
Matlab code for the adaptive annealed importance sampling based marginal likelihood estimator

This is a version of the adaptive annealed importance sampler (AAIS) [3] designed for user-friendliness. Below is a brief introduction to the algorithm's background, followed by instructions on how to use the code.

%% -------------------  Background ------------------- %%
Code developer：Bin Liu;   Algorithm contributor: Bin Liu, Jim Berger (Duke University), Tom Loredo (Cornell University), Merlise Clyde (Duke University)

The initial version of the algorithm was developed in 2009 by Bin during his postdoctoral work at the Department of Statistical Science, Duke University, under the support of NSF Grant AST-0507481. Its development was motivated by a challenging computational problem, namely marginal likelihood estimation, in the field of exoplanet detection, detailed in [1]. Tom introduced Bin to the astrophysical background; Jim and Merlise introduced previous efforts related to this problem. Jim provided valuable insights through many collaborative discussions. He also contributed a crucial idea to the algorithm: adding new mixture components at the sample with the biggest importance weight. Bin devised techniques for combining the annealing strategy with adaptive importance sampling and mixture modeling, with the added feature of adjusting the number of mixing components adaptively. In summary, this algorithm is the result of close collaboration between Bin, Jim, Tom, and Merlise.

Bin initially presented the initial version of the algorithm in an invited talk titled 'Adaptive T-mixture Importance Sampling Method' at the Sequential Monte Carlo (SMC) Method Transition Workshop, hosted by the Statistical and Applied Mathematical Sciences Institute (SAMSI) on November 9-10, 2009. Tom later provided a brief introduction to the algorithm, then called annealed adaptive importance sampling, in one of his conference presentations in 2010 [2]. Preliminary simulation results for this algorithm were reported in [1]. Later, Bin made a lot of efforts to improve the algorithm by testing and updating it on simulated and real datasets. A comprehensive introduction to the final version of this algorithm, along with its applications to model selection and extrasolar planet detection, was published in [3].

Reference:
[1] Loredo, T.J., Berger, J.O., Chernoff, D.F., Clyde, M.A. and Liu, B., Bayesian Methods for Analysis and Adaptive Scheduling of Exoplanet Observations, Statistical Methodology, vol.9, no.1-2, pp.101-114, 2012. 
[2] Loredo, T.J., Bayesian methods for exoplanet science: Planet detection, orbit estimation, and adaptive observing, ADA 6-Sixth Conference on Astronomical Data Analysis, May 6, 2010. https://www.researchgate.net/publication/230631076_Bayesian_methods_for_exoplanet_science_Planet_detection_orbit_estimation_and_adaptive_observing
[3]  Liu, B., Adaptive Annealed Importance Sampling for Multimodal Posterior Exploration and Model Selection with Application to Extrasolar Planet Detection, The Astrophysical Journal Supplement Series, vol. 213, no.14, pp.1-16, 2014. doi:10.1088/0067-0049/213/1/14.

%% -------------------  How to use this code ------------------- %%
The code here shows an example application of our algorithm AAIS for multimodal posterior sampling and marginal likelihood estimation.
The likelihood function is set as an outer product of 7 univariate distributions as below
(1) 3/5*G(10+x|2,3)+2/5*G(10-x|2,5);
(2) 3/4*skN(x|3,1,5)+1/4*skN(x|-3,3,-6);
(3) S(x|0,9,4);
(4) 1/2*B(x+3|3,3)+1/2*N(x|0,1);
(5) 1/2*e(x|1)+1/2*e(-x|1);
(6) skN(x|0,8,-3);
(7) 1/8*N(x|-10,0.1)+1/4*N(x|0,0.15)+5/8*N(x|7,0.2),
where G denotes the gamma distribution, skN the skew-normal distribution, B the beta distribution, and e the exponential distribution. The 2nd dimension has two modes bracketing a deep ravine, the 4th dimension has one low, broad mode that overlaps a second sharper mode, and the 7th dimension has 3 distinct, well-separated modes. Only the 5th dimension is symmetric. There is a range of tail behaviors as well, from Gaussian to heavy-tailed. This likelihood function was first used in [4], then used in [2][3].

In the code here, AAIS is used to estimate the marliginal likelihood (i.e., model evidence), which is an integral of the above likelihood function over the parameter space.  Just run 'main.m' to see the results.

You can replace the likelihood function used here with your own one. Then run the algorithm to get the corresponding estimation results (You likely need to re-initialize the algorithm according to your problem setting, e.g., by specifying the dimension of your likelihood function, the sample size, the initial proposal function, your anneling schedule). 

[4] Crooks, J. L., Berger, J. O.,and Loredo, T. J., Posterior-Guided Importance Sampling for Calculating Marginal Likelihoods with Application to Bayesian Exoplanet Searches, Discussion Paper Series of Dept. of Statitical Science, 2007.
