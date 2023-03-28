# MarginalLikelihoodEstimator
This is matlab code for the adaptive annealed importance sampling (AAIS) based marginal likelihood estimator, reported in [3] https://iopscience.iop.org/article/10.1088/0067-0049/213/1/14 .

This is a version of the code designed for user-friendliness. Another more complex version using the real astro model and datasets can be found at: https://github.com/robinlau1981/AAIS .

Below is a brief introduction to the algorithm's background, followed by an instruction on how to use the code.

%% -------------------  Background ------------------- %%

Code developer：Bin Liu;   

Algorithm contributors: Bin Liu, Jim Berger (Duke University), Tom Loredo (Cornell University), Merlise Clyde (Duke University)

The initial version of this algorithm was developed in 2009 by Bin during his postdoctoral work at the Department of Statistical Science, Duke University, under the support of NSF Grant AST-0507481. Its development was motivated by a challenging computational problem, namely marginal likelihood estimation, in the field of exoplanet detection, detailed in [1]. Tom introduced Bin to the astrophysical background; Jim, Tom, and Merlise introduced previous efforts that have been made related to this problem, especially [4]. Jim provided valuable insights through many collaborative discussions with Bin. He contributed a crucial idea to the algorithm: adding new mixture components at the sample with the biggest importance weight. Bin devised techniques for combining the annealing strategy with adaptive importance sampling and mixture modeling, with the added feature of adjusting the number of mixing components adaptively. In summary, the initial version of this algorithm is the result of close collaboration between Bin, Jim, Tom, and Merlise.

Bin initially presented the algorithm in an invited talk titled 'Adaptive T-mixture Importance Sampling Method' at the Transition Workshop (on November 9-10, 2009), of a year-long project called Sequential Monte Carlo (SMC) Methodology hosted by the Statistical and Applied Mathematical Sciences Institute (SAMSI). Tom later provided a brief introduction to the algorithm, then called annealed adaptive importance sampling, in one of his conference presentations in 2010 [2]. Preliminary simulation results for this algorithm were reported in [1]. Later, Bin made efforts in testing the algorithm on more simulated and especially real astro datasets and improving it by solving some technical issues. The final version of the algorithm, along with details on how to use it for model selection and exoplanet detection, was published in [3].

Reference:

[1] Loredo, T.J., Berger, J.O., Chernoff, D.F., Clyde, M.A. and Liu, B., Bayesian Methods for Analysis and Adaptive Scheduling of Exoplanet Observations, Statistical Methodology, vol.9, no.1-2, pp.101-114, 2012. 

[2] Loredo, T.J., Bayesian methods for exoplanet science: Planet detection, orbit estimation, and adaptive observing, ADA 6-Sixth Conference on Astronomical Data Analysis, May 6, 2010. https://www.researchgate.net/publication/230631076_Bayesian_methods_for_exoplanet_science_Planet_detection_orbit_estimation_and_adaptive_observing

[3]  Liu, B., Adaptive Annealed Importance Sampling for Multimodal Posterior Exploration and Model Selection with Application to Extrasolar Planet Detection, The Astrophysical Journal Supplement Series, vol. 213, no.14, pp.1-16, 2014. doi:10.1088/0067-0049/213/1/14.

[4] Crooks, J. L., Berger, J. O.,and Loredo, T. J., Posterior-Guided Importance Sampling for Calculating Marginal Likelihoods with Application to Bayesian Exoplanet Searches, Discussion Paper Series of Dept. of Statitical Science, Duke University, 2007.

%% -------------------  How to use this code ------------------- %%

The code released here shows an example application of the AAIS algorithm for multimodal posterior sampling and marginal likelihood estimation. The marliginal likelihood (i.e., model evidence) is defined to be an integral of a target likelihood function over its parameter space. The target likelihood function involved here is an outer product of 7 univariate distributions as below

(1) 3/5 * G(10+x|2,3)+2/5 * G(10-x|2,5);

(2) 3/4 * skN(x|3,1,5)+1/4 * skN(x|-3,3,-6);

(3) S(x|0,9,4);

(4) 1/2 * B(x+3|3,3)+1/2 * N(x|0,1);

(5) 1/2 * e(x|1)+1/2 * e(-x|1);

(6) skN(x|0,8,-3);

(7) 1/8 * N(x|-10,0.1)+1/4 * N(x|0,0.15)+5/8 * N(x|7,0.2),

where G denotes the gamma distribution, skN the skew-normal distribution, B the beta distribution, and e the exponential distribution. This function was first designed in [4], then adopted in [2][3].

Just run 'main.m' to see the estimation results. You can replace the likelihood function used here with your own one. Then run 'main.m' to get the corresponding estimation results (you likely need to re-initialize the algorithm according to your problem setting, e.g., by specifying the dimension of your likelihood function, the sample size, the initial proposal function, your annealing schedule). 

If you find this work useful, please kindly cite following papers:

@article{liu2014adaptive,

  title={Adaptive annealed importance sampling for multimodal posterior exploration and model selection with application to extrasolar planet detection},
  
  author={Liu, Bin},
  
  journal={The Astrophysical Journal Supplement Series},
  
  volume={213},
  
  number={1},
  
  pages={14},
  
  year={2014},
  
  publisher={IOP Publishing}
  
}

@article{loredo2012bayesian,

  title={Bayesian methods for analysis and adaptive scheduling of exoplanet observations},
  
  author={Loredo, Thomas J and Berger, James O and Chernoff, David F and Clyde, Merlise A and Liu, Bin},
  
  journal={Statistical Methodology},
  
  volume={9},
  
  number={1-2},
  
  pages={101--114},
  
  year={2012},
  
  publisher={Elsevier}
  
}
