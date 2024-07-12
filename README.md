# Perturbed-factor-model
1) LRFpertgrp1.R and LRFpertgrp2.R run the codes for the two examples from the paper (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9624461/) simulation case 2.
2) FBPFA-PFA.R has the function PFA() that can fit both Full-Bayes version or you can specify alpha too.
3) There are two sets of functions FBPFA-PFA with fixed latent dim.R and FBPFA-PFA.R The first one will need a pre-specified latent variable dimension, but cumulative shrinkage prior will still be applied. The second one will set the dimension based on either PCA (if ini.PCA=T) or set it as 'p', the dimension of the data.
4) applicationOF_FBPFA.R provides an example usage of above fn.
5) Under the heteroscedastic latent factor model of PFA, the covariance matrix is LSigma_1L^T + Sigma_e, where Sigma_1 is the variance of the latent variables. Then fit$Loading[i]%%diag(fit$Latentsigma[i]^2)%%t(fit$Loading[i]) + diag(fit$Errorsigma[i]^2) is the i-th posterior samples of Covariance.
Similarly fit$Loading[i]%%diag(fit$Latentsigma[i]) are the posterior samples of the loading matrix which is comparable to other papers. So, one can apply the Varimax or the one discussed in the paper to get a point estimate of the loading matrix. 
(Reduce('+', fit$Loading[200:499])/length(fit$Latentsigma[200:499])) %*% diag(sigma2p) also gives a point estimate of the loading matrix without any adjustment. For a simple simulation setting like mine, it worked fine.
6) The plotload.R file contains a function, built on IMIFA package is useful to make nice plots for loading matrices.
