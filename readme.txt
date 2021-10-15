
Develop a soft K-means algorithm where each of the Gaussians has a
diagonal covariance matrix given by,

E = [sigmak^2     0
        0    sigmak^r^2]

and the four variance parameters sigma1^2, sigma2^2, sigma1^r^2 
and sigma2^r^2, as well as the two mixing coeficients pi1 and pi2 
are now updated at each iteration using the maximum likelihood 
results derived in lectures for the univariate Gaussian. In this case,
simply replace the scalar values xn and k in each of the formulas, 
with the vector values xn and muk. Note that the scalar variance ^2
k then becomes a vector: ^sigmak^2 = [sigmak^2, sigmak^r^2]