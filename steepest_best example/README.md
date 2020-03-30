P(L|x,theta)=1/(1+exp(-L*theta*X))

min(theta)0.5*||theta||^2-sum(logP(L|X,theta))
