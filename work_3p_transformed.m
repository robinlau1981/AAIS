function [logprior,loglike]=work_3p_transformed(thetas,data)
loglike=loglikelihood_transformed_17d(thetas,data);
logprior=log_prior_transformed_v5(thetas);