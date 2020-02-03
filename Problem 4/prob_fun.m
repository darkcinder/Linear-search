function p = prob_fun(label, x, theta)
% Binary probability distribution, must learn parameter theta to 
% maximize the probability P(label = L | x, theta).
    p = 1/(1 + exp(-label*theta'*x));
end
