function rtnorm(a, b, mu, sigma)
%%
%   Function for samples truncated normal random variables:
%   f(x, a, b, mu, sigma) \propto I(x \in [a,b]) e^{-\frac{(x-\mu)^2}{2\sigma^2}} 
%
%   a     - (nx1) left limit
%   b     - (nx1) right limit
%   mu    - (nx1) means
%   sigma - (nx1) standrad deviations
%   wrapper funciton till rtnorm (mex) function
%
%%