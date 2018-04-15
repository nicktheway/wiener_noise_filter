function [R, p] = computePR(u, d, m)
% COMPUTEPR - Compute auto-correlation and cross-correlation
%   
% SYNTAX
%
%   [R, P] = COMPUTEPR( U, D, M )
%
% INPUT
%
%   U   Input signal                    [n-vector]
%   D   Desired signal                  [n-vector]
%   M   Filter length (tap)             [scalar]
%
% OUTPUT
%
%   R   Auto-correlation matrix         [M-by-M]
%   P   Cross-correlation vector       [M-vector]
%
% DESCRIPTION
%
%   COMPUTEPR returns the autocorrelation matrix and
%   crosscorelation vector, respecively. It uses xcorr.
%
% DEPENDENCIES
%
%   <none>
%
%
% See also      autocorr, xcorr
%
  
  % Make sure the inputs are column vectors
  u = u(:);
  d = d(:);
  
  % compute the auto-correlation vector with the needed lags
  r = xcorr(u, m-1, 'biased');
  r = r(m:end);
  
  % generate the toeplitz R matrix from the above vector
  R = toeplitz(r);
  
  % compute the cross-correlation vector
  pp = xcorr(d, u, m-1, 'biased');
  p = pp(m:end); 
  
end


%%------------------------------------------------------------
%
% AUTHOR
%
%   Nikolaos Katomeris, 8551, ngkatomer@auth.gr
%
% VERSION
%
%   1.0 - April 13, 2018
%
% CHANGELOG
% 
%   1.0 (April 13, 2018) - Nikolaos Katomeris
%       * initial implementation
%
% ------------------------------------------------------------

