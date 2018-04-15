function [y, w, wt] = gradientDescent(s, R, p, mu, e)
% GRADIENTDESCENT - Perform gradient descent to compute Wiener coefficients
%   
% SYNTAX
%
%   [Y, W, WT] = GRADIENTDESCENT( S, P, R, MU )
%   [Y, W, WT] = GRADIENTDESCENT( S, P, R, MU, E )
%
% INPUT
%
%   S   Input singal to adaptive filter         [n-vector]
%   R   Auto-correlation matrix                 [m-by-m]
%   P   Cross-correlation vector                [m-vector]
%   MU  Step value                              [scalar]
%   E   Error threshold (optional)              [scalar]
%
% OUTPUT
%
%   Y   Adaptive filter output signal           [n-vector]
%   W   Final Wiener coefficients               [m-vector]
%   WT  Evolution/adaptation of Wiener coeff    [m-by-n]
%
% DESCRIPTION
%
%   GRADIENTDESCENT 
%
% DEPENDENCIES
%
%   <none>
%
%
% See also      
%
  if nargin == 4
      n = length(s);
      m = length(p);

      w = rand(m, 1); 

      wt = zeros([m,n]); wt(:,1) = w;
      y = zeros(n, 1);

      for i=m:n
        w = w + mu*(p-R*w); % Adaptation steps
        wt(:,i) = w;
        y(i) = s(i:-1:i-m+1)' * w; % filter
      end
      return;
  elseif nargin == 5 && e > 0
      n = length(s);
      m = length(p);

      w = rand(m, 1); 

      wt = zeros([m,n]); wt(:,1) = w;
      y = zeros(n, 1);

      for i=m:n
        w = w + mu*(p-R*w); % Adaptation steps
        if norm(p-R*w) < e
            for j=i:n
                wt(:,j) = w;
                y(j) = s(j:-1:j-m+1)' * w; % filter
            end
            break;
        end
        wt(:,i) = w;
        y(i) = s(i:-1:i-m+1)' * w; % filter
      end
      return;
  else
      disp('Unsupported input parameters');
  end
  
end


%%------------------------------------------------------------
%
% AUTHORS
%
%   Dimitris Floros                         fcdimitr@auth.gr
%   Nikolaos Katomeris  8551                ngkatomer@auth.gr
%
% VERSION
%
%   0.2 - April 13, 2018
%
% CHANGELOG
%   0.2 (Apr 13, 2018) - Nikolaos
%       * added optional error threshold parameter
%
%   0.1 (Mar 12, 2018) - Dimitris
%       * initial implementation
%
% ------------------------------------------------------------

