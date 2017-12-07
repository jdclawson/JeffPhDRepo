function  [yt,yd] = hp(y,lambda)
  % Sparse implentation
  % A larger lambda results in a smoother trend series. 
  % For quarterly data Hodrick and Prescott(1980) use phi=1600.
  % Also returns the series with no trend: yd(:,i)=y(:,i)-yt(:,i)
  %                                       i=1,...k
  %
  % References:
  %    Measurement'' QUARTERLY REVIEW, Federal Reserve Bank 
  %    of Minneapolis, Fall 1986.
  %
  % Jesus Fernandez-Villaverde, 4-26-01
  T = length(y);
	K = spdiags([ones(T-2,1) -2*ones(T-2,1) ones(T-2,1)],0:2,T-2,T);
  A = sparse(eye(T)+lambda*(K'*K));
 	yt = A\y;
  yd = y-yt;