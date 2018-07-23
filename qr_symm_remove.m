function [Q,R] =  qr_symm_remove( k, Qin, Rin, Arow_k, Acol_k )
% [Q,R] =  qr_symm_remove( k, Qin, Rin, [Arow_k, Acol_k ])
%
%  Essentially remove the k-th row and k-th column by
%  zeroing out row    A(:,k) = 0, 
%  zeroing out column A(k,:) = 0 
%  placing 1 on diagonal A(k,k) = 1
%
%  ---------------------------------------------
%  update QR factorization as two rank-1 updates
%  ---------------------------------------------
%

% -----------------------------------------
% ek is the k-th column of identity matrix
% -----------------------------------------
n = size(Qin,1);
ek = zeros(n,1); ek(k) = 1;

if (nargin == 3),
  % -----------------------------------------------
  % either pass in k-th row and k-th column
  % or
  % need to regenerate the k-th row and k-th column
  % from Q and R in O(n^2) work
  % -----------------------------------------------

  use_multiply = 0;
  if (use_multiply),
    Arow_k = (ek'*Qin)*Rin;
    Acol_k = Qin*(Rin*ek);
  else
    Arow_k = Qin(k,:) * Rin;
    Acol_k = Qin * Rin(:,k);
  end;
  
end;
Akk = Acol_k(k);


% --------------------------------------------
% perform two rank-1 updates
% 
% A <-   A + (-ek) * A(k,:)' + (-A(:,k)+(A(k,k)+1)*ek) * ek' ;
% --------------------------------------------

% -------------------
% first rank-1 update
% -------------------
u = -ek; 
v = Arow_k(:)';
v = reshape(v, numel(v),1);

% -------------
% Note save some work by
% w_in = Qin'*u 
%      = Qin'*(-ek) 
%      = -Qin'*ek 
%      = -Qin(k,:)';
% -------------
w_in = -Qin(k,:)';

[Qin2,Rin2] = qrupdate_rank_1( Qin, Rin, u, v, w_in );

% --------------------
% second rank-1 update
% --------------------
u2 = -Acol_k ;

% --------------------
% u2 = u2 + (Akk+1)*ek;
% --------------------
u2(k) = 1;

v2 = ek;
[Q,R] = qrupdate_rank_1( Qin2, Rin2, u2, v2 );

end


