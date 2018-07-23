% ---------------------
% simple script to test
% qrupdate_rank_1
% ---------------------

% ---------------------------------------
% should also work for rectangular n ~= m
% ---------------------------------------
n = 10;
m = 10;
A = rand(n,m);

use_complex = 1;
if (use_complex),
  A = A + rand(n,m)*sqrt(-1);
end;

[Qin,Rin] = qr(A);
u = rand(n,1);
v = rand(m,1);

[Q,R] = qrupdate_rank_1(Qin,Rin, u, v );

nrows_Q = size(Q,1);
ncols_Q = size(Q,2);


nerrors = 0;
tol = 10^(-7);

% -----------------------------------
% double check, R is upper triangular
% -----------------------------------
err_R = norm( triu(R) - R, 1);
is_ok_R = (err_R < tol );
if (~is_ok_R),
  disp(sprintf('qrupdate_rank_1: R is not upper triangular, err_R=%g', ...
	  err_R));
  nerrors = nerrors + 1;
end;

% ----------------------------
% double check Q is orthogonal
% ----------------------------
err_Q = norm(Q'*Q - eye(ncols_Q,ncols_Q),1);
is_ok_Q = (err_Q < tol );
if (~is_ok_Q),
  disp(sprintf('qrupdate_rank_1: Q is not orthogonal, err_Q=%g', ...
	   err_Q));
  nerrors = nerrors + 1;
end;

A = Qin * Rin;
err_A = norm( (A + u*v') - Q*R,1);
is_ok_A = (err_A < tol );
if (~is_ok_A),
  disp(sprintf('qrupdate_rank_1: invalid QR, err_A=%g', ...
	     err_A));
  nerrors = nerrors + 1;
end;


if (nerrors == 0),
  disp(sprintf('qrupdate_rank_1 passed simple test, n=%d, m=%d', ...
               n, m ));
end;
