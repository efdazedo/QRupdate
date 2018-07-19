% --------------------------------------
% simple script to test qr_symm_remove()
% --------------------------------------
n = 5;
m = n;
A = rand(n,m);
[Qin,Rin] = qr(A);
k = 3;

[Q,R] = qr_symm_remove(k, Qin, Rin);

Anew = A;
Anew(:,k) = 0;
Anew(k,:) = 0;
Anew(k,k) = 1;

tol = 10^(-7);
nerrors  = 0;

err_R = norm(tril( R,-1),1);
is_ok_R = (err_R < tol);
if (~is_ok_R),
  disp(sprintf('R is not upper triangular, err_R=%g', err_R));
  nerror = nerror + 1;
end;

ncols_Q = size(Q,2);
err_Q = norm( Q'*Q - eye(ncols_Q,ncols_Q),1);
is_ok_Q = (err_Q < tol);
if (~is_ok_Q),
   disp(sprintf('Q is not orthogonal, err_Q=%g', err_Q));
   nerror = nerror + 1;
end;


err_Anew = norm(Q*R - Anew,1);
is_ok_Anew = (err_Anew < tol);
if (~is_ok_Anew),
  disp(sprintf('error in QR of Anew, err_Anew=%g', err_Anew));
  nerrors = nerrors + 1;
end;



if (nerrors == 0),
  disp(sprintf('qr_symm_remove passed simple test, n=%d',n));
end;
