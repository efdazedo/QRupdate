% simple test of qrupdate_rc
%
n = 10;
A = rand(n,n);

% -----------------------------------
% create QR factorization by adding one row and one column 
% using qrupdate_rc()
% -----------------------------------

% --------------------------------------------------
% start with trivial  QR factorization of 1x1 matrix
% --------------------------------------------------
Q = eye(1,1);
R = A(1,1);
for k=1:(n-1),
   kp1  = k + 1;
   Qin = eye( kp1, kp1);
   Qin(1:k,1:k) = Q(1:k,1:k);
   Qin(kp1,kp1) = 1;

   Rin = zeros( kp1, kp1);
   Rin(1:k,1:k) = R(1:k,1:k);
 
   % ---------------------------
   % augment extra row and colum
   % ---------------------------
   Rin(1:k, kp1) = A(1:k,kp1);
   Rin(kp1, 1:k) = A(kp1,1:k);
   Rin(kp1, kp1) = A(kp1,kp1);
   
   [Q,R] = qrupdate_rc(k, Qin, Rin); 
end;


nerrors = 0;
% -----------------------------
% check  Q should be orthogonal
% -----------------------------
tol = 1e-7;
err_Q = norm(Q'*Q-eye(n,n),1);
is_ok_Q = (err_Q < tol );
if (~is_ok_Q),
   disp(sprintf('Q is not orthogonal, err_Q=%g', err_Q));
   nerrors = nerrors+1;
end;

% ---------------------------
% check R is upper triangular
% ---------------------------
err_R = norm( R - triu(R),1);
is_ok_R = (err_R < tol);
if (~is_ok_R),
  disp(sprintf('R is not upper triangular, err_R=%g',err_R));
  nerrors = nerrors + 1;
end;



% -------------
% check A = Q*R
% -------------
err_A = norm(A - Q*R,1);
is_ok_A = (err_A < tol );
if (~is_ok_A),
   disp(sprintf('A not equal Q*R, err_A=%g', err_A));
   nerrors  = nerrors + 1;
end;
   
   
if (nerrors == 0),
  disp(sprintf('n=%d, qrupdate_rc passed simple test',n));
end;

