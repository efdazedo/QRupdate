function [Q,R] = qrupdate_rank_1(Qin,Rin, u,v)
% [Q,R] = qrupdate_rank_1(Qin,Rin, u,v)
%
% ---------------------------------------
% update QR factorization
% A = Qin * Rin
% Q * R = (A + u * v')
%
% algorithm from section 12.5
% in book "Matrix Computations" by Golub and Van Lan
% ---------------------------------------
idebug = 0;

% ------------------------------------------
% note reuse storage of Qin, Rin  in C code
% ------------------------------------------
Q = Qin;
R = Rin;

nrows_Q = size(Q,1);
ncols_Q = size(Q,2);

nrows_R = size(R,1);
ncols_R = size(R,2);

n = numel(u);
m = numel(v);
u = reshape( u, n, 1 );
v = reshape( v, m, 1 );

isok = (n == ncols_Q) && (m == ncols_R);
if (~isok),
  error(sprintf('qrupdate_rank_1: invalid sizes, m=%d,n=%d',m,n));
  return;
end;


% -------------------------------------------
% Let A + u*v' = Q*(R + w*v'), where w = Q'*w
% -------------------------------------------
w = Q'*u;

% ----------------------------------------------------------- 
% Our goal is to apply Given's rotations to make w = norm(w) * e1
% e1 is 1st column of identity matrix
% upper triangular matrix R will be modified but 
% will be upper Hessenberg
%
% the Given's rotation are also applied to Q
%
% then finally note that
% Rnew + norm(w) * e1 * v'
% modifies ONLY the 1st row of Rnew and after
% the update, Rnew  is still upper Hessenberg
%
% then perform Given's rotation to make Rnew from
% upper Hessenberg to upper triangular
% --------------------------------------------------------- 

nw = numel(w);
for k=(nw-1):-1:1, 
  sa = w(k);
  sb = w(k+1);
  if (sb  == 0),
     continue;
  end;


  [c,s] = srotg(sa,sb);
  G = [c, -s;  ...
       s, c];

  % ----------
  % apply to w
  % ----------
  w(k) = c*sa + (-s)*sb;
  w(k+1) = 0;

  % ------------------------
  % apply to R(:(k+1), k:m)
  % ------------------------
  R(k:(k+1),k:m) = G * R(k:(k+1),k:m);

  % ---------------------------
  % apply G' to Q(1:n, k:(k+1)) 
  % ---------------------------
  Q(1:n, k:(k+1)) = Q(1:n,k:(k+1)) * G';
end;

% ------------------
% rank-1 update to R
% ------------------
R(1,1:m) = R(1,1:m) + w(1) * v(1:m)';

% -------------------------------------------------------
% now R is upper Hessenberg form, with extra sub-diagonal 
% under diagonal, need Given's rotation to
% reduce back to upper triangular
% -------------------------------------------------------
  
nk = min( nrows_R, ncols_R);
for k=1:ncols_R,
  if (k+1 >  nrows_R),
      break;
  end;

  sa = R(k,k);
  sb = R(k+1,k);
  if (sb == 0),
     continue;
  end;

  [c,s] = srotg( sa, sb);
  G = [c,  -s; ...
       s,   c];

  % ------------------------
  % apply to R(k:(k+1), k:n)
  % ------------------------
  R(k:(k+1), k:m)  = G * R(k:(k+1), k:m);

  % ------------------------
  % apply to Q(1:n, k:(k+1))
  % ------------------------
  Q(1:n, k:(k+1)) = Q(1:n, k:(k+1)) * G';
end;


if (idebug >= 1),
  % --------------------------------
  % just for debugging  and checking
  % --------------------------------
  tol = 10^(-7);

  % -----------------------------------
  % double check, R is upper triangular
  % -----------------------------------
  err_R = norm( tril(R,-1) , 1);
  is_ok_R = (err_R < tol );
  if (~is_ok_R),
    R
    error(sprintf('qrupdate_rank_1: R is not upper triangular, err_R=%g', ...
                  err_R));

    return;
  end;

  % ----------------------------
  % double check Q is orthogonal
  % ----------------------------
  err_Q = norm(Q'*Q - eye(ncols_Q,ncols_Q),1);
  is_ok_Q = (err_Q < tol );
  if (~is_ok_Q),
    error(sprintf('qrupdate_rank_1: Q is not orthogonal, err_Q=%g', ...
                   err_Q));
    return;
  end;

  if (idebug >= 2),
    A = Qin * Rin;
    err_A = norm( (A + u*v') - Q*R,1);
    is_ok_A = (err_A < tol );
    if (~is_ok_A),
      error(sprintf('qrupdate_rank_1: invalid QR, err_A=%g', ...
                     err_A));
    end;
   end;
end;


