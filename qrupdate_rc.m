function [Q,R] = qrupdate_rc(ksave, Qin,Rin)
%
% [Q,R] = qrupdate_rc(ksave, Qin,Rin)
% update QR factorization with extra row, extra col 
%
% on input
%
% Qin(1:ksave,1:ksave), R(1:ksave,1:ksave)
% is QR factorization
%
% Rin is augmented with new row and new column
%
% Rin(1:ksave,  ksave+1) is new column
% Rin(ksave+1, 1:ksave)  is new row
% Rin( ksave+1, ksave+1) is new diagonal
%
idebug = 0;

% -------------------------------
% use vector notation or for loop
% -------------------------------
use_vec = 0;

ksize = ksave + 1;

% ----------------------------------------   
% may reuse storage for Qin, Rin in C code
% ----------------------------------------   

R = zeros(ksize,ksize);
R(1:ksize,1:ksize) = Rin(1:ksize,1:ksize);

Q = eye( ksize,ksize);
Q(1:ksave,1:ksave) = Qin(1:ksave,1:ksave);
Q(ksize,1:ksave) = 0;
Q(1:ksave,ksize) = 0;
Q(ksize,ksize) = 1;

if (idebug >= 1),
  err = norm( Q(1:ksize,1:ksize)'*Q(1:ksize,1:ksize) - eye(ksize,ksize),1);
  disp(sprintf('qrupdate_rc: norm(Q*Q-eye) %g ', err ));
end;

if (idebug >= 1),
   % ----------------------
   % just for debug testing
   % ----------------------
   [Qtmp,Rtmp] = qr( Rin(1:ksize,1:ksize) );
   Q(1:ksize,1:ksize) = Q(1:ksize,1:ksize)*Qtmp(1:ksize,1:ksize)';
   R(1:ksize,1:ksize) = Rtmp(1:ksize,1:ksize);
   return;
end;

%
% ---------------------------------------
% column (ksave+1) in R 
% is modified by Q'
% ---------------------------------------
R(1:ksave, ksave+1) = Qin(1:ksave,1:ksave)'*Rin(1:ksave,ksave+1);

  %  ------------------------
  %  r11  r12  r13 c14
  %       r22  r23 c24
  %            r33 c34
  %  c41  c42  c43 c44
  %  ------------------------

k = ksave;
kp1 = k + 1;
if (use_vec),
  r1 = zeros(kp1,1);
  r2 = zeros(kp1,1);
  q1 = zeros(kp1,1);
  q2 = zeros(kp1,1);
end;

for j=1:k,
    % ---------------------
    % setup Givens rotation
    % ---------------------
    sa = R(j,j);
    sb = R(kp1,j);

    is_real = ((imag(sa) == 0) && (imag(sb) == 0));
    if (is_real),
      [c,s] = srotg( sa, sb );
      G = [c, -s; ...
           s,  c];
    else
      [c,s] = zrotg( sa, sb);
      G = [c,       -s; ...
           conj(s), c];
    end;
 
    % ---------------------
    % apply Givens rotation
    % ---------------------
    % [ c  -s] * [sa] = [dnorm]
    % [ s   c]   [sb]   [0]
    %
    % where dnorm = norm( [sa,sb],2)
    % ------------------------------
    
    if (use_vec),
      r1(j:kp1) = R(j,  j:kp1);
      r2(j:kp1) = R(kp1,j:kp1);
      R(j,  j:kp1) =   G(1,1)  * r1(j:kp1) + G(1,2) * r2(j:kp1);
      R(kp1,j:kp1) =   G(2,1)  * r1(j:kp1) + G(2,2) * r2(j:kp1);
    else
      for i=j:kp1,
          r_j_i   = R(j,i);
          r_kp1_i = R(kp1,i);
          R(j,  i) = G(1,1) * r_j_i + G(1,2) * r_kp1_i;
          R(kp1,i) = G(2,1) * r_j_i + G(2,2) * r_kp1_i;
      end;
    end;

    % -----------------------------------------
    % transpose of Givens apply to columns of Q
    % [q1 q2] * Gt,   where Gt = G'
    % -----------------------------------------
    Gt = G';
    if (use_vec),
      q1(1:kp1) = Q(1:kp1,j);
      q2(1:kp1) = Q(1:kp1,kp1);
      Q(1:kp1,  j) = q1(1:kp1) * Gt(1,1) + q2(1:kp1) * Gt(2,1);
      Q(1:kp1,kp1) = q1(1:kp1) * Gt(1,2) + q2(1:kp1) * Gt(2,2);
    else
      for i=1:kp1,
        q_i_j    = Q(i,j);
        q_i_kp1  = Q(i,kp1);
        Q(i,j)   = q_i_j * Gt(1,1) + q_i_kp1 * Gt(2,1);
        Q(i,kp1) = q_i_j * Gt(1,2) + q_i_kp1 * Gt(2,2);
      end;
    end;

  end;

end;

