function [c,s] = srotg(sa,sb)
% -------------------------
% [c,s] = srotg(sa,sb)
% setup Givens rotation
% ---------------------
% [ c  -s] * [sa] = [dnorm]
% [ s   c]   [sb]   [0]
%
% where dnorm = norm( [sa,sb],2)
% ------------------------------
%
% Note this is essentially xROTG in BLAS library
% -------------------------
use_hypot = 1;
if (use_hypot),
  % ----------------------------------------------------------------------
  % hypot is safe way of computing sqrt( sa*sa + sb*sb ) without overflow
  % ----------------------------------------------------------------------
  % r = hypot(sa,sb);
  scale = abs(sa) + abs(sb);
  r = scale * sqrt( (sa/scale)^2 + (sb/scale)^2 );
  c = sa/r;
  s = -sb/r;
else
 % ------------------------------------------------------------- 
 % IEEE copysign(x,y) function can be emulated as abs(x)*sign(y)
 % ------------------------------------------------------------- 
 if (sb == 0),
    c = copysign(1,sa);
    s = 0;
    r = abs(sa);
 elseif (sa == 0),
    c = 0;
    s = -copysign(1,sb);
    r = abs(sb);
 elseif (abs(sb) > abs(sa)),
    t = sa/sb;
    u = copysign(sqrt(1+t*t),sb);
    s = -1/u;
    c = -s*t;
    r = sb * u;
  else
    t = sb/sa;
    u = copysign(sqrt(1+t*t),sa);
    c = 1/u;
    s = -c*t;
    r = sa*u;
  end
end
