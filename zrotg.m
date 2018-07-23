function [c,s] = zrotg(ca,cb)
% -------------------------
% [c,s] = zrotg(ca,bb)
% setup Givens rotation
% ---------------------
% 
% [ c        -s] * [ca] = [dnorm]
% [ conj(s)   c]   [cb]   [0]
%
% where dnorm = norm( [ca,cb],2)
% ------------------------------
% Note this is essentially xROTG in BLAS library
% -------------------------

if (abs(ca) == 0),
   c = 0;
   s = 1;
else
   scale = abs(ca) + abs(cb);
   % ---------------------------------------------------------
   % note avoid expensive complex division
   %
   % norm = scale * sqrt( abs(ca/scale)^2 + abs(cb/scale)^2 );
   % ---------------------------------------------------------
   inv_scale = 1/scale;   % use real arithmetic
   real_ca = real(ca); imag_ca = imag(ca);
   real_cb = real(cb); imag_cb = imag(cb);

   ca_div_scale = complex( real_ca*inv_scale, imag_ca*inv_scale );
   cb_div_scale = complex( real_cb*inv_scale, imag_cb*inv_scale );
   norm = scale * sqrt(   abs(ca_div_scale)^2 + abs(cb_div_scale)^2 );

   alpha  = ca/abs(ca);
   c = abs(ca)/norm;
   s = -alpha * conj(cb)/norm;
end;
