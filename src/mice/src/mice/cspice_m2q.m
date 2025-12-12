function [q] = cspice_m2q(r)
   switch nargin
      case 1
         r = zzmice_dp(r);
      otherwise
         error ( 'Usage: [_q(4)_] = cspice_m2q( _r(3,3)_ )' )
   end
   try
      [q] = mice('m2q_c', r);
   catch spiceerr
      rethrow(spiceerr)
   end