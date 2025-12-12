function [r] = cspice_q2m( q )
   switch nargin
      case 1
         q = zzmice_dp(q);
      otherwise
         error ( 'Usage: [_r(3,3)_] = cspice_q2m( _q(4)_ )' )
   end
   try
      [r] = mice('q2m_c', q );
   catch spiceerr
      rethrow(spiceerr)
   end