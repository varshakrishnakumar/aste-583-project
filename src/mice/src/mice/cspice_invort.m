function [mit] = cspice_invort( m )
   switch nargin
      case 1
         m = zzmice_dp(m);
      otherwise
         error( 'Usage: [mit(3,3)] = cspice_rotmat( m(3,3) )' )
   end
   try
      [mit] = mice('invort_c', m );
   catch spiceerr
      rethrow(spiceerr)
   end