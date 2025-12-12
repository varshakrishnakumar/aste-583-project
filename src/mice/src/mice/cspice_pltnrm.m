function [normal] = cspice_pltnrm( v1, v2, v3)
   switch nargin
      case 3
         v1    = zzmice_dp(v1);
         v2    = zzmice_dp(v2);
         v3    = zzmice_dp(v3);
      otherwise
         error ( 'Usage: [normal] = cspice_pltnrm( v1, v2, v3)' )
   end
   try
      [normal] = mice( 'pltnrm_c', v1, v2, v3);
   catch spiceerr
      rethrow(spiceerr)
   end