function [r] = cspice_axisar( axis, angle)
   switch nargin
      case 2
         axis  = zzmice_dp(axis);
         angle = zzmice_dp(angle);
      otherwise
         error ( 'Usage: [r(3,3)] = cspice_axisar( axis(3), angle)' )
   end
   try
      [r] = mice( 'axisar_c', axis, angle );
   catch spiceerr
      rethrow(spiceerr)
   end