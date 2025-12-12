function [vout] = cspice_vprjp( vin, plane )
   switch nargin
      case 2
         vin   = zzmice_dp( vin );
         plane = zzmice_pln( plane );
      otherwise
         error ( 'Usage: [vout(3)] = cspice_vprjp( vin(3), plane )' )
   end
   try
      [vout] = mice('vprjp_c', vin, plane );
   catch spiceerr
      rethrow(spiceerr)
   end