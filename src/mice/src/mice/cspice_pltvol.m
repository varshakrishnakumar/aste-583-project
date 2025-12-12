function [pltvol] = cspice_pltvol( vrtces, plates )
   switch nargin
      case 2
         vrtces  = zzmice_dp(vrtces);
         plates  = zzmice_int(plates);
      otherwise
         error ( 'Usage: [pltvol] = cspice_pltvol( vrtces, plates )' )
   end
   try
      [pltvol] = mice('pltvol_c', vrtces, plates);
   catch spiceerr
      rethrow(spiceerr)
   end