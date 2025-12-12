function [pltar] = cspice_pltar( vrtces, plates )
   switch nargin
      case 2
         vrtces  = zzmice_dp(vrtces);
         plates  = zzmice_int(plates);
      otherwise
         error ( 'Usage: [pltar] = cspice_pltar( vrtces, plates )' )
   end
   try
      [pltar] = mice('pltar_c', vrtces, plates);
   catch spiceerr
      rethrow(spiceerr)
   end