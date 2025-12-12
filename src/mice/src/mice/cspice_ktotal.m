function [count] = cspice_ktotal( kind )
   switch nargin
      case 1
         kind = zzmice_str(kind);
      otherwise
         error ( 'Usage: count = cspice_ktotal(`kind`)' )
   end
   try
      [count] = mice( 'ktotal_c', kind );
      count = zzmice_dp(count);
   catch spiceerr
      rethrow(spiceerr)
   end