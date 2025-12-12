function [normal, konst] = cspice_pl2nvc( plane )
   switch nargin
      case 1
         plane = zzmice_pln( plane );
      otherwise
         error( ['Usage: [normal(3), konst] = cspice_pl2nvc( plane )'] )
   end
   try
      [normal, konst] = mice('pl2nvc_c', plane );
   catch spiceerr
      rethrow(spiceerr)
   end