function [normal, point] = cspice_pl2nvp( plane )
   switch nargin
      case 1
         plane = zzmice_pln( plane );
      otherwise
         error ( ['Usage: [normal(3), point(3)] = cspice_pl2nvp( plane )'] )
   end
   try
      [normal, point] = mice('pl2nvp_c', plane );
   catch spiceerr
      rethrow(spiceerr)
   end