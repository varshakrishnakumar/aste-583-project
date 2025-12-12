function [ellips, found] = cspice_inedpl( a, b, c, plane )
   switch nargin
      case 4
         a     = zzmice_dp(a);
         b     = zzmice_dp(b);
         c     = zzmice_dp(c);
         plane = zzmice_pln(plane);
      otherwise
         error ( 'Usage: [ellips, found] = cspice_inedpl( a, b, c, plane )' )
   end
   try
      [ellips, found] = mice( 'inedpl_c', a, b, c, plane );
      found = zzmice_logical(found);
   catch spiceerr
      rethrow(spiceerr)
   end