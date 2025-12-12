function [point, span1, span2] = cspice_pl2psv( plane )
   switch nargin
      case 1
         plane = zzmice_pln( plane );
      otherwise
         error ( ['Usage: [point(3), span1(3), span2(3)] ' ...
                  '= cspice_pl2psv( plane )'] )
   end
   try
      [point, span1, span2] = mice('pl2psv_c', plane );
   catch spiceerr
      rethrow(spiceerr)
   end