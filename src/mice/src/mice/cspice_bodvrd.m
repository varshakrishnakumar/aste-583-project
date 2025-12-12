function [values] = cspice_bodvrd(bodynm, item, maxn)
   switch nargin
      case 3
         bodynm = zzmice_str(bodynm);
         item   = zzmice_str(item);
         maxn   = zzmice_int(maxn);
      otherwise
         error ( ['Usage:  [values] = cspice_bodvrd( `bodynm`, ', ...
                                                    '`item`, ',   ...
                                                    'maxn )' ] )
   end
   try
      [values] = mice( 'bodvrd_c', bodynm, item, maxn);
   catch spiceerr
      rethrow(spiceerr)
   end