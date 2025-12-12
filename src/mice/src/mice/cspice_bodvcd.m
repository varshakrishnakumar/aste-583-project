function [values] = cspice_bodvcd (bodyid, item, maxn)
   switch nargin
      case 3
         bodyid = zzmice_int(bodyid);
         item   = zzmice_str(item);
         maxn   = zzmice_int(maxn);
      otherwise
         error ( 'Usage:  [values] = cspice_bodvcd(bodyid, `item`, maxn)' )
   end
   try
      [values] = mice( 'bodvcd_c', bodyid, item, maxn);
   catch spiceerr
      rethrow(spiceerr)
   end