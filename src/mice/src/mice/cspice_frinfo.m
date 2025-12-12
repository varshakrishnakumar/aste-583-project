function [cent, frclss, clssid, found] = cspice_frinfo( frcode )
   switch nargin
      case 1
         frcode = zzmice_int(frcode);
      otherwise
         error( ['Usage: [_cent_, _frclss_, _clssid_, _found_] = ' ...
                                           'cspice_frinfo(_frcode_)'] )
   end
   try
      [frinfo] = mice('frinfo_s', frcode);
      cent   = reshape( [frinfo.center],   1, [] );
      frclss = reshape( [frinfo.class],    1, [] );
      clssid = reshape( [frinfo.class_ID], 1, [] );
      found  = reshape( [frinfo.found],    1, [] );
   catch spiceerr
      rethrow(spiceerr)
   end