function [frcode, frname, cent, found] = cspice_ccifrm( frclss, clssid )
   switch nargin
      case 2
         frclss = zzmice_int(frclss);
         clssid = zzmice_int(clssid);
      otherwise
         error( ['Usage: [ frcode, `frname`, cent, found] = ' ...
                'cspice_ccifrm( frclss, clssid)'] )
   end
   try
      [ccifrm, cent] = mice( 'ccifrm_s', frclss, clssid );
      frcode   = reshape( [ccifrm.code],  1, [] );
      frname   = char( ccifrm.name );
      found    = reshape( [ccifrm.found], 1, [] );
   catch spiceerr
      rethrow(spiceerr)
   end