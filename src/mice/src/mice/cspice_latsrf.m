function [srfpts] = cspice_latsrf( method, target, et, fixref, lonlat )
   switch nargin
      case 5
         method = zzmice_str(method);
         target = zzmice_str(target);
         et     = zzmice_dp(et);
         fixref = zzmice_str(fixref);
         lonlat = zzmice_dp(lonlat);
      otherwise
         error ( ['Usage: [srfpts] = cspice_latsrf( `method`, `target`, ' ...
                                                'et, `fixref`, lonlat )'] )
   end
   try
      [srfpts] = mice( 'latsrf_c', method, target, et, fixref, lonlat );
   catch spiceerr
      rethrow(spiceerr)
   end