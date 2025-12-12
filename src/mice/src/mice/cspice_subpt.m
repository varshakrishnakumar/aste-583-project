function [spoint, alt] = cspice_subpt( method, target, et, abcorr, obsrvr )
   switch nargin
      case 5
         method = zzmice_str(method);
         target = zzmice_str(target);
         et     = zzmice_dp(et);
         abcorr = zzmice_str(abcorr);
         obsrvr = zzmice_str(obsrvr);
      otherwise
         error ( ['Usage: [_spoint(3)_, _alt_] = '     ...
                  'cspice_subpt( `method`, `target`, ' ...
                  '_et_, `abcorr`, `obsrvr`)']  )
   end
   try
      [subpt] = mice('subpt_s', method, target, et, abcorr, obsrvr);
      spoint   = reshape( [subpt.pos], 3, [] );
      alt      = reshape( [subpt.alt], 1, [] );
   catch spiceerr
      rethrow(spiceerr)
   end