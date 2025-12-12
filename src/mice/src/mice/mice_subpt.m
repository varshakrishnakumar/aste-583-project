function [spoint] = mice_subpt( method, target, et, abcorr, obsrvr )
   switch nargin
      case 5
         method = zzmice_str(method);
         target = zzmice_str(target);
         et     = zzmice_dp(et);
         abcorr = zzmice_str(abcorr);
         obsrvr = zzmice_str(obsrvr);
      otherwise
         error ( ['Usage: [_spoint_] = mice_subpt( `method`, `target`,'    ...
                  ' _et_, `abcorr`, `obsrvr`)'] )
   end
   try
      [spoint] = mice('subpt_s', method, target, et, abcorr, obsrvr);
   catch spiceerr
      rethrow(spiceerr)
   end