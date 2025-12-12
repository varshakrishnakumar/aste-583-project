function [subslr] = mice_subslr( method, target, et, fixref, abcorr, obsrvr )
   switch nargin
      case 6
         method = zzmice_str(method);
         target = zzmice_str(target);
         et     = zzmice_dp(et);
         fixref = zzmice_str(fixref);
         abcorr = zzmice_str(abcorr);
         obsrvr = zzmice_str(obsrvr);
      otherwise
         error ( ['Usage: [_subslr_] = '                                   ...
                  'mice_subslr( `method`, `target`, _et_,'                 ...
                  ' `fixref`, `abcorr`, `obsrvr` )'])
   end
   try
      [subslr] = mice('subslr_s', method, target, et, fixref, abcorr, obsrvr);
   catch spiceerr
      rethrow(spiceerr)
   end