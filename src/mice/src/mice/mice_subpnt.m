function [subpnt] = mice_subpnt( method, target, et, fixref, abcorr, obsrvr )
   switch nargin
      case 6
         method = zzmice_str(method);
         target = zzmice_str(target);
         et     = zzmice_dp(et);
         fixref = zzmice_str(fixref);
         abcorr = zzmice_str(abcorr);
         obsrvr = zzmice_str(obsrvr);
      otherwise
         error ( ['Usage: [_subnt_] = '                                    ...
                  'mice_subpnt( `method`, `target`, _et_,'                 ...
                  ' `fixref`, `abcorr`, `obsrvr`)'])
   end
   try
      [subpnt] = mice('subpnt_s', method, target, et, fixref, abcorr, obsrvr);
   catch spiceerr
      rethrow(spiceerr)
   end