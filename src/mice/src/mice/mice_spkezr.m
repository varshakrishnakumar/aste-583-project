function [starg] = mice_spkezr(targ, et, ref, abcorr, obs)
   switch nargin
      case 5
         targ   = zzmice_str(targ);
         et     = zzmice_dp(et);
         ref    = zzmice_str(ref);
         abcorr = zzmice_str(abcorr);
         obs    = zzmice_str(obs);
      otherwise
         error ( ['Usage: [_starg_] = '                                    ...
                  'mice_spkezr( `targ`, _et_, `ref`, `abcorr`, `obs`)'] )
   end
   try
      [starg] = mice('spkezr_s', targ, et, ref, abcorr, obs);
   catch spiceerr
      rethrow(spiceerr)
   end