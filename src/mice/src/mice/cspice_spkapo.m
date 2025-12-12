function [ptarg, lt] = cspice_spkapo( targ, et, ref, sobs, abcorr )
   switch nargin
      case 5
         targ = zzmice_int(targ);
         et = zzmice_dp(et);
         ref = zzmice_str(ref);
         sobs = zzmice_dp(sobs);
         abcorr = zzmice_str(abcorr);
      otherwise
         error ( ['Usage: [ptarg(3), lt] = ',                             ...
                  'cspice_spkapo( targ, et, `ref`, sobs(6), `abcorr` )'] )
   end
   try
      [ptarg_s] = mice('spkapo_s', targ, et, ref, sobs, abcorr);
      ptarg  = reshape( [ptarg_s.pos   ], 3, [] );
      lt     = reshape( [ptarg_s.lt    ], 1, [] );
   catch spiceerr
      rethrow(spiceerr)
   end