function [starg, lt] = cspice_spkez( targ, et, ref, abcorr, obs )
   switch nargin
      case 5
         targ   = zzmice_int(targ);
         et     = zzmice_dp(et);
         ref    = zzmice_str(ref);
         abcorr = zzmice_str(abcorr);
         obs    = zzmice_int(obs);
      otherwise
         error ( ['Usage: [starg(6), lt] = cspice_spkez( targ, et,' ...
                  ' `ref`, `abcorr`, obs )'] )
   end
   try
      [starg_s] = mice('spkez_s', targ, et, ref, abcorr, obs);
      starg  = reshape( [starg_s.state ], 6, [] );
      lt     = reshape( [starg_s.lt    ], 1, [] );
   catch spiceerr
      rethrow(spiceerr)
   end