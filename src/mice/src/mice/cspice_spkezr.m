function [starg, lt] = cspice_spkezr(targ, et, ref, abcorr, obs)
   switch nargin
      case 5
         targ   = zzmice_str(targ);
         et     = zzmice_dp(et);
         ref    = zzmice_str(ref);
         abcorr = zzmice_str(abcorr);
         obs    = zzmice_str(obs);
      otherwise
         error ( ['Usage: [_starg(6)_, _lt_] = ' ...
                  'cspice_spkezr( `targ`, _et_, `ref`, `abcorr`, `obs`)'] )
   end
   try
      [starg_s] = mice('spkezr_s',targ,et,ref,abcorr,obs);
      starg     = reshape( [starg_s.state], 6, [] );
      lt        = reshape( [starg_s.lt   ], 1, [] );
   catch spiceerr
      rethrow(spiceerr)
   end