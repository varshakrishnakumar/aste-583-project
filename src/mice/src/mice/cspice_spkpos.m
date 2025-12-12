function [ptarg, lt] = cspice_spkpos(targ, et, ref, abcorr, obs)
   switch nargin
      case 5
         targ   = zzmice_str(targ);
         et     = zzmice_dp(et);
         ref    = zzmice_str(ref);
         abcorr = zzmice_str(abcorr);
         obs    = zzmice_str(obs);
      otherwise
         error ( ['Usage: [_ptarg(3)_, _lt_] = ' ...
                  'cspice_spkpos( `targ`, _et_, `ref`, `abcorr`, `obs`)'] )
   end
   try
      [ptarg_s] = mice('spkpos_s', targ, et,ref, abcorr, obs);
      ptarg     = reshape( [ptarg_s.pos], 3, [] );
      lt        = reshape( [ptarg_s.lt],  1, [] );
   catch spiceerr
      rethrow(spiceerr)
   end