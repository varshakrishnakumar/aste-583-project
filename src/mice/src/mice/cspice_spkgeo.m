function [state, lt] = cspice_spkgeo( targ, et, ref, obs )
   switch nargin
      case 4
         targ = zzmice_int(targ);
         et = zzmice_dp(et);
         ref = zzmice_str(ref);
         obs = zzmice_int(obs);
      otherwise
         error ( ['Usage: [state(6), lt] = ',               ...
                  'cspice_spkgeo( targ, et, `ref`, obs )'] )
   end
   try
      [state_s] = mice('spkgeo_s', targ, et, ref, obs);
      state  = reshape( [state_s.state ], 6, [] );
      lt     = reshape( [state_s.lt    ], 1, [] );
   catch spiceerr
      rethrow(spiceerr)
   end