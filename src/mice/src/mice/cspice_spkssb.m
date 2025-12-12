function [starg] = cspice_spkssb( targ, et, ref )
   switch nargin
      case 3
         targ = zzmice_int(targ);
         et = zzmice_dp(et);
         ref = zzmice_str(ref);
      otherwise
         error ( 'Usage: [starg(6)] = cspice_spkssb( targ, et, `ref` )' )
   end
   try
      [starg] = mice('spkssb_c', targ, et, ref);
   catch spiceerr
      rethrow(spiceerr)
   end