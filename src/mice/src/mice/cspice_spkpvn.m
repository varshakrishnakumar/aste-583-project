function [ref, state, center] = cspice_spkpvn( handle, descr, et)
   switch nargin
      case 3
         handle = zzmice_int(handle);
         descr  = zzmice_dp(descr);
         et     = zzmice_dp(et);
      otherwise
         error ( ['Usage: [ref, state(6), center] = ' ...
                          'cspice_spkpvn( handle, descr(5), et)'] )
   end
   try
      spkpvn = mice( 'spkpvn_s', handle, descr, et );
      ref    = reshape( [spkpvn.ref   ], 1, [] );
      state  = reshape( [spkpvn.state ], 6, [] );
      center = reshape( [spkpvn.center], 1, [] );
   catch spiceerr
      rethrow(spiceerr)
   end