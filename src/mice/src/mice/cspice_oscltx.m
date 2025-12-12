function [elts] = cspice_oscltx( state, et, mu )
   switch nargin
      case 3
         state = zzmice_dp(state);
         et = zzmice_dp(et);
         mu = zzmice_dp(mu);
      otherwise
         error ( ['Usage: [elts(SPICE_OSCLTX_NELTS)] = ',                  ...
                                    'cspice_oscltx( state(6), et, mu )'] )
   end
   try
      [elts] = mice('oscltx_c', state, et, mu);
   catch spiceerr
      rethrow(spiceerr)
   end