function [clkstr] = cspice_scfmt( sc, ticks )
   switch nargin
      case 2
         sc    = zzmice_int(sc);
         ticks = zzmice_dp(ticks);
      otherwise
         error ( 'Usage: [`clkstr`] = cspice_scfmt( sc, ticks )' )
   end
   try
      [clkstr] = mice('scfmt_c', sc, ticks);
   catch spiceerr
      rethrow(spiceerr)
   end