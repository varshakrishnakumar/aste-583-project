function [ticks] = cspice_sctiks(sc, clkstr)
   switch nargin
      case 2
         sc     = zzmice_int(sc);
         clkstr = zzmice_str(clkstr);
      otherwise
         error ( 'Usage: [_ticks_] = cspice_sctiks(sc, _`clkstr`_)' )
   end
   try
      [ticks] = mice('sctiks_c', sc, clkstr);
   catch spiceerr
      rethrow(spiceerr)
   end