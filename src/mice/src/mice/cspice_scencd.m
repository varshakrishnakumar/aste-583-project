function [sclkdp] = cspice_scencd(sc, sclkch)
   switch nargin
      case 2
         sc     = zzmice_int(sc);
         sclkch = zzmice_str(sclkch);
      otherwise
         error ( 'Usage: [_sclkdp_] = cspice_scencd(sc, _`sclkch`_)' )
   end
   try
      [sclkdp] = mice('scencd_c',sc, sclkch);
   catch spiceerr
      rethrow(spiceerr)
   end