function [sclkch] = cspice_scdecd(sc, sclkdp)
   switch nargin
      case 2
        sc     = zzmice_int(sc);
        sclkdp = zzmice_dp(sclkdp);
      otherwise
         error ( 'Usage: [_`sclkch`_] = cspice_scdecd(sc, _sclkdp_)' )
   end
   try
      [sclkch] = mice('scdecd_c',sc,sclkdp);
   catch spiceerr
      rethrow(spiceerr)
   end