function [et] = cspice_sct2e(sc,sclkdp)
   switch nargin
      case 2
         sc     = zzmice_int(sc);
         sclkdp = zzmice_dp(sclkdp);
      otherwise
         error( 'Usage: [_et_] = cspice_sct2e(sc, _sclkdp_)' )
   end
   try
      [et] = mice('sct2e_c',sc, sclkdp);
   catch spiceerr
      rethrow(spiceerr)
   end