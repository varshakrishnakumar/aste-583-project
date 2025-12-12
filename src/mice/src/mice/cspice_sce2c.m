function [sclkdp] = cspice_sce2c(sc, et)
   switch nargin
      case 2
         sc = zzmice_int(sc);
         et = zzmice_dp(et);
      otherwise
         error ( 'Usage: [_sclkdp_] = cspice_sce2c(sc, _et_)' )
   end
   try
      [sclkdp] = mice('sce2c_c',sc, et);
   catch spiceerr
      rethrow(spiceerr)
   end