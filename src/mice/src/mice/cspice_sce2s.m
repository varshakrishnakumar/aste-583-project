function [sclkch] = cspice_sce2s(sc, et)
   switch nargin
      case 2
         sc = zzmice_int(sc);
         et = zzmice_dp(et);
      otherwise
         error ( 'Usage: [_`sclkch`_] = cspice_sce2s(sc, _et_)' )
   end
   try
      [sclkch] = mice('sce2s_c',sc, et);
   catch spiceerr
      rethrow(spiceerr)
   end