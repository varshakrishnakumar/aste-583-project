function [et] = cspice_scs2e(sc, sclkch)
   switch nargin
      case 2
         sc     = zzmice_int(sc);
         sclkch = zzmice_str(sclkch);
      otherwise
         error ( 'Usage: [_et_] = cspice_scs2e(sc, _`sclkch`_)' )
   end
   try
      [et] = mice('scs2e_c',sc, sclkch);
   catch spiceerr
      rethrow(spiceerr)
   end