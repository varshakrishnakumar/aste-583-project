function [frcode] = cspice_namfrm(frname)
   switch nargin
      case 1
         frname = zzmice_str(frname);
      otherwise
         error ( 'Usage: [_frcode_] = cspice_namfrm(_`frname`_)' )
   end
   try
      [frcode] = mice('namfrm_c',frname);
   catch spiceerr
      rethrow(spiceerr)
   end