function [frname] = cspice_frmnam( frcode )
   switch nargin
      case 1
         frcode = zzmice_int(frcode);
      otherwise
         error( 'Usage: [_`frname`_] = cspice_frmnam(_frcode_)' )
   end
   try
      [frname] = mice('frmnam_c', frcode);
   catch spiceerr
      rethrow(spiceerr)
   end