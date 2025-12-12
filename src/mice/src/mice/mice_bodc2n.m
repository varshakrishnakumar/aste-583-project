function [name] = mice_bodc2n(code)
   switch nargin
      case 1
         code = zzmice_int(code);
      otherwise
         error ( 'Usage: [_`name`_] = mice_bodc2n(_code_)' )
   end
   try
      [name] = mice('bodc2n_s', code);
   catch spiceerr
      rethrow(spiceerr)
   end