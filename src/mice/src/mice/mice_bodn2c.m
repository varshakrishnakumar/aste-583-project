function [code] = mice_bodn2c(name)
   switch nargin
      case 1
         name = zzmice_str(name);
      otherwise
         error ( 'Usage: [_code_] = mice_bodn2c(_`name`_)' )
   end
   try
      [code] = mice('bodn2c_s',name);
   catch spiceerr
      rethrow(spiceerr)
   end