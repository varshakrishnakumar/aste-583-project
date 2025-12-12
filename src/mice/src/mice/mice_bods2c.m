function [code] = mice_bods2c(name)
   switch nargin
      case 1
         name = zzmice_str(name);
      otherwise
         error ( 'Usage: [_code_] = mice_bods2c(_`name`_)' )
   end
   try
      [code] = mice('bods2c_s',name);
   catch spiceerr
      rethrow(spiceerr)
   end