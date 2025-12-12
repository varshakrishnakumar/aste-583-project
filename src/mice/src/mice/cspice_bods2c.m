function [code,found] = cspice_bods2c(name)
   switch nargin
      case 1
         name = zzmice_str(name);
      otherwise
         error ( 'Usage: [_code_, _found_] = cspice_bods2c(_`name`_)' )
   end
   try
      [bods2c] = mice('bods2c_s',name);
      code     = reshape( [bods2c.code], 1, [] );
      found    = reshape( [bods2c.found], 1, [] );
   catch spiceerr
      rethrow(spiceerr)
   end