function [code,found] = cspice_bodn2c(name)
   switch nargin
      case 1
         name = zzmice_str(name);
      otherwise
         error ( 'Usage: [_code_, _found_] = cspice_bodn2c(_`name`_)' )
   end
   try
      [bodn2c] = mice('bodn2c_s',name);
      code     = reshape( [bodn2c.code],  1, [] );
      found    = reshape( [bodn2c.found], 1, [] );
   catch spiceerr
      rethrow(spiceerr)
   end