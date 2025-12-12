function [name] = cspice_bodc2s(code)
   switch nargin
      case 1
         code = zzmice_int(code);
      otherwise
         error ( 'Usage: [_`name`_] = cspice_bodc2s(_code_)' )
   end
   try
      [bodc2s] = mice('bodc2s_s', code);
      name     = char( bodc2s.name );
   catch spiceerr
      rethrow(spiceerr)
   end