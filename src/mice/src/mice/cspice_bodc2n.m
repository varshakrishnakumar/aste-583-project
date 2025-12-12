function [name, found] = cspice_bodc2n(code)
   switch nargin
      case 1
         code = zzmice_int(code);
      otherwise
         error ( 'Usage: [_`name`_, _found_] = cspice_bodc2n(_code_)' )
   end
   try
      [bodc2n] = mice('bodc2n_s', code);
      name     = char( bodc2n.name );
      found    = reshape( [bodc2n.found], 1, [] );
   catch spiceerr
      rethrow(spiceerr)
   end