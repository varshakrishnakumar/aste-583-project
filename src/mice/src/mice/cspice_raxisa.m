function [axis, angle] = cspice_raxisa(matrix)
   switch nargin
      case 1
         matrix = zzmice_dp(matrix);
      otherwise
         error ( [ 'Usage: [ axis(3), angle] = ' ...
                   'cspice_raxisa( matrix(3,3) )' ] )
   end
   try
      [axis, angle] = mice('raxisa_c', matrix);
   catch spiceerr
      rethrow(spiceerr)
   end