function [x, y, z] = cspice_frame( x )
   switch nargin
      case 1
         x = zzmice_dp( x );
      otherwise
         error ( 'Usage: [x(3), y(3), z(3)] = cspice_frame( x(3) )' )
   end
   try
      [x, y, z] = mice( 'frame_c', x );
   catch spiceerr
      rethrow(spiceerr)
   end