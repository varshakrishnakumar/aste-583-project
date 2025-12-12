function [nxpts, xpt] = cspice_inrypl( vertex, dir, plane )
   switch nargin
      case 3
         vertex  = zzmice_dp(vertex);
         dir     = zzmice_dp(dir);
         plane   = zzmice_pln(plane);
      otherwise
         error ( 'Usage: [nxpts, xpt] = cspice_inrypl( vertex, dir, plane )' )
   end
   try
      [nxpts, xpt] = mice( 'inrypl_c', vertex, dir, plane );
   catch spiceerr
      rethrow(spiceerr)
   end