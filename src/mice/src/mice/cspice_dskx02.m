function [ plid, xpt, found] = cspice_dskx02( handle, dladsc, vertex, raydir)
   switch nargin
      case 4
         handle = zzmice_int( handle );
         dladsc = zzmice_int( dladsc );
         vertex = zzmice_dp( vertex );
         raydir = zzmice_dp( raydir );
      otherwise
         error ( [ 'Usage: [ plid, xpt, found] = ' ...
                   'cspice_dskx02( handle, dladsc, vertex, raydir)' ] )
   end
   try
      [ plid, xpt, found] = mice( 'dskx02_c', ...
                                  handle, dladsc, vertex, raydir );
      found = zzmice_logical(found);
   catch spiceerr
      rethrow(spiceerr)
   end