function [ivals, found] = cspice_gipool( name, start, room )
   switch nargin
      case 3
         name  = zzmice_str(name);
         start = zzmice_int(start);
         room  = zzmice_int(room);
      otherwise
         error ( [ 'Usage: [ivals(), found] = '...
                   'cspice_gipool( `name`, start, room )' ] )
   end
   try
      [ivals, found] = mice( 'gipool_c', name, start, room);
      ivals = zzmice_dp(ivals);
      found = zzmice_logical(found);
   catch spiceerr
      rethrow(spiceerr)
   end