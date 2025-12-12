function [cvals,found] = cspice_gnpool( name, start, room )
   switch nargin
      case 3
         name  = zzmice_str(name);
         start = zzmice_int(start);
         room  = zzmice_int(room);
      otherwise
         error ( ['Usage: [cvals(), found] = ' ...
                  'cspice_gnpool( `name`, start, room )' ] )
   end
   try
      [cvals, found] = mice( 'gnpool_c', name, start, room );
      found = zzmice_logical(found);
   catch spiceerr
      rethrow(spiceerr)
   end