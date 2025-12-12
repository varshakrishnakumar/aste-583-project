function [shape, frame, bsight, bounds] = cspice_getfov( instid, room )
   switch nargin
      case 2
         instid = zzmice_int( instid );
         room   = zzmice_int( room   );
      otherwise
         error ( ['Usage: [`shape`, `frame`, bsight(3), bounds(3,N)] = ' ...
                          'cspice_getfov( instid, room )']  )
   end
   try
      [shape, frame, bsight, bounds] = mice( 'getfov_c', instid, room );
   catch spiceerr
      rethrow(spiceerr)
   end