function [shape, frame, bsight, bounds] = cspice_getfvn( inst, room )
   switch nargin
      case 2
         inst = zzmice_str(inst);
         room = zzmice_int(room);
      otherwise
         error ( [ 'Usage: [`shape`, `frame`, bsight(3), bounds(3,N)] = '   ...
                   'cspice_getfvn( `inst`, room )' ] )
   end
   try
      [shape, frame, bsight, bounds] = mice('getfvn_c', inst, room);
   catch spiceerr
      rethrow(spiceerr)
   end