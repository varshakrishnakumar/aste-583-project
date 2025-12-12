function [pstart, pstop] = cspice_scpart( sc )
   switch nargin
      case 1
         sc = zzmice_int(sc);
      otherwise
         error ( [ 'Usage: [pstart(nparts), pstop(nparts)] = '              ...
                   'cspice_scpart( sc )' ] )
   end
   try
      [pstart, pstop] = mice('scpart_c', sc);
   catch spiceerr
      rethrow(spiceerr)
   end