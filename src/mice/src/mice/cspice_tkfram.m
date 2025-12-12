function [rot, frame, found] = cspice_tkfram( frcode )
   switch nargin
      case 1
         frcode = zzmice_int(frcode);
      otherwise
         error ( 'Usage: [rot(3,3), frame, found] = cspice_tkfram( frcode )' )
   end
   try
      [rot, frame, found] = mice('tkfram_c', frcode);
      found = zzmice_logical(found);
   catch spiceerr
      rethrow(spiceerr)
   end