function [next, rest] = cspice_nextwd( string )
   switch nargin
      case 1
         string = zzmice_str(string, true);
      otherwise
         error ( 'Usage: [`next`, `rest`] = cspice_nextwd( `string` )' )
   end
   try
      [next, rest] = mice('nextwd_c', string);
   catch spiceerr
      rethrow(spiceerr)
   end