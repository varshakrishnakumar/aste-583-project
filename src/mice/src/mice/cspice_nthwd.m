function [word, loc] = cspice_nthwd( string, nth )
   switch nargin
      case 2
         string = zzmice_str(string, true);
         nth    = zzmice_int(nth);
      otherwise
         error ( 'Usage: [`word`, loc] = cspice_nthwd( `string`, nth )' )
   end
   try
      [word, loc] = mice('nthwd_c', string, nth);
   catch spiceerr
      rethrow(spiceerr)
   end