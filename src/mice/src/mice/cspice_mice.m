function [value] = cspice_mice( item )
   switch nargin
      case 1
         item = zzmice_str(item);
      otherwise
         error ( 'Usage: [`value`] = cspice_mice( `item` )' )
   end
   try
      [value] = mice( 'cspice_mice', item );
   catch spiceerr
      rethrow(spiceerr)
   end