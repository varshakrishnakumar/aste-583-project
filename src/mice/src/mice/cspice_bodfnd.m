function [bodfnd] = cspice_bodfnd( body, item )
   switch nargin
      case 2
         body   = zzmice_int(body);
         item   = zzmice_str(item);
      otherwise
         error ( ['Usage: [bodfnd] = cspice_bodfnd( body, `item` )' ] )
   end
   try
      [bodfnd] = mice( 'bodfnd_c', body, item );
      [bodfnd] = zzmice_logical(bodfnd);
   catch spiceerr
      rethrow(spiceerr)
   end