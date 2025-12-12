function [output] = cspice_timout( et, pictur )
   switch nargin
      case 2
         et     = zzmice_dp(et);
         pictur = zzmice_str(pictur);
      otherwise
         error( 'Usage: [_`output`_] = cspice_timout( _et_, `pictur` )' )
   end
   try
      [output] = mice('timout_c', et, pictur);
   catch spiceerr
      rethrow(spiceerr)
   end