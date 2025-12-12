function [state] = cspice_conics( elts, et )
   switch nargin
      case 2
         elts = zzmice_dp(elts);
         et   = zzmice_dp(et);
      otherwise
         error ( 'Usage: [_state(6)_] = cspice_conics( _elts(8)_, _et_ )' )
   end
   try
      [state] = mice('conics_c', elts, et);
   catch spiceerr
      rethrow(spiceerr)
   end