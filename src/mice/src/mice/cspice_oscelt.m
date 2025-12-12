function [elts] = cspice_oscelt( state, et, mu )
   switch nargin
      case 3
         state = zzmice_dp(state);
         et    = zzmice_dp(et);
         mu    = zzmice_dp(mu);
      otherwise
         error( 'Usage: [_elts(8)_] = cspice_oscelt( _state(6)_, _et_, mu )' )
   end
   try
      [elts] = mice('oscelt_c', state, et, mu );
   catch spiceerr
      rethrow(spiceerr)
   end