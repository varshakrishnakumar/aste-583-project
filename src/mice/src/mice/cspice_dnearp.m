function [dnear, dalt, found] = cspice_dnearp( state, a, b, c )
   switch nargin
      case 4
         state = zzmice_dp(state);
         a     = zzmice_dp(a);
         b     = zzmice_dp(b);
         c     = zzmice_dp(c);
      otherwise
         error ( [ 'Usage: [dnear(6), dalt(2), found] = '                   ...
                   'cspice_dnearp( state(6), a, b, c )' ] )
   end
   try
      [dnear, dalt, found] = mice('dnearp_c', state, a, b, c);
      found = zzmice_logical(found);
   catch spiceerr
      rethrow(spiceerr)
   end