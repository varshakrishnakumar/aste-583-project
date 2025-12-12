function [qout] = cspice_qxq( q1, q2 )
   switch nargin
      case 2
         q1 = zzmice_dp(q1);
         q2 = zzmice_dp(q2);
      otherwise
         error ( 'Usage: [qout(4)] = cspice_qxq( q1(4), q2(4) )' )
   end
   try
      [qout] = mice('qxq_c', q1, q2);
   catch spiceerr
      rethrow(spiceerr)
   end