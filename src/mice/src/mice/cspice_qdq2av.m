function [av] = cspice_qdq2av( q, dq )
   switch nargin
      case 2
         q = zzmice_dp(q);
         dq = zzmice_dp(dq);
      otherwise
         error ( 'Usage: [av(3)] = cspice_qdq2av( q(4), dq(4) )' )
   end
   try
      [av] = mice('qdq2av_c', q, dq);
   catch spiceerr
      rethrow(spiceerr)
   end