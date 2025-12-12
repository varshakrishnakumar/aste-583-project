function [wnreld] = cspice_wnreld( a, op, b )
   switch nargin
      case 3
         a  = zzmice_win(a);
         op = zzmice_str(op);
         b  = zzmice_win(b);
      otherwise
         error( 'Usage: [wnreld] = cspice_wnreld( a, `op`, b )' )
      end
   try
      [wnreld] = mice( 'wnreld_c', [zeros(6,1); a], op, [zeros(6,1); b] );
      [wnreld] = zzmice_logical(wnreld);
   catch spiceerr
      rethrow(spiceerr)
   end