function [wnelmd] = cspice_wnelmd( point, window )
   switch nargin
      case 2
         point  = zzmice_dp(point);
         window = zzmice_win(window);
      otherwise
         error( 'Usage: [wnelmd] = cspice_wnelmd( point, window )' )
      end
   try
      [wnelmd] = mice( 'wnelmd_c', point, [zeros(6,1); window] );
      [wnelmd] = zzmice_logical(wnelmd);
   catch spiceerr
      rethrow(spiceerr)
   end