function  [meas, avg, stddev, idxsml, idxlon] = cspice_wnsumd( window )
   switch nargin
      case 1
         window = zzmice_win( window );
      otherwise
         error( [ 'Usage: [meas, avg, stddev, idxsml, idxlon] ' ...
                  ' = cspice_wnsumd( window )' ] )
      end
   try
      [wnsumd] = mice( 'wnsumd_s',  [zeros(6,1);window] );
      meas     = reshape( [wnsumd.meas],   1, [] );
      avg      = reshape( [wnsumd.avg],    1, [] );
      stddev   = reshape( [wnsumd.stddev], 1, [] );
      idxsml   = reshape( [wnsumd.idxsml], 1, [] );
      idxlon   = reshape( [wnsumd.idxlon], 1, [] );
   catch spiceerr
      rethrow(spiceerr)
   end