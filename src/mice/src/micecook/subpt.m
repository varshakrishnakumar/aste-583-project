function subpt()
   SPICETRUE   = logical(1);
   SPICEFALSE  = logical(0);
   abcorr      = 'LT+S';
   answer      = 'n';
   disp ( ' '                                                        )
   disp ( '             Welcome to SUBPT'                            )
   disp ( ' '                                                        )
   disp ( 'This program demonstrates the use of CSPICE in computing' )
   disp ( 'the apparent sub-observer point on a target body. The'    )
   disp ( 'computations use light time and stellar aberration'       )
   disp ( 'corrections.'                                             )
   disp ( ' '                                                        )
   disp( ' ' )
   leap = input( 'Enter the name of a leapseconds kernel file: ', 's');
   cspice_furnsh( leap )
   disp( ' ' )
   pck = input('Enter the name of a planetary constants kernel: ', 's' );
   cspice_furnsh( pck )
   disp( ' ' )
   spk = input( 'Enter the name of a binary SPK ephemeris file: ', 's');
   cspice_furnsh( spk )
   disp( ' ' )
   disp( 'Working ... Please wait' )
   disp( ' ' )
   cont = SPICETRUE;
   while ( cont == SPICETRUE )
      obs = input( 'Enter the name of the observing body: ', 's');
      disp( ' ' )
      targ = input('Enter the name of a target body: ', 's');
      disp( ' ' )
      fixfrm = input('Enter the name of the target body-fixed frame: ', 's');
      disp( ' ' )
      maxpts = input( 'Enter the number of points to calculate: ' );
      disp( ' ' )
      if ( maxpts <= 0 )
         maxpts = 1;
      end
      if ( maxpts == 1 )
         utcbeg = input( 'Enter the UTC time: ', 's');
         disp(' ')
         etbeg = cspice_str2et( utcbeg );
         delta = 0.;
         epoch  = etbeg;
      else
         utcbeg = input( 'Enter the beginning UTC time: ', 's');
         disp(' ')
         utcend = input( 'Enter the ending UTC time: ', 's');
         disp(' ')
         etbeg = cspice_str2et( utcbeg );
         etend = cspice_str2et( utcend );
         delta  = ( etend - etbeg ) / (maxpts - 1. );
         epoch = [0:(maxpts-1)]*delta + etbeg;
      end
      disp( 'Planetocentric coordinates for the nearest point' )
      disp( 'on the target body to the observing body (deg).'  )
      txt = sprintf( 'Target body: %s          Observing body: %s', ...
                                                              targ, ...
                                                              obs );
      disp( txt )
      disp( ' ' )
      disp( '       UTC Time            Lat         Lon')
      disp( '----------------------------------------------')
      npts   = 1;
      while ( npts <= maxpts )
         [spoint, trgepc, srfvec] = cspice_subpnt( 'Near point: ellipsoid', ...
                                                  targ, epoch(npts),        ...
                                                  fixfrm, abcorr, obs );
         [ radius, lon, lat] = cspice_reclat( spoint );
         lon = lon * cspice_dpr;
         lat = lat * cspice_dpr;
         utcout = cspice_et2utc( epoch(npts), 'C', 3 );
         txt = sprintf ('  %.20s  %9.5f    %9.5f', utcout, ...
                                                   lat,    ...
                                                   lon  );
         disp( txt )
         npts  = npts + 1;
      end
      disp( ' ' )
      answer = input( 'Continue? (Enter Y or N): ', 's');
      if ( strcmp( 'N', answer) || strcmp( 'n', answer) )
         cont = SPICEFALSE;
      end
   end
   cspice_kclear