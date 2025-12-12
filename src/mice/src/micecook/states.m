function states()
   format      = 'c';
   prec        = 0;
   maxpts      = 0;
   SPICETRUE   = logical(1);
   SPICEFALSE  = logical(0);
   disp( ' '                                                  )
   disp( '                Welcome to STATES'                  )
   disp( ' '                                                  )
   disp( 'This program demonstrates the use of NAIF S- and P-')
   disp( 'Kernel (SPK) files and subroutines by computing the')
   disp( 'state of a target body as seen from an observing'   )
   disp( 'body at a number of epochs within a given time'     )
   disp( 'interval.'                                          )
   disp( ' '                                                  )
   disp( ' ' )
   leap = input( 'Enter the name of a leapseconds kernel file: ', 's');
   disp( ' ' )
   cspice_furnsh( leap )
   spk = input( 'Enter the name of a binary SPK ephemeris file: ', 's');
   disp( ' ' )
   cspice_furnsh( spk  )
   obs = input( 'Enter the name of the observing body: ', 's');
   disp( ' ' )
   targ = input('Enter the name of a target body: ', 's');
   disp( ' ' )
   while ( maxpts <= 0 )
      maxpts = input( 'Enter the number of states to be calculated: ');
      if ( maxpts <= 0 )
         disp( 'The number of states must be greater than 0.')
      end
      disp( ' ' )
   end
   if ( maxpts == 1 )
      utcbeg = input('Enter the UTC time: ', 's');
      disp( ' ' )
   elseif ( maxpts > 1 )
      utcbeg = input( 'Enter the beginning UTC time: ', 's' );
      disp( ' ' )
      utcend = input( 'Enter the ending UTC time: ', 's' );
      disp( ' ' )
   end
   frame  = input( 'Enter the inertial reference frame (e.g.:J2000): ', 's' );
   disp( ' ' )
   disp( 'Type of correction                              Type of state    ')
   disp( '-------------------------------------------------------------    ')
   disp( '''LT+S''    Light-time and stellar aberration     Apparent state ')
   disp( '''LT''      Light-time only                       True state     ')
   disp( '''NONE''    No correction                         Geometric state')
   disp( ' ' )
   abcorr = input( 'Enter ''LT+S'', ''LT'', or ''NONE'': ', 's');
   disp( ' ' )
   disp( 'Working ... Please wait' )
   disp( ' ' )
   if ( maxpts == 1 )
      etbeg = cspice_str2et( utcbeg );
   elseif ( maxpts > 1 )
      etbeg = cspice_str2et ( utcbeg );
      etend = cspice_str2et ( utcend );
   end
   if ( maxpts > 1 )
      delta = ( etend - etbeg ) / ( maxpts - 1.);
      et    = [0:(maxpts-1)]*delta + etbeg;
   elseif( maxpts == 1 )
      et    = etbeg;
   end
   [state, lt] = cspice_spkezr( targ, et, frame, abcorr, obs );
   utc = cspice_et2utc( et, format, prec );
   cont = SPICETRUE;
   i    = 1;
   while ( (i <= maxpts)  &&  (cont == SPICETRUE) )
      txt = sprintf( 'For time %d of %d, the state of:', i, maxpts );
      disp( txt)
      txt = sprintf( 'Body            : %s', targ );
      disp( txt)
      txt = sprintf( 'Relative to body: %s', obs );
      disp( txt)
      txt = sprintf( 'In Frame        : %s', frame );
      disp( txt)
      txt = sprintf( 'At UTC time     : %s', utc(i,:) );
      disp( txt)
      disp( ' ' )
      disp('                 Position (km)              Velocity (km/s)'    )
      disp('            -----------------------     -----------------------')
      state_ele = state(:,i);
      txt = sprintf( '          X: %23.16e     %26.16e', state_ele(1), ...
                                                         state_ele(4) );
      disp( txt)
      txt = sprintf( '          Y: %23.16e     %26.16e', state_ele(2), ...
                                                         state_ele(5) );
      disp( txt)
      txt = sprintf( '          Z: %23.16e     %26.16e', state_ele(3), ...
                                                         state_ele(6) );
      disp( txt)
      txt = sprintf( '  MAGNITUDE: %23.16e     %26.16e', ...
                                   norm(state_ele(1:3)), ...
                                    norm(state_ele(4:6)) );
      disp( txt)
      disp( ' ' )
      if ( i < maxpts )
         disp( ' ' )
         answer = input( 'Continue? (Enter Y or N): ', 's');
      end
      if ( strcmp( 'N', answer) || strcmp( 'n', answer) )
         cont = SPICEFALSE;
      end
      i  = i + 1;
   end
   cspice_kclear