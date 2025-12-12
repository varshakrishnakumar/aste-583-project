function simple()
   SPICETRUE   = logical(1);
   SPICEFALSE  = logical(0);
   abcorr      = 'LT';
   answer      = 'n';
   MAXPTS      = 10;
   disp(' ')
   disp(' ')
   disp('                    Welcome to SIMPLE'                  )
   disp(' ')
   disp('This program calculates the angular separation of two'  )
   disp('target bodies as seen from an observing body.'          )
   disp(' ')
   disp('The angular separations are calculated for each of 10'  )
   disp('equally spaced times in a given time interval. A table' )
   disp('of the results is presented.')
   disp(' ')
   ref  = 'J2000';
   corr = 'LT+S';
   leap = input( 'Enter the name of a leapseconds kernel file: ', 's');
   cspice_furnsh( leap )
   disp( ' ' )
   spk = input( 'Enter the name of a binary SPK ephemeris file: ', 's');
   cspice_furnsh( spk )
   disp( ' ' )
   cont = SPICETRUE;
   while ( cont == SPICETRUE )
      obs = input( 'Enter the name of the observing body: ', 's');
      disp( ' ' )
      targ1 = input('Enter the name of the first target body: ', 's');
      disp( ' ' )
      targ2 = input('Enter the name of the second target body: ', 's');
      disp( ' ' )
      utcbeg = input( 'Enter the beginning UTC time: ', 's');
      disp(' ')
      utcend = input( 'Enter the ending UTC time: ', 's');
      disp(' ')
      disp( ' ' )
      disp( 'Working ... Please wait' )
      disp( ' ' )
      etbeg = cspice_str2et( utcbeg );
      etend = cspice_str2et( utcend );
      delta = ( etend - etbeg ) / (MAXPTS - 1. );
      et = [0:(MAXPTS-1)]*delta + etbeg;
      [ pos1, lt1] = cspice_spkpos( targ1, et, ref, corr, obs );
      [ pos2, lt2] = cspice_spkpos( targ2, et, ref, corr, obs );
      y = cspice_vsep( pos1, pos2);
      y = y * cspice_dpr;
      disp( ' ' )
      txt = sprintf( 'The angular separation between bodies %s and %s,', ...
                                                         targ1, targ2 );
      disp( txt )
      txt = sprintf( 'as seen from body %s.', obs );
      disp( txt )
      disp( ' ')
      utcbeg = cspice_et2utc( etbeg, 'C', 0 );
      txt = sprintf( 'From: %s', utcbeg );
      disp( txt )
      utcend = cspice_et2utc( etend, 'C', 0 );
      txt = sprintf( 'To  : %s', utcend );
      disp( txt )
      utctim = cspice_et2utc( et, 'C', 0 );
      disp( ' ' )
      disp( '       UTC Time                 Separation' )
      disp( '----------------------------------------------' )
      for i=1:MAXPTS
         txt = sprintf( '  %.20s  %15.8f deg', utctim(i,:), y(i) );
         disp( txt )
      end
      disp( ' ')
      answer = input( 'Continue? (Enter Y or N): ', 's');
      if ( strcmp( 'N', answer) || strcmp( 'n', answer) )
         cont = SPICEFALSE;
      end
   end
   cspice_kclear