   geomError = !solver_.CalculateAngle(target_state_, angles);
   if( geomError)
   {
    //2. taking 0 and max, subtract 0.38 of its difference
    double goldB = 0.618, goldS = 0.382;
    int iteration = 10;
    math::PlatformState stateLeft;
    stateLeft.heihgt = state.heihgt;
    stateLeft.tangage = goldS * state.tangage;
    stateLeft.roll = goldS * state.roll;
    math::PlatformState stateRight;
    stateRight.heihgt = state.heihgt;
    stateRight.tangage = state.tangage - goldS * state.tangage;
    stateRight.roll = state.roll - goldS * state.roll;
    math::PlatformState stateBorderLeft;
    stateBorderLeft.heihgt = state.heihgt;
    stateBorderLeft.tangage = 0;
    stateBorderLeft.roll = 0;
    math::PlatformState stateBorderRight;
    stateBorderRight.heihgt = state.heihgt;
    stateBorderRight.tangage = state.tangage;
    stateBorderRight.roll = state.roll;
    double eps = 0.1;
    if( fabs( state.tangage) > fabs( state.roll))
    {
     eps = fabs( eps * state.tangage);
    }
    else
    {
     eps = fabs( eps * state.roll);
    }
    for( int i = iteration; i > 0; i--)
    {
     bool resultLeft = solver_.CalculateAngle( stateLeft, angles);
     bool resultRight = solver_.CalculateAngle( stateRight, angles);
     if( resultLeft && resultRight)
     {
      //3.Point is in right segment
      //3.1.Checking values of angles, whether difference of largest of them is small enough
      bool processNew = false;
      if( fabs( state.tangage) > fabs( state.roll))
      {
       if( fabs( fabs( stateBorderRight.tangage) - fabs( stateRight.tangage)) > eps)
       {
        processNew = true;
       }
      }
      else
      {
       if( fabs( fabs( stateBorderRight.roll) - fabs( stateRight.roll)) > eps)
       {
        processNew = true;
       }
      }
      if( processNew)
      {
       stateBorderLeft.tangage = stateRight.tangage;
       stateBorderLeft.roll = stateRight.roll;
       
       stateLeft.tangage = stateRight.tangage + ( stateBorderRight.tangage - stateRight.tangage) * goldS;
       stateLeft.roll = stateRight.roll + ( stateBorderRight.roll - stateRight.roll) * goldS;
       
       stateRight.tangage = stateBorderRight.tangage - ( stateBorderRight.tangage - stateRight.tangage) * goldS;
       stateRight.roll = stateBorderRight.roll - ( stateBorderRight.roll - stateRight.roll) * goldS;
       continue;
      }
      else
      {
       state.heihgt = stateRight.heihgt;
       state.tangage = stateRight.tangage;
       state.roll = stateRight.roll;
       break;
      }
     }
     else if( !resultLeft && !resultRight)
     {
      //4.Point is in left segment
      bool processNew = false;
      if( fabs( state.tangage) > fabs( state.roll))
      {
       if( fabs( fabs( stateBorderRight.tangage) - fabs( stateRight.tangage)) > eps)
       {
        processNew = true;
       }
      }
      else
      {
       if( fabs( fabs( stateBorderRight.roll) - fabs( stateRight.roll)) > eps)
       {
        processNew = true;
       }
      }
      if( processNew)
      {
       stateBorderRight.tangage = stateLeft.tangage;
       stateBorderRight.roll = stateLeft.roll;
       
       stateRight.tangage = stateLeft.tangage - ( stateLeft.tangage - stateBorderLeft.tangage) * goldS;
       stateRight.roll = stateLeft.roll - ( stateLeft.roll - stateBorderLeft.roll) * goldS;
  
       stateLeft.tangage = stateBorderLeft.tangage + ( stateLeft.tangage - stateBorderLeft.tangage) * goldS;
       stateLeft.roll = stateBorderLeft.roll + ( stateLeft.roll - stateBorderLeft.roll) * goldS;
       continue;
      }
      else
      {
       state.heihgt = stateLeft.heihgt;
       state.tangage = stateLeft.tangage;
       state.roll = stateLeft.roll;
       break;
      }
     }
     else
     {
      //5.Point is in mid segment
       stateBorderLeft.tangage = stateLeft.tangage;
       stateBorderLeft.roll = stateLeft.roll;
       
       stateBorderRight.tangage = stateRight.tangage;
       stateBorderRight.roll = stateRight.roll;
       
       stateLeft.tangage = stateBorderLeft.tangage + ( stateBorderRight.tangage - stateBorderLeft.tangage) * goldS;
       stateLeft.roll = stateBorderLeft.roll + ( stateBorderRight.roll - stateBorderLeft.roll) * goldS;
       
       stateRight.tangage = stateBorderRight.tangage - ( stateBorderRight.tangage - stateBorderLeft.tangage) * goldS;
       stateRight.roll = stateBorderRight.roll - ( stateBorderRight.roll - stateBorderLeft.roll) * goldS;
       continue;
     }
    }
    /*fprintf( stdout, "Final result is\n");
    fflush( stdout);*/
    solver_.CalculateAngle( state, angles);
   }