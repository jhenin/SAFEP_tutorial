# genLambdas
# Generates a linear sequence of num values between 0 and 1 inclusive
# Arguments: num, an integer
# Results: returns a string of the form "0 x1 x2 x3 x4 ... xnum=1"
proc genLambdas { num } {
    set lambdaSched ""
    for {set k 0} {$k <= $num} {set k [expr $k+1.0]} {
      append lambdaSched "[expr $k/$num] "
    }
    return $lambdaSched
}

# setTI
# Takes some basic arguments and builds a changing bias for DBC TI calculations
# Arguments:
#    cvName: the name of the colvar to be biased
#    biasType: the type of the bias
#    forceConst0: the force constant at lambda=0. (forceConstant)
#    forceConst1: the force constant at lambda=1. (targetForceConstant)
#    upperWalls: the upperWall of the restraint. 
#    nWindows: the number of lambda windows
#    equilSteps: the number of steps before computing TI gradients
#    stepsPerWindow: the number of steps for each lambda value
#    releaseFlag: a flag that indicates whether the restraint is being imposed 
#        or released. 0=imposed (force increased over time) 
#        1=released (force decreases over time)
# Results:
#    Creates a new bias based on the inputs and passes it to the Colvars module

proc setTI { cvName biasType forceConst0 forceConst1 forceExp upperWalls nWindows equilSteps stepsPerWindow releaseFlag} {
	
    set lambdaSched [genLambdas $nWindows]
    
    if { $releaseFlag == "True" } {
    	puts "Scaling $cvName from $forceConst1 to $forceConst0 over $nWindows windows."
    	set lambdaSched [lreverse $lambdaSched]
    }
    
    set TIbias "$biasType {\n \
    	colvars               $cvName \n \  
	targetForceConstant   $forceConst1   \n \
	targetForceExponent   $forceExp     \n \
    	upperWalls            $upperWalls     \n \
    	forceConstant         $forceConst1     \n \
	targetEquilSteps      $equilSteps     \n \
    	lambdaSchedule        $lambdaSched      \n \
	targetNumSteps        $stepsPerWindow  \n \
        }"
        
    cv config $TIbias

    return
}
