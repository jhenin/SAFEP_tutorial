proc genLambdas { num } {
    set lambdaSched ""
    for {set k 0} {$k <= $num} {set k [expr $k+1.0]} {
      append lambdaSched "[expr $k/$num] "
    }
    return $lambdaSched
}

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
}
