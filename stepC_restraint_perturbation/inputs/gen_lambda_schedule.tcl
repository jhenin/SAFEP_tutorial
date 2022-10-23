proc genLambdas { num } {
  for {set k 0.0} {$k <= $num} {set k [expr $k+1.0]} {
    puts -nonewline "[expr $k/$num] "
  }
  puts " "
}
