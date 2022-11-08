#!/bin/bash

for x in RFEP*; do
	cd $x
	echo $x
	nohup namd2 +p4 run.namd > $x.log &!
	cd ..
done

