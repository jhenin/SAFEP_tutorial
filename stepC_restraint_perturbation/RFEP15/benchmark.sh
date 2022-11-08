#!/bin/bash

timeout 300 nohup namd2 +p1 run.namd > p1.log
timeout 300 nohup namd2 +p2 run.namd > p2.log
timeout 300 nohup namd2 +p3 run.namd > p3.log
timeout 300 nohup namd2 +p4 run.namd > p4.log
timeout 300 nohup namd2 +p5 run.namd > p5.log
timeout 300 nohup namd2 +p6 run.namd > p6.log
timeout 300 nohup namd2 +p8 run.namd > p8.log
timeout 300 nohup namd2 +p7 run.namd > p7.log
