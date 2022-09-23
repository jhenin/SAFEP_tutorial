import numpy as np
import argparse

def print_lambda_schedule(nstages):
	step = 1./nstages
	lambda_schedule = "1.000"
	for i in range(1,nstages+1,1):
		next_lambda = 1-i*step
		lambda_schedule = lambda_schedule + f" {next_lambda:.3f}"
	print(lambda_schedule)
	return

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("n", type=int, help="The number of stages in your lambda schedule")
	args = parser.parse_args()
	if args.n > 0:
		print_lambda_schedule(args.n)
	else: 
		print("n must be a positive integer")