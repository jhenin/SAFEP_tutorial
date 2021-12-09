#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Input: fepout file from an interleaved double-wide sampling run from 0 to 1
# (with the last window from 1 to the previous lambda value)
# Output:
# - <base>_fwd.fepout (forward deltaE only)
# - <base>_bwd.fepout (backward deltaE only)
# - <base>_bwd_r.fepout (backward data, with order of windows reversed to simulate a
#   reverse sequential transformation (0 -> 1 -> 0)

# Note:
# The first window samples forward only, the last window backward only
#
# Jérôme Hénin <henin@ibpc.fr> (2018)
# with contributions from Tom Joseph (UPenn)

import sys
import os.path
import argparse
import re
import math
import fileinput

try:
    from scipy.interpolate import InterpolatedUnivariateSpline
    interp_avail = True
except ImportError:
    interp_avail = False


parser = argparse.ArgumentParser(description='Process fepout data from Interleaved Double-Wide Sampling')
parser.add_argument('filenames', nargs='+', help='fepout file name(s) with lambda from 0 to 1')
parser.add_argument('-T', '--temperature', type=float, default=300., help='temperature for FEP estimate, in K')
parser.add_argument('-r', '--reverse', type=bool, default=True, help='reverse order of backward-sampling windows')
parser.add_argument('-i', '--interpolate', action='store_true', help='interpolate data rather than subsampling it')
args = parser.parse_args()

# Store the output in the same directory as the first file
basename = os.path.splitext(args.filenames[0])[0]

interpolating = False
if args.interpolate:
    if interp_avail:
        interpolating = True
    else:
        print('Interpolation requires the scipy.interpolate package. Defaulting to subsampling.')



class Window:
    def __init__(self, l1, l2):
        self.deltaE = []
        self.steps = []
        self.T = []
        self.equilSteps = -1
        self.lambda1 = l1
        self.lambda2 = l2
        self.alchOutFreq = -1

    def interpolate(self):
        equil = True
        acc_exp = 0
        acc_N = 0
        buf = ''
        # Create interpolation function for deltaE data
        func_dE = InterpolatedUnivariateSpline(self.steps, self.deltaE)
        for t in range(len(global_steps)):
            s = global_steps[t]
            # # Skip equilibration steps altogether
            # if s <= self.equilSteps:
            #     continue
            if equil and s > self.equilSteps:
                buf += ('#{} STEPS OF EQUILIBRATION AT LAMBDA {} COMPLETED\n#STARTING COLLECTION OF ENSEMBLE AVERAGE\n').format(self.equilSteps, self.lambda1)
                equil = False
                acc_exp = 0
                acc_N = 0
            dE = float(func_dE(s))
            acc_exp += math.exp(-beta * dE)
            if acc_exp == 0: # catch floating-point underflow due to very high dE (clash)
                acc_exp = 1e-16
            acc_N += 1
            buf += 'FepEnergy: {:6d}  {:14.4f} {:14.4f} {:14.4f} {:14.4f} {:14.4f} {:14.4f} {:14.4f} {:14.4f}\n'.format(s, 0, 0, 0, 0, dE, 0., args.temperature, -RT*math.log(acc_exp / acc_N))
        return buf, -RT*math.log(acc_exp / acc_N)

    def subsample(self):
        equil = True
        acc_exp = 0
        acc_N = 0
        buf = ''

        stride = 1
        # # find first gap in data
        # for i in range(len(self.steps)-1):
        #     diff = self.steps[i+1] - self.steps[i]
        #     if diff != self.alchOutFreq:
        #         if diff % self.alchOutFreq != 0:
        #             sys.exit('Error: gap between data points is ' + str(diff) + ' instead of a multiple of ' + str(alchOutFreq))
        #         else:
        #             stride = (diff - 1) // self.alchOutFreq
        #         break

        # if stride == 0:
        #     stride = 1

        for t in range(0, len(self.steps), stride):
            s = self.steps[t]
            # # Skip equilibration steps altogether
            # if s <= self.equilSteps:
            #     continue
            if equil and s > self.equilSteps:
                buf += ('#{} STEPS OF EQUILIBRATION AT LAMBDA {} COMPLETED\n#STARTING COLLECTION OF ENSEMBLE AVERAGE\n').format(self.equilSteps, self.lambda1)
                equil = False
                acc_exp = 0
                acc_N = 0
            dE = self.deltaE[t]
            acc_exp += math.exp(-beta * dE)
            if acc_exp == 0: # catch floating-point underflow due to very high dE (clash)
                acc_exp = 1e-16
            acc_N += 1
            buf += 'FepEnergy: {:6d}  {:14.4f} {:14.4f} {:14.4f} {:14.4f} {:14.4f} {:14.4f} {:14.4f} {:14.4f}\n'.format(t, 0, 0, 0, 0, dE, 0., self.T[t], -RT*math.log(acc_exp / acc_N))
        return buf, -RT*math.log(acc_exp / acc_N)


    if interpolating:
        resample = interpolate
    else:
        resample = subsample


    def process(self, outfile, is_backward, dG_in):
        if len(self.steps) == 0:
            return
        if is_backward:
            print("Backward data: " + str(self.lambda1) + " to " + str(self.lambda2) + " (" + str(len(self.steps)) + " points)")
        else:
            print("Forward data: " + str(self.lambda1) + " to " + str(self.lambda2) + " (" + str(len(self.steps)) + " points)")

        outfile.write("#NEW FEP WINDOW: LAMBDA SET TO {} LAMBDA2 {}\n".format(self.lambda1, self.lambda2))
        buf, dG = self.resample()
        dG_out = dG_in + dG
        outfile.write(buf)
        outfile.write('#Free energy change for lambda window [ {} {} ] is {} ; net change until now is {}\n'.format(self.lambda1, self.lambda2, dG, dG_out))
        return dG_out


# Lists of forward and backward windows data
windows_f = []
windows_b = []
IDWS_on = False
# Steps for all data, for interpolation
global_steps = []

first_window = True
last_window = False
RT = 1.98720650096e-3 * args.temperature #in kcal/mol
beta = 1. / RT

# Global step counter, because we don't depend on the fepout files to have consistent step numbers
step_counter = -1
alchOutFreq = -1
equilSteps = 0
equil_re = re.compile(r'#(\d*) STEPS OF EQUILIBRATION')

for line in fileinput.input(args.filenames):
    fields = re.split(' +', line.strip())

    # Most lines in the file are energy values, so start with them
    # Collect all timestep numbers (both fwd and bwd) for interpolation
    if line.startswith('FepE'):
        l = len(fields)
        if l == 8:
            i_temp = 7
        elif l == 10:
            i_temp = 8
        else:
            print("Wrong number of fields (expected 8 or 10) in line:\n" + line)

        # Detect first step and value of alchOutFreq
        if step_counter == -1:
            first_step = int(fields[1])
            step_counter = first_step
        elif step_counter == first_step:
            wf.alchOutFreq = int(fields[1]) - first_step
            if IDWS_on:
                wb.alchOutFreq = wf.alchOutFreq
            step_counter += wf.alchOutFreq

        if first_window:
            global_steps.append(step_counter)

        # Discard every other sample if not doing IDWS or interpolating
        if line.startswith('FepEnergy:') and (interpolating or IDWS_on or (step_counter//alchOutFreq)%2 == 0):
            wf.steps.append(step_counter)
            wf.deltaE.append(float(fields[6]))
            wf.T.append(float(fields[8]))

        if line.startswith('FepE_back:'):
            wb.steps.append(step_counter)
            wb.deltaE.append(float(fields[6]))
            wb.T.append(float(fields[8]))

        if wf.alchOutFreq > 0:
            step_counter += wf.alchOutFreq

    elif line.startswith('#NEW FEP WINDOW:'):

        if step_counter >= 0:
            first_window = False
        step_counter = -1

        # Extract lambda values
        cur_lambda1 = float(fields[6])
        cur_lambda2 = float(fields[8])

        # create new "forward" window
        wf = Window(cur_lambda1, cur_lambda2)

        if len(fields) == 11:
            cur_lambdaIDWS = float(fields[10])
            IDWS_on = True
            wb = Window(cur_lambda1, cur_lambdaIDWS)
        else:
            cur_lambdaIDWS = -1.
            IDWS_on = False

        # Check that we start from an end-point
        if first_window:
            if cur_lambda1 == 1:
                start_from_one = True
            elif cur_lambda1 == 0:
                start_from_one = False
            else:
                sys.exit('Lambda should start at zero or one, found ' + str(cur_lambda1))

        # Skip the last window for the "forward" output, as it is backward from 1
        if ((not start_from_one and cur_lambda1 == 1) or (start_from_one and cur_lambda1 == 0)):
            last_window = True
            # append the "forward" data to the list of backward windows
            windows_b.append(wf)
        else:
            windows_f.append(wf)

        if IDWS_on:
            windows_b.append(wb)

    match = equil_re.match(line)
    if match != None:
        equilSteps = int(match.group(1))
        wf.equilSteps = equilSteps
        if IDWS_on:
            wb.equilSteps = equilSteps


outf = open(basename + '_fwd.fepout', 'w')
dG_acc = 0.
for w in windows_f:
    dG_acc = w.process(outf, False, dG_acc)
outf.close()

if args.reverse:
    windows_b.reverse()

outb = open(basename + '_bwd.fepout', 'w')
dG_acc = 0.
for w in windows_b:
    dG_acc = w.process(outb, True, dG_acc)
outb.close()
