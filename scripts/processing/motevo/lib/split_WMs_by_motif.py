#!/usr/bin/env python
'''
DESCRIPTION

    Downloaded motifs are in a single textfile, split them up in a script

FOR HELP

    python split_WMs_by_motif.py --help

AUTHOR:      Jake Yeung (jake.yeung@epfl.ch)
LAB:         Computational Systems Biology Lab (http://naef-lab.epfl.ch)
CREATED ON:  2017-06-14
LAST CHANGE: see git log
LICENSE:     MIT License (see: http://opensource.org/licenses/MIT)
'''

import sys, argparse, datetime

def write_pwm_to_file(line, outdir, readfile):
    '''
    Write PWMs to file, make it look like Saeed's PWMs
    line: first line of PWM file
    readfile: contains next several lines of PWM files: stop when you hit //
    '''
    motif=line.split("  ")[1].strip()
    outf='/'.join([outdir, ''.join([motif, ".pwm"])])
    with open(outf, "wb") as outf:
        outf.write("//\n")  # init
        outf.write(line)  # write first line
        while not line.startswith("//"):
            line = readfile.next()
            outf.write(line)  # finishes with //
    return None


def main():
    parser = argparse.ArgumentParser(description='Downloaded motifs are in a single textfile, split them up in a script')
    parser.add_argument('infile', metavar='INFILE',
                        help='mm10_weight_matrices')
    parser.add_argument('outdir', metavar='OUTDIR',
                        help='outdir to write individual PWMs')
    parser.add_argument('--quiet', '-q', action='store_true',
                        help='Suppress some print statements')
    args = parser.parse_args()

    # store command line arguments for reproducibility
    CMD_INPUTS = ' '.join(['python'] + sys.argv)    # easy printing later
    # store argparse inputs for reproducibility / debugging purposes
    args_dic = vars(args)
    ARG_INPUTS = ['%s=%s' % (key, val) for key, val in args_dic.iteritems()]
    ARG_INPUTS = ' '.join(ARG_INPUTS)

    # Print arguments supplied by user
    if not args.quiet:
        print(datetime.datetime.now().strftime('Code output on %c'))
        print('Command line inputs:')
        print(CMD_INPUTS)
        print ('Argparse variables:')
        print(ARG_INPUTS)

    with open(args.infile, "rb") as inf:
        for line in inf:
            # look for motif, then write PWM to file
            if line.startswith("NA"):
                write_pwm_to_file(line, args.outdir, inf)


if __name__ == '__main__':
    main()
