#!/usr/bin/env python
'''
DESCRIPTION

    Make trackhubs for UCSC browser using bigbed files.
    Tutorials for making trackhubs with trackhub python package from:
    https://pythonhosted.org/trackhub/tutorial.html#creating-a-composite-track

    Motevo Motif sitecounts: bigbed files

    rewrite in python 3 2019-02-01

FOR HELP

    python make_trackhubs.py --help

AUTHOR:      Jake Yeung (jake.yeung@epfl.ch)
LAB:         Computational Systems Biology Lab (http://naef-lab.epfl.ch)
CREATED ON:  2016-01-06
LAST CHANGE: see git log
LICENSE:     MIT License (see: http://opensource.org/licenses/MIT)

EXAMPLE

'''

import warnings
import sys, argparse, os
from trackhub import Hub, GenomesFile, Genome, TrackDb, Track
from trackhub import CompositeTrack
from trackhub.track import SubGroupDefinition
from trackhub import ViewTrack

from trackhub.upload import stage_hub, upload_hub

# import numpy as np

def get_files_from_dir(inputdir, ext=".bb"):
    '''
    Get *.bb files
    Output as dic with {name: name.bb}
    '''
    fnames = os.listdir(inputdir)
    fnames_dic = {}
    for f in fnames:
        if f.endswith(ext):
            # sample.bw -> sample
            sample = f.split(".")[0]
            # handle sample if contains { or }
            for badchar in ["{", "}", ","]:
                sample = sample.replace(badchar, "")
            if sample not in fnames_dic:
                fnames_dic[sample] = []
            # make new fname
            newf = ''.join([sample, ext])
            fnames_dic[sample].append(os.path.join(inputdir, newf))
    return(fnames_dic)

def main():
    parser = argparse.ArgumentParser(description='Make trackhubs for UCSC browser using bigBed files. \
                                     Outputs to CURRENT DIRECTORY.')
    parser.add_argument('inputdir', metavar='INDIR',
                        help='Directory containing bigBed files .bb ending')
    parser.add_argument('outdir', metavar='OUTDIR',
                        help='Directory for staging files')
    parser.add_argument('--quiet', '-q', action='store_true',
                         help='Suppress some print statements')
    parser.add_argument('--render', '-r', action='store_true',
                         help='Render file to current dir')
    parser.add_argument('--upload', '-u', action='store_true',
                         help='Upload file to webserver')
    parser.add_argument('--mm9', '-m', action='store_true',
                         help='Switch from mm10 to mm9')
    parser.add_argument('--has_strand', '-s', action='store_true',
                         help='Bed has strand (changes from 5 columns to 6)')
    parser.add_argument('--suffix', '-S', metavar="trackhub label suffix", default = "",
                         help='Suffix to label, for example H3K4me1')
    args = parser.parse_args()

    # store command line arguments for reproducibility
    CMD_INPUTS = ' '.join(['python'] + sys.argv)    # easy printing later
    # store argparse inputs for reproducibility / debugging purposes
    args_dic = vars(args)
    ARG_INPUTS = ['%s=%s' % (key, val) for key, val in args_dic.iteritems()]
    # ARG_INPUTS = ['%s=%s' % (key, val) for key, val in args_dic.items()]
    ARG_INPUTS = ' '.join(ARG_INPUTS)

    # Print arguments supplied by user
    if not args.quiet:
        print('Command line inputs:')
        print(CMD_INPUTS)
        print ('Argparse variables:')
        print(ARG_INPUTS)

    # define constants (hard coded)
    if args.mm9:
        genobuild = "mm9"
    else:
        genobuild = "mm10"
    jsuffix = "%s_%s" % (genobuild, args.suffix)
    print("Assigning prefix: %s" %jsuffix)
    # dirname: motevo_from_peaks/H3K4me1_peaks
    dirname = "motevo_from_peaks/%s_peaks/motevo_motifs_%s" % (args.suffix, jsuffix)
    hubname = "motevo_motifs_%s" % jsuffix
    shortlab = "motevo_%s" % jsuffix
    longlab = "Motevo motifs %s" % jsuffix
    email = "jake.yeung@epfl.ch"
    # url = "http://upnaepc2.epfl.ch"
    url = "http://upnaesrv1.epfl.ch"
    assay = "bigbed"
    jvis = "dense"
    # bigbed options loaded into ViewTrack
    jspectrum = "on"
    scoremax = 1000
    scoremin = 500

    # define URLs
    url_main = "%s/%s" % (url, dirname)
    url_base = "%s/%s/data" % (url, dirname)
    # upload_main = "~/Sites/%s" % dirname
    # upload_base = "~/Sites/%s/data" % dirname
    upload_main = "%s" % hubname
    upload_base = "%s/data" % hubname
    if not args.has_strand:
        ftype = "bigBed 5"
    else:
        ftype = "bigBed 6"
    # host = "circadian.epfl.ch"
    # user = "web"
    host = "upnaesrv1.epfl.ch"
    user = "websync"

    # define constants
    genomebuild = genobuild

    files_dic = get_files_from_dir(args.inputdir, ext=".bb")

    samples_dic = {}
    for sample in files_dic.keys():
        samples_dic[sample] = sample

    # init hub genomes file genome trackdb
    # Make my hub
    hub = Hub(hub = hubname,
              short_label = shortlab,
              long_label = longlab,
              email = email)
              # url = "%s/%s" % (url, dirname))

    hub.url = os.path.join(url_main, "%s.hub.txt" % hub.hub)

    genomes_file = GenomesFile()
    genome = Genome(genomebuild)
    trackdb = TrackDb()

    # add remote fn
    # hub.remote_fn = os.path.join(upload_main, "hub.txt")
    # genomes_file.remote_fn = os.path.join(upload_main, "genomes.txt")
    hub.remote_fn = upload_main
    genomes_file.remote_fn = upload_main
    trackdb.remote_fn = os.path.join(upload_main, genomebuild, "trackDb.txt")

    hub.add_genomes_file(genomes_file)
    genome.add_trackdb(trackdb)
    genomes_file.add_genome(genome)

    # init composite
    composite = CompositeTrack(name = hubname,
                               short_label = shortlab,
                               long_label = longlab,
                               tracktype = ftype)
    # make subgroups
    subgroups = [

        SubGroupDefinition(name = "sample",
                           label = "sample",
                           mapping = samples_dic),
    ]
    composite.add_subgroups(subgroups)
    # make viewTrack, a hierarchy containing my files, for example
    view = ViewTrack(name = "%sViewTrack" % assay,
                     view = "%s" % assay,
                     visibility = jvis,
                     tracktype = ftype,
                     short_label = "%s" % assay,
                     long_label = "%s assay" % assay,
                     # big bed labels
                     spectrum = jspectrum,
                     scoreMin = scoremin,
                     scoreMax = scoremax)
    composite.add_view(view)

    # make track
    for sample, wfs in files_dic.iteritems():
        for wf in wfs:
            sampname = os.path.basename(wf)
            bname = sampname.split(".")[0]
            track = Track(name = bname,
                            tracktype = ftype,
                            url = os.path.join(url_base, "%s" % sampname),
                            local_fn = os.path.abspath(wf),
                            remote_fn = os.path.join(upload_base, "%s" % sampname),
                            visibility = jvis,
                            shortLabel = bname,
                            longLabel = bname,
                            spectrum = jspectrum,
                            scoreMin = scoremin,
                            scoreMax = scoremax,
                            subgroups = {"sample" : sample})
            view.add_tracks(track)
    trackdb.add_tracks(composite)

    print('Track looks like this:')
    print(trackdb)

    if args.render:
        # print('Rendering to %s' % hub.local_fn)
        # results = hub.render()
        # upload_hub(hub=hub, host='localhost', remote_dir='example_grouping_hub')
        stage_hub(hub, staging=args.outdir)

    if args.upload:
        print('Uploading to web@circadian.epfl.ch')
        # for track in trackdb.tracks:
        #     upload_track(track = track, host = host, user = user)
        upload_hub(hub = hub, host = host, user = user, remote_dir="/data/web/sites/motevo_from_peaks")

    print('Subgroups:')
    for sg in subgroups:
        print(sg)
    print("Staging to path: %s" % args.outdir)
    # print('Hub URL if you uploaded to server: %s' % hub.url)

if __name__ == '__main__':
    main()
