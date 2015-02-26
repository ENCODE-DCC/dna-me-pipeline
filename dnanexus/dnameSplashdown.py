#!/usr/bin/env python
import argparse
import os
import sys
import itertools

import dxpy
from dxencode import dxencode as dxencode
import json

# NOTES: This command-line utility will run the DNA methylation for a single replicate
#      - All results will be written to a folder /RESULT_FOLDER_DEFAULT/<expId>/<biorep>_<techrep>.
#      - If any results are already in that directory, then the steps that created those results
#        will not be rerun.
#      - If any jobs for the experiment and replicate are already running, nothing new will be
#        launched.
#      - Most (BUT NOT ALL) of the code is generic, relying on the hard-coded JSON below and
#        relying on hard-coded tokens for step, file and parameter names:
#        - STEP_ORDER is the list of steps for single or paired-end reads
#        Long RNA-seq specific code is marked with '### LRNA specific'

GENOME_DEFAULT = 'hg19'
''' This the default Genome that human experiments are mapped to.'''

PROJECT_DEFAULT = 'dna-me-pipeline'
''' This the default DNA Nexus project to use for the long RNA-seq pipeline.'''

REF_PROJECT_DEFAULT = 'dna-me-pipeline'  ## 'ENCODE Reference Project'
''' This the default DNA Nexus project to find reference files in.'''

REF_FOLDER_DEFAULT = '/ref'
''' This the default folder that reference files are found in.'''

#RESULT_FOLDER_DEFAULT = '/runs'  Not sure which is better.
RESULT_FOLDER_DEFAULT = '/runs'

''' This the default location to place results folders for each experiment.'''

RUNS_LAUNCHED_FILE = "launchedRuns.txt"

STEP_ORDER = {
    # for SE or PE the list in order of steps to run
    'se': [ 'trim-se', 'map-se', 'extract' ],
    'pe': [ 'trim-pe', 'map-pe', 'extract' ]
    }

GENOME_REFERENCES = {
    # For looking up reference file names.
    # TODO: should use ACCESSION based fileNames
    'meIndex':  {
                    'hg19': {
                            'female':  "hg19.female.bisulfite.tgz",
                            'male': "hg19.female.bisulfite.tgz"
                            },
                    'mm10': {
                            'female': "female.mm10.chrom.fa.gz.bisulfite.tgz",
                            'male': "male.mm10.chrom.fa.gz.bisulfite.tgz"
                    },
                    'lambda': {
                            'female': "chrL.fa.gz.bisulfite.tgz",
                            'male': "chrL.fa.gz.bisulfite.tgz"
                    }
    },
    'genome':  {
                    'hg19': {
                            'female':  "female.hg19.fa.gz",
                            'male': "male.hg19.fa.gz"
                            },
                    'mm10': {
                            'female': "female.mm10.chrom.fa.gz",
                            'male': "male.mm10.chrom.fa.gz"
                    },
                    'lambda': {
                            'female': "chrL.fa.gz",
                            'male': "chrL.fa.gz"
                    }

    },
    'chrom_sizes':   {
                    'hg19': {
                            'female':   'female.hg19.chrom.sizes',
                            'male':     'male.hg19.chrom.sizes'
                            },
                    'mm10': {
                            'female':   'female.mm10.chrom.sizes',
                            'male':     'male.mm10.chrom.sizes'
                            },
                    'lambda': {
                            'female': "chrL.chrom.sizes",
                            'male': "chrL.chrom.sizes"
                    }
    }
}
GENOME_MAPPING = {
    "human": "hg19",
    "mouse": "mm10"

}

POST_TEMPLATES = {
    "mapped_reads": {
        "file_format": "bam",
        "output_type": "alignments",
        "derived_from": ["all_reads"]
    },
    "CG": {
        "file_format": "bed_bedMethyl",
        "output_type": "methyl CG",
        "derived_from": ["mapped_reads"]
    },
    "CHG": {
        "file_format": "bed_bedMethyl",
        "output_type": "methyl CHG",
        "derived_from": ["mapped_reads"]
    },
    "CHH": {
        "file_format": "bed_bedMethyl",
        "output_type": "methyl CHH",
        "derived_from": ["mapped_reads"]
    },
    "CGbb": {
        "file_format": "bedMethyl",
        "output_type": "methyl CG",
        "derived_from": ["CG"]
    },
    "CHGbb": {
        "file_format": "bedMethyl",
        "output_type": "methyl CHG",
       "derived_from": ["CHG"]
    },
    "CHHbb": {
        "file_format": "bedMethyl",
        "output_type": "methyl CHH",
        "derived_from": ["CHH"]
    }
}

def get_args():
    '''Parse the input arguments.'''
    ### LRNA specific
    ap = argparse.ArgumentParser(description="Launches long RNA-seq pipeline analysis for " +
                "one replicate on single or paired-end reads. Can be run repeatedly and will " +
                "launch only the steps that are needed to finished the pipeline. All results " +
                "will be placed in the folder /<resultsLoc>/<experiment>/<replicate>.")
    ### LRNA specific

    ap.add_argument('-e', '--experiment',
                    help='ENCODED experiment accession',
                    required=True)

    ap.add_argument('--br', '--biological-replicate',
                    help="Biological Replicate number (default: 1)",
                    default=1,
                    type=int,
                    required=True)

    ap.add_argument('--tr', '--technical-replicate',
                    help="Biological Replicate number (default: 1)",
                    default=1,
                    type=int,
                    required=True)

    ap.add_argument('--project',
                    help="Project to run analysis in (default: '" + PROJECT_DEFAULT + "')",
                    default=PROJECT_DEFAULT,
                    required=False)

    ap.add_argument('--refLoc',
                    help="The location to find reference files (default: '" + \
                                            REF_PROJECT_DEFAULT + ":" + REF_FOLDER_DEFAULT + "')",
                    default=REF_FOLDER_DEFAULT,
                    required=False)

    ap.add_argument('--resultsLoc',
                    help="The location to to place results folders (default: '<project>:" + \
                                                                    RESULT_FOLDER_DEFAULT + "')",
                    default=RESULT_FOLDER_DEFAULT,
                    required=False)

    ap.add_argument('--testserver',
                    help="Use the test server designated in keypairs.json",
                    action='store_true',
                    required=False)

    ap.add_argument('--test',
                    help='Test run only, do not launch anything.',
                    action='store_true',
                    required=False)

    ap.add_argument('--skipvalidate',
                    help='Skip running Validate Files',
                    action='store_true',
                    required=False)

    return ap.parse_args()

STEPS = {}  ## back backwards compatibility

def calculate_steps(applet):
    ''' create input object for a step '''
    try:
        dxapp = json.load(open(applet+'/dxapp.json'))
    except IOError:
        print "Cannot open applet json for %s" % applet
        sys.exit(0)

    params = {}
    inputs = {}
    results = {}
    inps = dxapp.get('inputSpec') or []
    for inp in inps:
        if inp['class'] == 'file' or inp['class'] == 'array:file':
            inputs[inp['name']] = inp['name']
        else:
            params[inp['name']] = inp['name']

    outs = dxapp.get('outputSpec') or []
    for out in outs:
        if out['class'] == 'file' or out['class'] == 'array:file':
            results[out['name']] = '/'+out['patterns'][0]
        else:
            pass
            # TODO not sure what to do with these

    return {
        'app': applet,
        'params': params,
        'inputs': inputs,
        'results': results
    }


def pipeline_specific_vars(args, mapping, pairedEnd, gzip=False):
    '''Adds pipeline specific variables to a dict, for use building the workflow.'''
    extras = {}
    extras['genome']   = mapping['genome']
    extras['gender']     = mapping['sex']
    extras['experiment'] = args.experiment
    extras['replicate']  = mapping['replicate']
    extras['library']    = mapping['library']
    extras['pairedEnd']  = pairedEnd
    extras['gzip'] = gzip
    if pairedEnd:
        extras['max_insert'] = 500
        extras['min_insert'] = 0
        # these might have other production values
    # workflow labeling
    genderToken = "XY"
    if extras['gender'] == 'female':
        genderToken = "XX"
    extras['description'] = "The ENCODE DNA Methylation Pipeline for WGBS"
    extras['title'] = "DNA Methylation single-end "
    extras['name'] = "dname_"+extras['genome']+genderToken+"SE_"
    if pairedEnd:
        extras['title'] = "DNA Methylation paired-end "
        extras['name'] = "dname_"+extras['genome']+genderToken+"PE_"
    extras['title'] += args.experiment+" - rep"+extras['replicate'] + " (library '"+extras['library']+"')"
    extras['name']  += args.experiment+"_rep"+extras['replicate']
    extras['subTitle'] = extras['genome']+", "+extras['gender']+"'."

    args.resultsLoc = RESULT_FOLDER_DEFAULT  # not sure we need genome
    extras['resultsFolder'] = args.resultsLoc + '/' + args.experiment + '/' + mapping['replicate']


    return extras

def get_software():

    ## NOTE bedToBigBed 2.6 version must be added!!!
    ## This data structure is wrong as well, should be "software_versions":[ { version: v, software: sw } ]
    return {
        "software_versions": [
            {"bismark": "v0.10.0"},
            {"samtools": "0.1.19-44428cd"},
            {"cxrepo-bed.py": "0.2"}
        ]

    }
#######################
def main():
    args = get_args()

    (AUTHID,AUTHPW,SERVER) = dxencode.processkey('default')
    url = SERVER + 'experiments/%s/?format=json&frame=embedded' %(args.experiment)
    response = dxencode.encoded_get(url, AUTHID, AUTHPW)
    exp = response.json()

    if not exp.get('replicates') or len(exp['replicates']) < 1:
        print "No replicates found in %s\n%s" % ( args.experiment, exp )
        sys.exit(1)

    #replicate = "rep%s_%s" % (args.br, args.tr)
    replicate = "%s_%s" % (args.br, args.tr)

    reps_mapping = dxencode.choose_mapping_for_experiment(exp)
    # could try to do all replicates here
    try:
        mapping = reps_mapping[(args.br,args.tr)]
    except KeyError:
        print "Specified replicate: %s could not be found in mapping." % replicate
        print reps_mapping
        sys.exit(1)

    mapping['replicate'] = replicate

    try:
        mapping['genome'] = GENOME_MAPPING[mapping.get('organism', "Not Found")]

    except KeyError:
        print "Organism %s not currently supported" % mapping['organism']
        sys.exit(1)

    if mapping['unpaired'] and not mapping['paired']:
        pairedEnd = False
    elif mapping['paired'] and not mapping['unpaired']:
        pairedEnd = True
    elif not mapping['unpaired'] and not mapping['paired']:
        print "Replicate has no reads either paired or unpaired"
        print mapping
        sys.exit(1)
    else:
        print "Replicate has both paired(%s) and unpaired(%s) reads, quitting." % (len(mapping['paired'], len(mapping['unpaired'])))
        print mapping
        sys.exit(1)

    psv = pipeline_specific_vars(args, mapping, pairedEnd)
    project = dxencode.get_project(args.project)
    projectId = project.get_id()


    ## TODO this is a bunch of ugly
    if pairedEnd:
        paired_fqs = {
            '1': [],
            '2': []
        }
        read1s = []
        read2s = []
        for (p1, p2) in mapping['paired']:
            paired_fqs[p1['paired_end']].append(p1['accession']+".fastq.gz")
            paired_fqs[p2['paired_end']].append(p2['accession']+".fastq.gz")
            read1s.append(p1['accession'])
            read2s.append(p2['accession'])
        pipePath = STEP_ORDER['pe']
        print "Generating workflow steps (paired-end)..."
    else:
        unpaired_fqs = [ f['accession']+".fastq.gz" for f in mapping['unpaired'] ]
        pipePath = STEP_ORDER['se']

    for step in pipePath:
        STEPS[step] = calculate_steps(step)

    pipeSteps = STEPS
    ## warning ugly kludge here
    file_globs = {}
    for app in STEPS.keys():
        for token in STEPS[app]['results'].keys():
            file_globs[token] = STEPS[app]['results'][token]

    print "Checking for prior results..."

    priors = dxencode.find_prior_results(pipePath,pipeSteps,psv['resultsFolder'],file_globs, projectId)

    if pairedEnd:
        priors['pair1_reads'] = dxencode.find_file_set(paired_fqs["1"], projectId)
        priors['pair2_reads'] = dxencode.find_file_set(paired_fqs["2"], projectId)
        priors['all_reads'] = priors['pair1_reads'] + priors['pair2_reads']
        submitted = {
            'all_reads': read1s + read2s
        }
    else:
        priors['reads'] = dxencode.find_file_set(unpaired_fqs, projectId)
        priors['all_reads'] = priors['reads']
        submitted = {
            'all_reads': [ f['accession'] for f in mapping['unpaired']],
        }


    print "Determining steps to run..."
    #print priors
    #sys.exit(1)
    # NOTE: stepsToDo is an ordered list of steps that need to be run
    deprecateFiles = [] # old results will need to be moved/removed if step is rerun
    stepsToDo = dxencode.determine_steps_to_run(pipePath,pipeSteps, priors, deprecateFiles, projectId, verbose=True)

    print "Checking for currently running analyses..."
    dxencode.check_run_log(psv['resultsFolder'],projectId, verbose=True)

    if len(stepsToDo):
        print "Pipeline incomplete, please resubmit jobs: %s" % stepsToDo
        sys.exit(0)

    print priors
    to_submit = [ k for k in priors.keys() if POST_TEMPLATES.get(k) ]
    n = 0 # skip reads
    print "Attempting to submit %s files to args.experiment" % len(to_submit)
    while(to_submit):
        if n > len(priors) * len(priors):
            print "Too many itereations: %s" % priors
            break
        token = to_submit.pop(0)
        print "%s %s - %s" % (token, priors[token], n)
        f_ob = POST_TEMPLATES.get(token, None)
        n += 1
        if f_ob:
            derive_check = f_ob.get('derived_from', [])
            if derive_check:
                derived = [ submitted[f] for f in derive_check if submitted.get(f) ]
                if not derived:
                    to_submit.append(token)
                    continue
                else:
                    f_ob['derived_from'] = list(itertools.chain(*derived))
            dxFile = dxpy.DXFile(dxid=priors[token])
            print "Post File: %s %s" % (token, dxFile.name)
            f_ob['dataset'] = args.experiment
            f_ob['lab'] = '/labs/j-michael-cherry/'
            f_ob['award'] = '/awards/U41HG006992/'
            f_ob['assembly'] = mapping['genome']
            ## temporary haxors until file display works
            f_ob['replicate'] = mapping['replicate_id']
            f_ob['notes'] = json.dumps(dxencode.create_notes(dxFile, get_software()))
            print json.dumps(f_ob, sort_keys=True, indent=4, separators=(',',': '))
            if args.testserver:
                server = 'test'
            else:
                server = 'www'

            if args.test:
                fake_acc = 'ENCFF%03dAAA' % n
                print "Fake submission: %s" % fake_acc
                submitted[token] = [ fake_acc ]
            else:
                applet = dxencode.find_applet_by_name('validate-post', projectId )
                job = applet.run({
                    "pipe_file": dxpy.dxlink(dxFile),
                    "file_meta": f_ob,
                    "key": server,
                    "debug": True,
                    "skipvalidate": args.skipvalidate or False
                    })
                print "Submitting %s" % job.id
                job.wait_on_done(interval=1)
                accession = job.describe()['output'].get('accession', "Unknown Acc")
                error = job.describe()['output'].get('error', "Unknown Error")
                submitted[token] = [ accession ]
                print "Posted (%s): %s" % (error, accession)

    # Exit if test only
    if args.test:
        print "Fake submitted %s files." % n
    if args.test:
        sys.exit(0)


if __name__ == '__main__':
    main()


