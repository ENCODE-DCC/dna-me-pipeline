#!/usr/bin/env python
import argparse
import os
import sys
import subprocess
from datetime import datetime

import dxpy
import dxencode
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

PROJECT_DEFAULT = 'dna-me-pipline'
''' This the default DNA Nexus project to use for the long RNA-seq pipeline.'''

REF_PROJECT_DEFAULT = 'dna-me-pipeline'
''' This the default DNA Nexus project to find reference files in.'''

REF_FOLDER_DEFAULT = '/ref'
''' This the default folder that reference files are found in.'''

#RESULT_FOLDER_DEFAULT = '/runs'  Not sure which is better.
RESULT_FOLDER_DEFAULT = '/'

''' This the default location to place results folders for each experiment.'''

RUNS_LAUNCHED_FILE = "launchedRuns.txt"

STEP_ORDER = {
    # for SE or PE the list in order of steps to run
    'se': [ 'trim-se', 'map-se', 'extract' ],
    'pe': [ 'trim-pe', 'map-pe', 'extract' ]
    }

GENOME_REFERENCES = {
    # For looking up reference file names.
    # TODO: should remove annotation if only one per genome
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
    'chromSizes':   {
                    'hg19': {
                            'female':   'female_hg19.chrom.sizes',
                            'male':     'male_hg19.chrom.sizes'
                            },
                    'mm10': {
                            'female':   'female.mm10.chrom.sizes',
                            'male':     'male.mm10.chrom.sizes'
                            }
                    }
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

    ap.add_argument('-br', '--biological-replicate',
                    help="Biological Replicate number (default: 1)",
                    default='1',
                    required=True)

    ap.add_argument('-tr', '--technical-replicate',
                    help="Biological Replicate number (default: 1)",
                    default='1',
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

    ap.add_argument('--test',
                    help='Test run only, do not launch anything.',
                    action='store_true',
                    required=False)

    ap.add_argument('--force',
                    help='Force rerunning all steps.',
                    action='store_true',
                    required=False)

    return ap.parse_args()

STEPS = {}  ## back backwards compatibility

def calculate_steps(applet):
    ''' create input object for a step '''
    try:
        dxapp = json.load(open(applet+'dxapp.json'))
    except IOError:
        print "Cannot open applet json for %s" % applet
        sys.exit(0)

    params = {}
    inputs = {}
    results = {}
    inps = dxapp.get('inputSpec') or []
    for inp in inps:
        if inp['class'] == 'file':
            inputs[inp['name']] = inp['name']
        else:
            params[inp['name']] = inp['name']

    outs = dxapp.get('outputSpec') or []
    for out in outs:
        if inp['class'] == 'file':
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


def pipelineSpecificExtras(genome, gender, experiment, replicate, library, pairedEnd):
    '''Adds pipeline specific variables to a dict, for use building the workflow.'''
    extras = {}
    extras['genome']   = genome
    extras['gender']     = gender
    extras['experiment'] = experiment
    extras['replicate']  = replicate
    extras['library']    = library
    extras['pairedEnd']  = pairedEnd
    # workflow labeling
    genderToken = "XY"
    if gender == 'female':
        genderToken = "XX"
    extras['description'] = "The ENCODE DNA Methylation Pipeline for WGBS"
    extras['title'] = "DNA Methylation single-end "
    extras['name'] = "dname_"+genome+genderToken+"SE_"
    if pairedEnd:
        extras['title'] = "DNA Methylation paired-end "
        extras['name'] = "dname_"+genome+genderToken+"PE_"
    extras['title'] += experiment+" - rep"+replicate + " (library '"+library+"')"
    extras['name']  += experiment+"_rep"+replicate
    extras['subTitle'] = genome+", "+gender+"'."

    return extras


def findPriorResults(pairedEnd,resultsFolder,projectId):
    '''Looks for all result files in the results folder.'''
    priors = {}
    steps = []
    if pairedEnd:
        steps = STEP_ORDER['pe']
    else:
        steps = STEP_ORDER['se']
    for step in steps:
        for fileToken in STEPS[step]['results'].keys():
            fid = dxencode.find_file(resultsFolder + STEPS[step]['results'][fileToken],projectId)
            if fid != None:
                priors[fileToken] = fid
    return priors

def findReferenceFiles(priors,refLoc,extras):
    '''Locates all reference files based upon gender, genome and annotation.'''
    ### LRNA specific
    refLoc=refLoc+'/'+extras['genome']
    meIx = refLoc+'/'+GENOME_REFERENCES['meIndex'][extras['genome']][extras['gender']]
    meIxFid = dxencode.find_file(meIx,REF_PROJECT_DEFAULT)
    if meIxFid == None:
        sys.exit("ERROR: Unable to locate DNA Methylation index file '" + meIx + "'")
    else:
        priors['meIndex'] = meIxFid

    chromSizes = refLoc+'/'+GENOME_REFERENCES['chromSizes'][extras['genome']][extras['gender']]
    chromSizesFid = dxencode.find_file(chromSizes,REF_PROJECT_DEFAULT)
    if chromSizesFid == None:
        sys.exit("ERROR: Unable to locate Chrom Sizes file '" + chromSizes + "'")
    else:
        priors['chromSizes'] = chromSizesFid

def determineStepsToDo(pairedEnd, priors, deprecate, projectId, force=False):
    '''Determine what steps need to be done, base upon prior results.'''
    willCreate = []
    stepsToDo = []
    steps = []
    if pairedEnd:
        steps = STEP_ORDER['pe']
    else:
        steps = STEP_ORDER['se']
    for step in steps:
        # Force will include the first step with all its inputs
        # This should avoid forcing concat if it isn't needed
        #
        if force:
            inputs = STEPS[step]['inputs'].keys()
            count = 0
            for input in inputs:
                if input in priors:
                    count += 1
            if count == len(inputs):
                stepsToDo += [ step ]
        if step not in stepsToDo:
            results = STEPS[step]['results'].keys()
            for result in results:
                if result not in priors:
                    #print "- Adding step '"+step+"' because prior '"+result+"' was not found."
                    stepsToDo += [ step ]
                    break
        # If results are there but inputs are being recreated, then step must be rerun
        if step not in stepsToDo:
            inputs = STEPS[step]['inputs'].keys()
            for inp in inputs:
                if inp in willCreate:
                    #print "- Adding step '"+step+"' due to prior step dependency."
                    stepsToDo += [ step ]
                    break
        # Any step that is rerun, will cause prior results to be deprecated
        # NOTE: It is necessary to remove from 'priors' so succeeding steps are rerun
        # NOTE: It is also important to move prior results out of target folder to avoid confusion!
        if step in stepsToDo:
            results = STEPS[step]['results'].keys()
            for result in results:
                willCreate += [ result ]
                if result in priors:
                    deprecate += [ priors[result] ]
                    del priors[result]
                    # if results are in folder, then duplicate files cause a problem!
                    # So add to 'deprecate' to move or remove before launching

    # Now make sure the steps can be found, and error out if not.
    for step in stepsToDo:
        app = STEPS[step]['app']
        dxApp = dxpy.find_data_objects(classname='file', name=app, name_mode='exact',
                                                         project=projectId, return_handler=False)
        if dxApp == None:
            print "ERROR: failure to locate app '"+app+"'!"
            sys.exit(1)

    return stepsToDo

def checkRunsPreviouslyLaunched(resultsFolder,projectId):
    '''Checks for currently running jobs and will exit if found.'''
    launchFilePath = resultsFolder + '/' + RUNS_LAUNCHED_FILE
    launchFids = dxencode.find_file(launchFilePath,projectId,multiple=True)
    if launchFids == None:
        print "  No prior jobs launched."
    else:
        # NOTE: Appending to the one file, but just in case handle multiple files.
        for fid in launchFids:
            with dxpy.open_dxfile(fid) as fd:
                for line in fd:
                    #print "Looking for job ["+line+"]"
                    runId = line.split(None,1)
                    if not runId[0].startswith('analysis-'):
                        continue
                    analysis = dxpy.DXAnalysis(dxid=runId[0])
                    if analysis == None:
                        continue
                    state = analysis.describe()['state']
                    # states I have seen: in_progress, terminated, done, failed
                    if state not in [ "done", "failed", "terminated" ]:
                        msg="Exiting: Can't launch because prior run ["+runId[0]+"] "
                        if len(runId) > 1:
                            msg+="("+runId[1]+") "
                        msg+= "has not finished (currently '"+state+"')."
                        print msg
                        sys.exit(1)
                    else:
                        msg="  Prior run ["+runId[0]+"] "
                        if len(runId) > 1:
                            msg+="("+runId[1]+") "
                        msg+= "is '"+state+"'."
                        print msg

def logThisRun(runId,resultsFolder,projectId):
    '''Adds a runId to the runsLaunched file in resultsFolder.'''
    # NOTE: DX manual lies?!  Append not possible?!  Then write new/delete old
    launchFilePath = resultsFolder + '/' + RUNS_LAUNCHED_FILE
    oldFid = dxencode.find_file(launchFilePath,projectId)
    newFh = dxpy.new_dxfile('a',project=projectId,folder=resultsFolder,name=RUNS_LAUNCHED_FILE)
    newFh.write(runId+' started:'+str(datetime.now())+'\n')
    if oldFid is not None:
        with dxpy.open_dxfile(oldFid) as oldFh:
            for oldRunId in oldFh:
                newFh.write(oldRunId+'\n')
        proj = dxpy.DXProject(projectId)
        proj.remove_objects([oldFid])
    newFh.close()


def createWorkflow(stepsToDo, priors, extras, resultsFolder, projectId, appProjectId=None):
    '''This function will populate a workflow for the stepsToDo.'''

    if len(stepsToDo) < 1:
        return None
    if appProjectId == None:
        appProjectId = projectId

    # create a workflow object
    wf = dxpy.new_dxworkflow(title=extras['name'],name=extras['name'],folder=resultsFolder,
                                            project=projectId,description=extras['description'])

    # NOTE: prevStepResults dict contains links to result files to be generated by previous steps
    prevStepResults = {}
    for step in stepsToDo:
        appName = STEPS[step]['app']
        app = dxencode.find_applet_by_name(appName, appProjectId)
        appInputs = {}
        # file inputs
        for fileToken in STEPS[step]['inputs'].keys():
            appInp = STEPS[step]['inputs'][fileToken]
            if fileToken in prevStepResults:
                appInputs[ appInp ] = prevStepResults[fileToken]
            elif fileToken in priors:
                if isinstance(priors[fileToken], list):
                    appInputs[ appInp ] = []
                    for fid in priors[fileToken]:
                        appInputs[ appInp ] += [ dxencode.get_file_link(fid) ]
                else:
                    appInputs[ appInp ] = dxencode.get_file_link(priors[fileToken])
            else:
                print "ERROR: step '"+step+"' can't find input '"+fileToken+"'!"
                sys.exit(1)
        # Non-file app inputs
        if 'params' in STEPS[step]:
            for param in STEPS[step]['params'].keys():
                appParam = STEPS[step]['params'][param]
                if param in extras:
                    appInputs[ appParam ] = extras[param]
                else:
                    print "ERROR: unable to locate '"+param+"' in extras."
                    sys.exit(1)
        # Add wf stage
        stageId = wf.add_stage(app, stage_input=appInputs, folder=resultsFolder)
        # outputs, which we will need to link to
        for fileToken in STEPS[step]['results'].keys():
            appOut = STEPS[step]['results'][fileToken]
            prevStepResults[ fileToken ] = dxpy.dxlink({ 'stage': stageId,'outputField': appOut })
    wfRun = wf.run({})
    return wfRun.describe()

#######################
def main():

    args = get_args()

    (AUTHID,AUTHPW,SERVER) = dxencode.processkey('default')
    url = SERVER + 'experiments/%s/?format=json&frame=embedded' %(args.experiment)
    response = dxencode.encoded_get(url, AUTHID, AUTHPW)
    exp = response.json()

    if not exp.get('replicates') or len(exp['replicates'] < 1):
        print "No replicates found in %s\n%s" % ( args.experiment, exp )
        sys.exit(1)

    rep = [ r for r in exp['replicates'] if r['biological_replicate_number'] == args.br and
                                            r['technical_replicate_number'] == args.tr ]

    replicate = "%s_%s" % (args.br, args.tr)
    if not rep or len(rep):
        print "No matching replicate found for %s" % replicate
        sys.exit(1)

    pairedEnd = rep.get('paired')
    try:
        library = rep['library']['accession']
        gender = 'female'
        if rep['library']['biosample']['sex'] == 'male':
            gender = 'male'
        organism = rep['library']['biosample']['donor']['organism']['name']
    except KeyError:
        print "Error, experiment %s replicate %s missing info\n%s" % (args.experiment, replicate, rep)

    if organism == 'mouse':
        genome = 'mm10'
    elif organism == 'human':
        genome = 'hg19'

    extras = pipelineSpecificExtras(genome,gender, args.experiment, replicate, library, pairedEnd)
    project = dxencode.get_project(args.project)
    projectId = project.get_id()

    args.resultsLoc = RESULT_FOLDER_DEFAULT + '/' + genome
    resultsFolder = args.resultsLoc + '/' + args.experiment + '/' + replicate
    if not args.test:
        if not dxencode.project_has_folder(project, resultsFolder):
            project.new_folder(resultsFolder,parents=True)

    print "Checking for prior results..."
    # Check if there are previous results
    # Perhaps reads files are already there?
    # NOTE: priors is a dictionary of fileIds that will be used to determine stepsToDo
    #       and fill in inputs to workflow steps
    priors = findPriorResults(pairedEnd,resultsFolder,projectId)

    print "Checking for read files..."
    # Find all reads files and move into place
    # TODO: files could be in: dx (usual), remote (url e.g.https://www.encodeproject.org/...
    #       or possibly local, Currently only DX locations are supported.
    reads1 = dxencode.find_and_copy_read_files(priors, args.reads1, args.test, 'reads1', resultsFolder, projectId)
    reads2 = dxencode.find_and_copy_read_files(priors, args.reads2, args.test, 'reads2', resultsFolder, projectId)

    print "Looking for reference files..."
    findReferenceFiles(priors,args.refLoc,extras)

    print "Determining steps to run..."
    # NOTE: stepsToDo is an ordered list of steps that need to be run
    deprecateFiles = [] # old results will need to be moved/removed if step is rerun
    stepsToDo = determineStepsToDo(pairedEnd, priors, deprecateFiles, projectId, force=args.force)

    # Report the plans
    print "Running '"+extras['title']+"'"
    print "     on "+extras['subTitle']
    if pairedEnd:
        print "- Reads1: "
    else:
        print "- Reads: "
    for fid in reads1:
        print "  " + dxencode.file_path_from_fid(fid)
    if pairedEnd:
        print "- Reads2: "
        for fid in reads2:
            print "  " + dxencode.file_path_from_fid(fid)
    print "- Reference files:"
    for token in GENOME_REFERENCES.keys():
        print "  " + dxencode.file_path_from_fid(priors[token],True)
    print "- Results written to: " + args.project + ":" +resultsFolder
    if len(stepsToDo) == 0:
        print "* All expected results are in the results folder, so there is nothing to do."
        print "  If this experiment/replicate needs to be rerun, then use the --force flag to "
        print "  rerun all steps; or remove suspect results from the folder before launching."
        sys.exit(0)
    else:
        print "- Steps to run:"
        steps = []
        if pairedEnd:
            steps = STEP_ORDER['pe']
        else:
            steps = STEP_ORDER['se']
        for step in steps:
            STEPS[step] = calculate_steps(step)
            if step in stepsToDo:
                print "  * "+STEPS[step]['app']+" will be run"
            else:
                if not step.find('concat') == 0:
                    print "    "+STEPS[step]['app']+" has already been run"

    print "Checking for currently running analyses..."
    checkRunsPreviouslyLaunched(resultsFolder,projectId)

    if len(deprecateFiles) > 0:
        if args.test:
            print "Would move "+str(len(deprecateFiles))+" prior result file(s) to '" + \
                                                                    resultsFolder+"/deprecated'."
            for fid in deprecateFiles:
                print "  " + dxencode.file_path_from_fid(fid)
        else:
            print "Moving "+str(len(deprecateFiles))+" prior result file(s) to '" + \
                                                                resultsFolder+"/deprecated'..."
            dxencode.move_files(deprecateFiles,resultsFolder+"/deprecated",projectId)

    # Exit if test only
    if args.test:
        print "TEST ONLY - exiting."
        sys.exit(0)

    print "Launch sequence initiating..."
    wfRun = createWorkflow(stepsToDo, priors, extras, resultsFolder,projectId)

    print "  We have liftoff!"
    logThisRun(wfRun['id'],resultsFolder,projectId)

    print "  Launched " + wfRun['id']
    print "(success)"

if __name__ == '__main__':
    main()

