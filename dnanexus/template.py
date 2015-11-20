#!/usr/bin/env python2.7
# template.py 1.0.0

import argparse, os, sys, subprocess, json
import commands
from datetime import datetime
from collections import deque

import dxpy

# The purpose of this module is to provide an alternative to dxencode/launch.py which does not rely upon encode at all.
# By making launchers derived from template, one can build all necessary dx applets and create a template workflow for
# the given pipeline.
#
# TODO: Ideally template.py => launch.py => dmeLaunch.py (fully functional) AND template.py => dmeLaunch.py (templating only).
#
# For a more complete description of the json structures expected, look at notes in:
#        https://github.com/ENCODE-DCC/dxencode/blob/master/launch.py.

# NOTES: This command-line utility will run the a pipeline for a single replicate or combined reps.
#      - All results will be written to a folder /<results_folder>/<exp_id>/rep<#>_<#>.
#      - If any results are already in that directory, then the steps that created those results
#        will not be rerun.
#      - If any jobs for the experiment and replicate are already running, nothing new will be
#        launched.
#      - Almost all the code is generic and derived classes will only need to extend or override a 
#        a few things.  It relies upon JSON pipeline descriptions which in simple cases can be 
#        directly read from dx (e.g. srnaLaunch.py).  More complex examples (e.g. lrnaLaunch.py)
#        must hard-code the pipeline json to handle repeated apps and dx name conflicts.
#        Tokens are used to abstract dx.app input/outout file names to avoid collisions.
#        - PIPELINE_BRANCHES[branch_id]["ORDER"] is the list of steps in the pipeline branch
#        - PIPELINE_BRANCHES[branch_id]["STEPS"] contain step definitions and enforces dependencies by:
#          'input' file tokens matching to 'result' file tokens of earlier steps.
#        - FILE_GLOBS is needed for locating result files from prior runs.
#      - Combined replicate processing is optional and depends upon PIPELINE_BRANCHES["COMBINED_REPS"], etc.
#        - The standard combination expected PIPELINE_BRANCH_ORDER[ "REP", "COMBINED_REPS" ].
#        - More complex combinations can be acheived by overriding add_combining_reps() in the derived class.
#        - Tying results from one branch to inputs into another relies on input file tokens ending with '_A', '_B' or '_ABC'.
#      - Control files are supported but most likely will require function overrides in derived
#        classes (e.g. rampageLaunch.py)
#      - In PIPELINE_BRANCHES[branch_id]["STEPS"], both 'inputs' and 'results' are key:value lists where the key is used
#        to attach results of one step to the inputs of another, while values must match dxapp.json names.
#      - The FILE_GLOBS are matched by input/result keys, not values.
#        (e.g. "input": { "pooled_bam": "bam" } matches file_glob "pooled_bam" and dxapp.json input named "bam".
#        This example would allow the same applet to be used on "replicate bams" and "poold bams" in the same pipeline.)
#      - RARE cases using 'output_values' and 'param_links' can support non-file outputs as inputs to later steps,
#        even across rep branch boundaries.
#
#   A design note about replicates in the launchers:
#      - "rep" objects should contain (almost) all information necesary to launch processing against them.
#        A few minor things common to the entire experiment are stored directly in the psv (pipeline specific variables) object
#      - Each "rep" object will get processed by one PIPELINE_BRANCH which determines which STEPS are run.
#      - All ENCODEd replicates ready for processing should get loaded into psv["reps"], keyed with letters: "a","b","c"...
#      - These "simple reps" have a "rep_tech" like "rep2_1" (for bio_rep=2, tech_rep=1).
#      - Combined replicates have a "rep_tech" like "reps1_1-1_2" (for combining simple reps rep1_1 and rep1_2).
#      - Combined replicates should be keyed (in sort order) after their "tributary" reps (e.g. "a","b","b-combines-a-b").
#      - The final combined replicate (aka "SEA" into which all rivers flow) is keyed with self.SEA_ID ("zzz")
#      - In the 'standard combination model' (which can be overridden in descendent classes): 
#        - Simple reps are processed by the "REP" pipeline branch.
#        - The one combined rep (aka "SEA") is processed by the "COMBINED_REPS" pipeline branch
#
# Example derived classes: 
#    https://github.com/ENCODE-DCC/long-rna-seq-pipeline/blob/master/dnanexus/lrnaLaunch.py
#    https://github.com/ENCODE-DCC/long-rna-seq-pipeline/blob/master/dnanexus/small-rna/srnaLaunch.py
#    https://github.com/ENCODE-DCC/long-rna-seq-pipeline/blob/master/dnanexus/rampage/rampageLaunch.py
#    https://github.com/ENCODE-DCC/dna-me-pipeline/blob/master/dnanexus/dmeLaunch.py
#    https://github.com/ENCODE-DCC/dnase_pipeline/blob/master/dnanexus/dnaseLaunch.py

class Launch(object):
    '''
    Template module is a simplified version of the Launch module that does not rely upon dxencode.py and
    stubs encodeD functionality.  When launchers are derived from Template instead of Launch, they can only 
    build applets and a template workflow with no inputs.  They cannot actually launch the worklows.
    '''

    PIPELINE_NAME = "MUST REPLACE IN DREIVED CLASS"
    ''' This must match the assay type returned by dxencode.get_assay_type({exp_id}).'''

    PIPELINE_HELP = "Launches '"+PIPELINE_NAME+"' pipeline " + \
                    "analysis for one replicate or combined replicates. "
    ''' This help title should name pipline and whether combined replicates are supported.'''
                    
    FOLDER_DEFAULT = '/runs/'
    ''' This the default location to place results folders for each experiment.'''
    
    PIPELINE_BRANCH_ORDER = None # MUST replace in derived class. Example: [ "REP", "COMBINED_REPS" ]
    '''A pipeline is frequently made of branches that flow into each other, such as replicate level to combined replicate.'''
    
    PIPELINE_BRANCHES = None # MUST replace in derived class. 'ORDER' and "STEPS" are required for each branch,
    # Example:               # though simple "STEPS" may be discovered.
    #{    "REP": {
    #            "ORDER": [ "rep_app_one", ... ],
    #     # NOTE: In simple cases "STEPS" can be discovered from dnanexus apps (example: srnaLaunch.py)
    #     # HINT: Try to discover the "STEPS" usubg build_simple_steps(verbose) to print the json, then modify from there  
    #            "STEPS": {
    #                        "rep_app_one": {
    #                            "inputs":  { "reads":    "reads" },
    #                            "app":       "rep_app_1", 
    #                            "params":  { "nthreads": "nthreads"}, 
    #                            "results": { "rep_bam":  "bam_out" }
    #                        }, ...
    #    },
    #    "COMBINED_REPS": ...
    #}
    '''Each branch must define the 'steps' and their (artificially) linear order.'''

    # "STEPS" objects within PIPELINE_BRANCHES contain the following:
    #    "step-name":   { # Must match "ORDER" and may match dx applet
    #        "app":     "dx-applet-name", # Must match applet
    #        "params":  { "param_token": "dx_input_name", ...}, # Non-files: Token must match psv key.
    #        "inputs":  { "input_token": "dx_input_name", ...},  # Files: Token *may* match name
    #        "results": { "result_token": "dx_output_name", ...}, # Files: Token *may* match name
    #        "output_values": { "output_token": "dx_output_name" }, # RARE: non-file output used as input to another step
    #        "param_links": { "param_token": {"name":"output_token", "rep":"sister"} }} # RARE: link input to non-file output 
    # Note, tokens *may* be the same as dx names, but MUST MATCH other tokens to define dependencies as: 
    #     stepA:result_token == stepB:input_token
    # Also result_tokens must be in FILE_GLOB.keys()
    #
    # NOTE: Simple pipelines have no reused steps or dx name to param complications and the 
    #       "STEPS" objects can be generated directly from dx apps.  If unsure, print the json from
    #       self.build_simple_steps(verbose=True) to determine if it will work for your pipeline.
    
    PRUNE_STEPS = None
    '''Rare setting that allows dynamically pruning steps from a pipeline befpre running it.  Used by lrnaLaunch.py.'''
    
    FILE_GLOBS = {"MAY": {}, "NEED": {}, "TO": {}, "REPLACE": {} }
    '''Globs for discovering existing results. Can be discoverd in simple piplines.'''

    CONTROL_FILE_GLOB = None
    '''Glob for discovering control files if experiment has controls. Othewise: None'''

    REFERENCE_FILES = {
        # For looking up reference file names.
        }

    GENOMES_SUPPORTED = ['hg19', 'mm10']
    GENOME_DEFAULT = 'hg19'
    ''' This is the default Genome that experiments are mapped to.'''
    
    STANDARD_LABELS = {
        'long-rna-seq':  { 'long': 'Long-RNA-seq',  'short': 'lrna',  'se_and_pe': True  }, 
        'small-rna-seq': { 'long': 'Small-RNA-seq', 'short': 'srna',  'se_and_pe': False }, 
        'rampage':       { 'long': 'Rampage',       'short': 'ramp',  'se_and_pe': False }, 
        'dnase-seq':     { 'long': 'DNase-seq',     'short': 'dnase', 'se_and_pe': True  },
        'dna-me':        { 'long': 'DNA-methlyation','short': 'dme',   'se_and_pe': True  },
    }
    '''Standard labelling requires exp_type specific labels.  This can be overridden in descendent classes.'''
    
    SEA_ID = 'zzz'
    '''SEA is the final branch into which all tributaries flow.'''

    def __init__(self):
        '''
        Launch expects to be the base class for derived pipeline specific launch classes.
        As such it will not run stand alone.
        an experiment.
        '''
        self.args = {} # run time arguments
        self.template = True
        self.build_apps = False
        self.no_refs = False
        self.proj_name = None
        self.project = None
        self.proj_id = None
        self.exp = {}  # Will hold the encoded exp json
        self.psv = {} # will hold pipeline specific variables.
        self.multi_rep = False       # This run includes more than one replicate
        self.combined_reps = False   # This run includes combined replicate workflows
        self.link_later = None       # Rare: when 2 sister branches link to each other, it requires a final wf pass.
        self.deprecate = []          # Master list of files to deprecate.  Needed to cross rep boundaries in steps_to_run
        print # TEMPORARY: adds a newline to "while retrieving session configuration" unknown error
    
    def get_args(self,parse=True):
        '''Parse the input arguments.'''
        # MAY NEED TO EXTEND OR REPLACE in derived class
        # see long-rna-seq: lrnaLaunch.py for and example of extending
        
        # Change this description to PIPELINE SPECIFIC version
        ap = argparse.ArgumentParser(description=self.PIPELINE_HELP +
                    "Can be run repeatedly and will launch only the steps that are needed to " +
                    "finish the pipeline. All results will be placed in the folder: " +
                    "/<resultsLoc>/<experiment>/<replicate>.")

        ap.add_argument('-r','--reps','--replicates',
                        help="Request one or more replicates (e.g.'rep1_1 2 2_2').'",
                        nargs='+',
                        required=False)

        ap.add_argument('-br','--br', '--biological-replicate',
                        help="Request single biological replicate number (--tr==1 implied)",
                        type=int,
                        #default='1',
                        required=False)

        ap.add_argument('-tr','--tr', '--technical-replicate',
                        help="Select single technical replicate number (default: 1, must select --br as well)",
                        type=int,
                        default='1',
                        required=False)

        ap.add_argument('--project',
                        help="Project to run analysis in (default: '" + \
                                                        self.env_get_current_project() + "')",
                        required=False)

        ap.add_argument('-f','--folder',
                        help="The location to to place results folders (default: '<project>:" + \
                                                                  self.FOLDER_DEFAULT + "')",
                        default=self.FOLDER_DEFAULT,
                        required=False)

        ap.add_argument('-g','--genome',
                        help="The genome assembly to run on (default: '<project>:" + \
                                                                  self.GENOME_DEFAULTS['human'] + "')",
                        default=None,
                        required=False)

        ap.add_argument('--run',
                        help='Run the workflow after assembling it.',
                        action='store_true',
                        required=False)

        ap.add_argument('--build_apps',
                        help='Build applets using build_apps script before building workflow.',
                        action='store_true',
                        required=False)

        ap.add_argument('--template',
                        help='Build a template workflow only.',
                        action='store_true',
                        required=False)

        ap.add_argument('--pe',
                        help='For template only: make paired-end mapping template.',
                        action='store_true',
                        required=False)

        ap.add_argument('--no_refs',
                        help='For template only: do not use standard ENCODE reference files as inputs.',
                        action='store_true',
                        required=False)

        ap.add_argument('--test',
                        help='Test run only, do not launch anything.',
                        action='store_true',
                        required=False)

        ap.add_argument('--force',
                        help='Force rerunning all steps.',
                        action='store_true',
                        required=False)

        if parse:
            return ap.parse_args()
        else:
            return ap
            
########### dxencode.py transplants:
    PRODUCTION_PROJECT = "ENCODE - Production runs"

    REF_PROJECT_DEFAULT = 'ENCODE Reference Files'
    ''' This the default DNA Nexus project to find reference files in.'''

    REF_FOLDER_DEFAULT = '/'
    ''' This the default folder that reference files are found in.'''

    GENOME_DEFAULTS = { 'human': 'hg19', 'mouse': 'mm10' }
    ''' This the default genomes for each supported organism.'''

    def env_get_current_project(self):
        ''' Returns the current project name for the command-line environment '''
        err, proj_name = commands.getstatusoutput('cat ~/.dnanexus_config/DX_PROJECT_CONTEXT_NAME')
        if err != 0:
            return None
        return proj_name
        
    def get_project(self, projectName, level=None):
        '''Returns the DXProject by name or errors out if not found.'''
        try:
            project = dxpy.find_one_project(name=projectName, name_mode='exact',
                                            level=level, return_handler=False)
        except:
            print "Could not find 1 and only 1 project named '"+projectName+"'."
            sys.exit(1)

        return dxpy.DXProject(project['id'])

    def find_applet_by_name(self, applet_name, applets_project_id):
        '''Looks up an applet by name in the project that holds tools.  From Joe Dale's code.'''
        return dxpy.find_one_data_object(classname="applet", name=applet_name,project=applets_project_id,
                                                                    zero_ok=False, more_ok=False, return_handler=True)

    def folder_normalize(self,folder,starting=True,trailing=True):
        '''Normalizes a folder to always begin with and end with '/'.'''
        if starting:
            if not folder.startswith('/'):
                folder = '/' + folder
        else:
            if folder.startswith('/'):
                folder = folder[1:]
        if trailing:
            if not folder.endswith('/'):
                folder = folder + '/'
        else:
            if folder.endswith('/'):
                folder = folder[:-1]
        return folder

    def umbrella_folder(self,folder,default,proj_name=None,exp_type=None,genome=None,annotation=None):
        '''Returns a normalized umbrella folder (that holds the experiments of a given type).'''
        if folder != default:
            return self.folder_normalize(folder)
                
        # No change to default, so build from parts if available
        if exp_type == None:
            return self.folder_normalize(folder)

        if proj_name == self.PRODUCTION_PROJECT:
            if exp_type == "long-rna-seq":
                folder = "/long-RNA-seq/runs/"
            elif exp_type == "small-rna-seq":
                folder = "/small-RNA-seq/runs/"
            elif exp_type == "dnase-seq":
                folder = "/DNAse-seq/runs/"
            elif exp_type == "dna-me":
                folder = "/WG Bisulfite (Methylation)/runs/"
            elif exp_type == "chip-seq":
                folder = "/ChIP-seq/runs/"
            else:
                folder = "/" + exp_type + '/runs/'
        else:
            if exp_type == "long-rna-seq":
                folder = "/lrna/"
            elif exp_type == "small-rna-seq":
                folder = "/srna/"
            elif exp_type == "dnase-seq":
                folder = "/dnase/"
            elif exp_type == "dna-me":
                folder = "/dme/"
            else:
                folder = "/" + exp_type + '/'

        if genome != None:
            folder +=  genome + '/'
            if annotation != None and genome == 'mm10':
                folder += annotation + '/'
                
        return folder

    FILES = {} ## Dict to cache files
    
    def find_file(self,filePath,project=None,verbose=False,multiple=False, recurse=True):
        '''Using a DX style file path, find the file.'''
        proj = project
        path = filePath
        fileName = filePath
        if filePath.find(':') != -1:
            proj, path = filePath.split(':', 1)
        if path.rfind('/') != -1:
            path, fileName = path.rsplit('/', 1)
        else:
            fileName = path
            path = '/'
        if proj == None:
            if verbose:
                print "ERROR: Don't know what project to use for '" + path + "'."
            return None
        if proj.find('project-') == 0:
            projId = proj
        else:
            projId = self.get_project(proj, level='VIEW').get_id()
        mode = 'exact'
        if filePath.find('*') or filePath.find('?'):
            mode = 'glob'
        fileDicts = list(dxpy.find_data_objects(classname='file', folder=path, name=fileName, recurse=recurse,
                                                name_mode=mode, project=projId, return_handler=False))

        if fileDicts == None or len(fileDicts) == 0:
            #print "- Found 0 files from '" + proj + ":" + filePath + "'."
            if verbose:
                print "ERROR: Failed to find '" + proj + ":" + filePath + "'."
            return None
        elif len(fileDicts) > 1 or multiple:
            #print "- Found "+str(len(fileDict))+" files from '" + proj + ":" + filePath + "'."
            if not multiple:
                if verbose:
                    print "ERROR: Found "+str(len(fileDicts))+" files when expecting 1 '" + proj + ":" + filePath + "'."
                return None
            fids = []
            for fileDict in fileDicts:
                self.FILES[fileDict['id']] = dxpy.dxlink(fileDict)
                fids.append( fileDict['id'] )
            return fids
        else:
            #print "- FOUND '" + proj + ":" + filePath + "'."
            self.FILES[fileDicts[0]['id']] = dxpy.dxlink(fileDicts[0])
            return fileDicts[0]['id']

    def description_from_fid(self,fid,properties=False):
        '''Returns file description object from fid.'''
        try:
            dxlink = self.FILES[fid]
        except:
            #logger.error("File %s not cached, trying id" % fid)
            dxlink = fid

        return dxpy.describe(dxlink,incl_properties=properties)

    def file_path_from_fid(self, fid, projectToo=False):
        '''Returns full dx path to file from a file id.'''
        fileDict = self.description_from_fid( fid )
        if fileDict['folder'] == '/':
            path = '/' + fileDict['name']
        else:
            path = fileDict['folder'] + '/' + fileDict['name']
        if projectToo:
            projDict = dxpy.describe(fileDict['project'])
            path = projDict['name'] + ':' + path
        return path

    def project_has_folder(self, project, folder):
        ''' Checks for a folder in a given DX project '''
        try:
            found = project.list_folder(folder)
        except:
            return False
        return True

########### dxencode.py transplants:

    def pipeline_specific_vars(self,args,verbose=False):
        '''Adds pipeline specific variables to a dict, for use building the workflow.'''
        # PLEASE EXTEND in derived class

        # psv can contain any variables, but it must contain these at a minimum:
        # - Any non-file input param needed to launch the workflow
        # - 'resultFolder' - full dx path (without project) to the results folder of the specific run
        # - 'name' - A short name used for specific workflow run. 
        # - 'description' - generic description of the pipeline.
        # - 'title'/['subtitle'] for command line output announcing what will be done 
        
        # Start with dict containing common variables
        psv = self.common_variables(args)
        
        # Now add pipline specific variables here...

        if verbose:
            print "Pipeline Specific Vars:"
            print json.dumps(psv,indent=4)
        return psv


    def find_ref_files(self,priors):
        '''Locates all reference files based upon organism and gender.'''
        # PLEASE REPLACE in derived class
        assert "must replace!" and False
        
        # Add your pipeline specific referernce file detection here
        
        # Examples:
        # If all the ref files are in json keyed by ['genome']['gender'] then try:
        for ref in self.REFERENCES_FILES.keys():
            dxfile = self.psv['refLoc']+self.REFERENCES_FILES[ref][self.psv['genome']][self.psv['gender']]
            fid = self.find_file(dxfile,self.REF_PROJECT_DEFAULT)
            if fid == None:
                sys.exit("ERROR: Unable to locate ref file: '" + dxfile + "'")
            else:
                priors[ref] = fid
        
        # Or explicitly by key
        chrom_sizes = self.psv['refLoc']+self.REFERENCE_FILES['chrom_sizes'][self.psv['genome']][self.psv['gender']]
        chrom_sizes_fid = self.find_file(chrom_sizes,self.REF_PROJECT_DEFAULT)
        if chrom_sizes_fid == None:
            sys.exit("ERROR: Unable to locate Chrom Sizes file '" + chrom_sizes + "'")
        else:
            priors['chrom_sizes'] = chrom_sizes_fid
        self.psv['ref_files'] = self.REFERENCE_FILES.keys()
        

    def add_combining_reps(self, psv):
        '''Defines how replicated are combined.'''
        # MAY NEED TO REPLACE for pipelines with combined steps
        
        reps = psv['reps']
        # In the 'standard combining model' PIPELINE_BRANCH_ORDER = [ "REP", "COMBINED_REPS" ]
        # and all ENCODEd replicates are in psv['reps'] keyed as 'a','b',... and having rep['rep_tech'] = 'rep1_1'
        # All these simple reps will have rep['branch_id'] = "REP"
        
        # If there are non-standard combinations, such as tech_reps combining, then bio_reps combining, define it here...
        # HINT: dnaseLaunch.py has a working example

        # If requested with -cr, OR if exactly 2 reps from different bio_reps were found then self.combined_reps == True
        # In the 'standard combining model' the one combined replicate keyed with self.SEA_ID is created below:
        # OVERRIDE as you will!
        if 'COMBINED_REPS' in self.PIPELINE_BRANCH_ORDER and self.combined_reps:
            assert len(reps.keys()) == 2
            assert 'a' in reps.keys()
            assert 'b' in reps.keys()
            
            sea = {} # SEA is the final branch into which all tributaries flow
            sea['branch_id'] = 'COMBINED_REPS'
            sea['tributaries'] = ['a','b']
            sea['rep_tech'] = 'reps' + reps['a']['rep_tech'][3:] + \
                                 '-' + reps['b']['rep_tech'][3:]  # combined delimited by '-'
            # Special case of 2 allows for designating sisters
            psv['reps'][sea['tributaries'][0]]['sister'] = sea['tributaries'][1]
            psv['reps'][sea['tributaries'][1]]['sister'] = sea['tributaries'][0]
            psv['reps'][self.SEA_ID] = sea
            psv['rep_tech'] = sea['rep_tech']  # psv gets labelled the same as the sea
            
        
    def find_combined_inputs(self,rep,steps,file_globs,verbose=False):
        '''Finds the inputs for a combined run in the directories of the replicate runs.'''
        # MAY NEED TO REPLACE for pipelines with combined steps
        #verbose=True
        
        if not self.combined_reps or 'tributaries' not in rep:
            return

        print "Detecting combined-replicate inputs for %s..." % rep['rep_tech']
        inputs = {}
        for step in rep['path']:
            if verbose:
                print "DEBUG: step: "+ step 
            for file_token in steps[step]['inputs'].keys():
                if verbose:
                    print "DEBUG:   file_token: "+ file_token
                if file_token[-2] == '_':
                    letter = file_token[-1].lower()
                    if letter in rep['tributaries']:
                        tributaries = [ letter ]
                    elif rep['tributaries'][0][1] == '-': # second level combining? (e.g. tributaties = ['b-biorep','e-biorep'])
                        # Difficulty is that inp_a/inp_b are being matched to ['b-biorep','e-biorep'] by position
                        letters = 'abcdefghijklmnopqrstuvwxyz'
                        ix = letters.find(letter) # hopefully returns 0, then 1
                        if ix != -1 and ix < len(rep['tributaries']):
                            tributaries = [ rep['tributaries'][ix] ]  # Too clever by half?
                    if len(tributaries) == 0:
                        print "*** " + rep_key + " not found in rep['tributaries']"
                        continue
                elif file_token.endswith('_ABC'):  # Also too clever by 3/4ths: ending says all tributaries contribute to set.
                    tributaries = rep['tributaries']
                else:
                    continue
                if verbose:
                    print "DEBUG:   - found "+str(len(tributaries))+" tributaries."
                rep['priors'][file_token] = []
                inputs[file_token] = []
        
        if inputs:
            if verbose:
                print "DEBUG: inputs:"
                print json.dumps(inputs,indent=4,sort_keys=True)
                #sys.exit(1)
        if 'inputs' not in rep:
            rep['inputs'] = inputs
        else:
            rep['inputs'].update( inputs )            
    

    ############## NOT EXPECTED TO OVERIDE THE FOLLOWING METHODS ############
    def common_variables(self,args,fastqs=True):
        '''Initializes dict with common variables from args and encoded.'''
        # NOT EXPECTED TO OVERRIDE
        
        self.test = args.test
        if args.no_refs:
            self.no_refs = True
        if args.build_apps:
            self.build_apps = True
            
        cv = {}
        self.proj_name = self.env_get_current_project()
        if self.proj_name == None or args.project != None:
            self.proj_name = args.project
        if self.proj_name == None:
            print "Please enter a '--project' to run in."
            sys.exit(1)
        self.project = self.get_project(self.proj_name)
        self.proj_id = self.project.get_id()        

        cv['project']    = self.proj_name
        cv['experiment'] = 'template'
        
        print "Templating experiment specifics..."
        cv['exp_type'] = self.PIPELINE_NAME  # By definition, this should be correct
        self.load_template_reps(args, cv)
        
        assert 'a' in cv['reps']

        # Only supported genomes
        cv['gender'] = cv['reps']['a']['sex']
        organism = cv['reps']['a']['organism']
        if organism in self.GENOME_DEFAULTS:
            cv['genome'] = self.GENOME_DEFAULTS[organism]
            if args.genome != None and args.genome != self.GENOME_DEFAULTS['human']:
                if args.genome not in self.GENOMES_SUPPORTED:
                    print "Genome %s not currently supported" % args.genome
                    sys.exit(1)
                cv['genome'] = args.genome            
        else:
            print "Organism %s not currently supported" % organism
            sys.exit(1)

        # Paired ends?  Make sure combined reps reflect their tributaries and exp reflects all
        if self.combined_reps:
            cv['paired_end'] = False
            for rep_id in sorted( cv['reps'].keys() ):  # Must be in order
                if len(rep_id) == 1: # Simple reps are the tributaries and should have "paired_end" set already
                    assert "paired_end"      in cv['reps'][rep_id]
                    assert "tributaries" not in cv['reps'][rep_id]
                    continue
                river = cv['reps'][rep_id]
                assert "tributaries" in river
                paired = None
                for trib_id in river['tributaries']:
                    tributary = cv['reps'][trib_id]
                    if paired == None:
                        first_rep_tech = tributary['rep_tech']
                        paired = tributary['paired_end']
                        if 'paired_end' not in river:
                            river['paired_end'] = paired
                        if paired:
                            cv['paired_end'] = paired
                    elif paired != tributary['paired_end']:
                        print "- WARNING: Combining replicates '"+first_rep_tech+"' and '"+tributary['rep_tech']+"' " + \
                                        "are expected to be all paired-end or single-end!"
                        river['paired_end'] = False  # Mixes will be treated as single end!
                        #sys.exit(1)
        else:  # for independent simple replicates, top level doesn't really matter.
            cv['paired_end'] = cv['reps']['a']['paired_end']

        # Default locations
        cv['refLoc'] = self.REF_FOLDER_DEFAULT
        if not self.no_refs:
            cv['refLoc'] += cv['genome'] + '/'
        genome = cv['genome']
        if self.no_refs:
            genome = None # If no refs then hide genome
        cv['resultsLoc'] = self.umbrella_folder(args.folder,self.FOLDER_DEFAULT,self.proj_name,cv['exp_type'],genome)
        cv['resultsFolder'] = cv['resultsLoc']
        self.update_rep_result_folders(cv)
        
        # Standard labels:
        if cv['exp_type'] in self.STANDARD_LABELS:
            labels = self.STANDARD_LABELS[cv['exp_type']]
            cv['description'] = "The ENCODE "+labels['long']+" pipeline"
            cv['title'] = labels['long']
            cv['name'] = labels['short']
            if not self.no_refs:
                cv['name'] = labels['short'] + "_"+cv['genome']
                if not self.template:
                    if cv['gender'] == 'female':
                        cv['name'] += "XX"
                    else:
                        cv['name'] += "XY"
            if labels['se_and_pe']:
                if cv['paired_end']:
                    cv['title'] += " paired-end"
                    cv['name'] += "PE"
                else:
                    cv['title'] += " single-end"
                    cv['name']  += "SE"
            cv['title'] += " " + cv['experiment']
            cv['name']  += "_" + cv['experiment']
            for ltr in cv['reps'].keys():
                rep = cv['reps'][ltr]
                rep['title'] = cv['title'] + " - "+rep['rep_tech']
                if 'library_id' in rep:
                    rep['title'] += " (library '"+rep['library_id']+"')"
                if not self.no_refs: 
                    rep['title'] += " on "+cv['genome']
                    if not self.template:
                        rep['title'] += ", "+cv['gender']
                rep['name']  = cv['name']+"_"+rep['rep_tech']
            cv['title']   += " - "+cv['rep_tech']
            if not self.no_refs: 
                cv['title']   += " on "+cv['genome']
                if not self.template:
                    cv['title']   += ", "+cv['gender']
            cv['name']    += "_"+cv['rep_tech']

        return cv

    def load_template_reps(self, args, cv):
        '''
        Creates template replicates building the cv['reps'] list from requested or all.
        '''
        # NOT EXPECTED TO OVERRIDE
        
        # A design note about replicates in the launchers:
        # - "rep" objects should contain (almost) all information necesary to launch processing against them.
        #   A few minor things common to the entire experiment are stored directly in the psv (pipeline specific variables)
        # - Each "rep" object will get processed by one PIPELINE_BRANCH which determines which STEPS are run.
        # - All ENCODEd replicates ready for processing should get loaded into psv["reps"], keyed with letters: "a","b","c"...
        # - These "simple reps" have a "rep_tech" like "rep2_1" (for bio_rep=2, tech_rep=1).
        # - Combined replicates have a "rep_tech" like "reps1_1-2_1" (for combining simple reps rep1_1 and rep2_1).
        # - Combined replicates should be keyed (in sort order) after their "tributary" reps (e.g. "a","b","b-combines-a-b").
        # - The final combined replicate (aka "SEA" into which all rivers flow) is keyed with self.SEA_ID ("zzz")
        # - In the 'standard combination model' (which can be overridden in descendent classes): 
        #   - Simple reps are processed by the "REP" pipeline branch.
        #   - The one combined rep (aka "SEA") is processed by the "COMBINED_REPS" pipeline branch

        # Determine rep_techs requested
        rep_techs = []
        if 'reps' in args and args.reps != None:
            for rep_tech in args.reps:
                # Normalize rep_tech string
                if '_' not in rep_tech:
                    rep_tech = rep_tech + '_1' # default to tech_rep 1
                if not rep_tech.startswith('rep'):
                    rep_tech = 'rep' + rep_tech
                rep_techs.append( (rep_tech,int(rep_tech[3:].split('_')[0]),int(rep_tech[3:].split('_')[1])) )
        elif args.br != None and args.br != 0:
            rep_tech = 'rep' + str(args.br) + '_' + str(args.tr)
            rep_techs.append( (rep_tech,int(args.br),args(args.tr)) )
        else:
            rep_techs = [ ('rep1_1',1,1), ('rep2_1',2,1) ]  # Default is 2 biological replicate pipeline run
            
        # Create rep dicts that will be used to build pipeline ins and outs
        cv_reps = {}
        letters = deque('abcdefghijklmnopqrstuvwxyz')
        processed = []
        for (rep_tech,br,tr) in rep_techs:
            if rep_tech in processed:
                print >> sys.stderr, "ERROR: requested %s twice." % rep_tech
                sys.exit(1)
            processed.append(rep_tech)
            rep = { 'br': br, 'tr': tr,'rep_tech': rep_tech }
            rep['organism'] = 'human'
            if args.genome == 'mm10':
                rep['organism'] = 'mouse'
            rep['sex'] = 'male' # Sorry, but the patriarchy has ruled that we must include chrY.
            rep['paired_end'] = args.pe  # Defaults to se        
            rep['has_reads'] = False
                    
            # Load read and control files, only if requested.
            run_type = None  # The files should be consistent as all single-end or paired-end, but some mapping got it wrong
            rep['fastqs'] = { "1": [], "2": [] }
            rep['controls'] = []
            rep['priors'] = {} 
            rep['inputs'] = {}
            rep['inputs']['Reads1'] = []
            rep['inputs']['Reads2'] = []
                
            ltr = letters.popleft()
            cv_reps[ltr] = rep
            #print "* Added: cv_reps[%s] as %s" % (ltr, rep_tech)
            
        # Add combining replicates as appropriate
        # Assume combined reps if supported AND exactly 2 reps AND for different biological reps
        if self.PIPELINE_BRANCH_ORDER != None and 'COMBINED_REPS' in self.PIPELINE_BRANCH_ORDER:
            if len(cv_reps.keys()) == 2 and cv_reps['a']['br'] != cv_reps['b']['br']:
                self.combined_reps = True
        self.multi_rep = (len(cv_reps.keys()) > 1)
            
        # A little more rep tidying
        for ltr in cv_reps.keys():
            rep = cv_reps[ltr]
            rep['branch_id'] = 'REP'  # MUST override for non-standard pipelines which have ['REP','COMBINED_REPS']
            rep['concat_id'] = 'reads'
            if 'paired_end' in rep and rep['paired_end']: 
                rep['concat_id2'] = 'reads2'
                
        # mult-rep rep_tech: 
        if self.combined_reps:
            cv['rep_tech'] = 'reps' + cv_reps['a']['rep_tech'][3:] + \
                                '-' + cv_reps['b']['rep_tech'][3:]  # combined delimited by '-'
        elif self.multi_rep:
            cv['rep_tech'] = 'reps:' + cv_reps['a']['rep_tech'][3:]
            for ltr in sorted(cv_reps.keys()):
                if ltr != 'a':
                    cv['rep_tech'] += ','+cv_reps[ltr]['rep_tech'][3:] # not combined are delimited by ','
        else:
            cv['rep_tech'] = cv_reps['a']['rep_tech']
            
        cv['reps'] = cv_reps 

        # Call override-able function to define replicate combinations.
        self.add_combining_reps(cv)
        for rep_id in sorted( cv['reps'].keys() ):  # Everything will be initiated
            rep = cv['reps'][rep_id]
            if "priors" not in rep:
                rep['priors'] = {} 
                rep['inputs'] = {}


    def update_rep_result_folders(self, cv):
        ''' create input object for a step and extends the file_globs dict as appropriate.'''
        # NOT EXPECTED TO OVERRIDE
        for rep_id in sorted( cv['reps'].keys() ):
            rep = cv['reps'][rep_id]
            rep['resultsFolder'] = cv['resultsFolder']


    def build_applets_if_necessary(self):
        ''' create input object for a step and extends the file_globs dict as appropriate.'''
        # Build apps first (if necessary)...
        if self.build_apps:
            build_cmd = os.path.dirname(sys.argv[0]) + '/build_applets'
            # TODO: make sure to use the correct project self.proj_id
            if self.test:
                print "Requires building the applets with '%s' before this can succeed!" % build_cmd
            else:
                # TODO: Just request the applets being used.
                print "Requesting to build the applets with '%s' ... " % build_cmd
                subprocess.call([build_cmd])
                self.build_apps = False # Done with that!

    def build_a_step(self, applet, file_globs, proj_id):
        ''' create input object for a step and extends the file_globs dict as appropriate.'''
        # NOT EXPECTED TO OVERRIDE
        try:
            dxapp = dxpy.find_one_data_object(classname="applet", name=applet, project=proj_id,
                                              zero_ok=False, more_ok=False,describe=True)
        except IOError:
            print "Cannot find applet '"+applet+"' in "+proj_id
            sys.exit(1)

        params = {}
        inputs = {}
        results = {}
        inps = dxapp['describe'].get('inputSpec') or []
        for inp in inps:
            in_name = inp['name'].encode('ascii','ignore')
            if inp['class'] == 'file' or inp['class'] == 'array:file':
                inputs[in_name] = in_name
            else:
                params[in_name] = in_name

        outs = dxapp['describe'].get('outputSpec') or []
        for out in outs:
            if out['class'] == 'file' or out['class'] == 'array:file':
                out_name = out['name'].encode('ascii','ignore')
                results[out_name] = out_name
                if out_name not in file_globs and 'patterns' in out:
                    file_globs[out_name] = '/'+out['patterns'][0].encode('ascii','ignore')
            else:
                pass
                # TODO not sure what to do with these

        return {
            'app': applet,
            'params': params,
            'inputs': inputs,
            'results': results
        }


    def build_simple_steps(self,pipe_path, proj_id, verbose=False):
        '''
        Builds dict of steps for the apps in the pipeline and a dict of file_globs for look up.
        Only works for pipelines where every step in pipe_path is a distinct app,
        and each result glob is uniq for all pipeline results.
        '''
        # NOT EXPECTED TO OVERRIDE
        self.build_applets_if_necessary()
        steps = {}
        file_globs = {}
        for step in pipe_path:
            steps[step] = self.build_a_step(step, file_globs, proj_id)
        if verbose:
            print "STEPS = "
            print json.dumps(steps,indent=4)
            print "FILE_GLOBS = "
            print json.dumps(file_globs,indent=4)

        return [ steps, file_globs ]


    def assemble_steps_and_globs(self,branch_id,branch):
        '''Assemble the requested steps and file globs.'''
        # NOT EXPECTED TO OVERRIDE
        assert 'ORDER' in branch
        pipe_path = branch['ORDER']    

        if self.FILE_GLOBS != None and 'STEPS' in branch:
            # Use for complex pipelines:
            file_globs = self.FILE_GLOBS
            pipe_steps = branch['STEPS']
        else:
            # Use for simple pipelines
            pipe_steps, file_globs = self.build_simple_steps(pipe_path,self.proj_id)
            # NOTE: Useful for discovering pipeline json, which can be copied/modified for new launcher
            #pipe_steps, file_globs = self.build_simple_steps(pipe_path,self.proj_id,verbose=True)
            #sys.exit(1)

        for rep in self.psv['reps'].values():
            if rep["branch_id"] != branch_id:
                continue
            if type(pipe_path) == dict:
                assert 'se' in pipe_path and 'pe' in pipe_path
                if rep['paired_end']:
                    rep['path'] = pipe_path['pe']
                else:
                    rep['path'] = pipe_path['se']
            else:
                rep['path'] = pipe_path
            # In a rare case, pruning of steps is done.            
            if self.PRUNE_STEPS != None and len(self.PRUNE_STEPS) != 0:
                for step in self.PRUNE_STEPS:
                    if step in pipe_steps:
                        del pipe_steps[step]
                    if step in rep['path']:
                        rep['path'].remove(step)
            rep['steps'] = pipe_steps
            
        return [ pipe_steps, file_globs ]

    def find_inputs_and_priors(self,branch_id,globs):
        '''Finds the inputs and priors for all reps for a given branch of the pipeline.'''
        # NOT EXPECTED TO OVERRIDE
        
        # Templating only requires linking up to combined inputs            
        for rep in self.psv['reps'].values():
            # If rep has tributaries, then look for those priors!
            if rep['branch_id'] != branch_id:
                continue
            if 'tributaries' in rep:
                self.find_combined_inputs(rep,rep['steps'],globs)


    def find_all_ref_files(self):
        '''Locates all reference files based upon organism and gender.'''
        # NOT EXPECTED TO OVERRIDE
        
        if not self.no_refs:
            print "Looking for reference files..."
            ref_priors = {}
            self.find_ref_files(ref_priors)
        
            # Add into all rep priors
            for rep in self.psv['reps'].values():
                rep['priors'].update(ref_priors)


    def determine_steps_to_run(self,pipe_path, steps, priors, deprecate, force=False, verbose=False):
        '''Determine what steps need to be done, base upon prior results.'''
        # NOT EXPECTED TO OVERRIDE
        #verbose=True  # Very useful when verifying new/updated launcher 
        will_create = []
        steps_to_run = []
        for step in pipe_path:
            steps_to_run += [ step ]
            if verbose:
                print "- Adding step '"+step+"' because all steps are steps to run when templating."

        # Now make sure the steps can be found, and error out if not.
        self.build_applets_if_necessary()
        for step in steps_to_run:
            app = steps[step]['app']
            dxApp = dxpy.find_data_objects(classname='file', name=app, name_mode='exact',
                                                        project=self.proj_id, return_handler=False)
            if dxApp == None:
                print "ERROR: failure to locate app '"+app+"'!"
                sys.exit(1)

        return steps_to_run


    def determine_steps_needed(self, force=False):
        '''Determine steps needed for replicate(s) and combined, base upon prior results.'''
        # NOT EXPECTED TO OVERRIDE
        
        print "Determining steps to run..."
        # NOTE: stepsToDo is an ordered list of steps that need to be run
        for rep_id in sorted( self.psv['reps'].keys() ):
            rep = self.psv['reps'][rep_id]
            rep['deprecate'] = [] # old results will need to be moved/removed if step is rerun
            rep['stepsToDo'] = self.determine_steps_to_run(rep['path'], rep['steps'], \
                                                    rep['priors'], rep['deprecate'], force=force)


    def find_results_from_prev_branch(self, river, inp_token, expect_set):
        '''If combined_reps step, look for input from results of replicate level step.'''
        # NOT EXPECTED TO OVERRIDE
        
        if not self.combined_reps:    # Not even doing a wf with combining branches 
            return None
        if 'tributaries' not in river:  # Not a combining branch
            return None
        if inp_token not in self.FILE_GLOBS: # assertable
            return None
        # The combined_step input tokens should either end in '_a', '_b' or the like,
        # Or the step expects to combine N files so expect_set=True and all matching results should be combined.
        
        # Using letters since 2 tech_rep branches could be combined into 2 bio_rep branches and then 1 experiment level branch
        # Thus tr branch_ids would be: 'a','b','c','d'
        #  and br branch_ids might be 'b-bio_rep1','d-bio_rep2', and those branches flow into the sea (self.SEA_ID='zzz').
        # So if a BRANCH STEP for combined rep 'b-bio_rep1' expects 2 inputs ending '_a' and '_b then
        #       the same STEP for combined rep 'd-bio_rep2' would expect the same 2 endings, not '_c' and '_d'
        letters = deque('abcdefghijklmnopqrstuvwxyz')
        looking_for = None
        if expect_set:
            # Expectation that when combining results, dx app's input names either end in '_A', '_ABC', etc. or '_set'
            results_array = []
        elif inp_token[-2] == '_' and inp_token.lower()[-1] in letters:
            looking_for = inp_token.lower()[-1]
        elif looking_for == None: 
            # NOT a warning: no current means to distinguish between expecting a prior from not expecting one
            return None
            
        for tributary_id in river['tributaries']:
            ltr = letters.popleft()
            if looking_for != None and looking_for != ltr:
                continue
            assert tributary_id in self.psv['reps']
            tributary = self.psv['reps'][tributary_id]
            if 'prevStepResults' not in tributary:
                continue
            prev_results = tributary['prevStepResults']
            # How to match prevResults[inp_token='bwa_bam'] in prevResults to compare to step['inputs'][inp_token='bam_a']
            # The only way is by matching file globs which MUST tie branch leaps together.
            for result_token in prev_results.keys():
                if result_token in self.FILE_GLOBS:
                    # TODO: ultimately these are paired by matching GLOBS!  This frailty might be avoided
                    if self.FILE_GLOBS[result_token] == self.FILE_GLOBS[inp_token]:
                        if expect_set:
                            #print "- Expected set in branch priors and found inp:'"+inp_token + \
                            #                                                "' from prior branch result:'"+result_token+"'" 
                            results_array.append(prev_results[result_token])
                            break  # Done with this branch, but look for more results in next branches
                        if looking_for == ltr:
                            #print "- Expected '"+ltr+"' in branch priors and found inp:'"+inp_token + \
                            #                                                "' from prior branch result:'"+result_token+"'" 
                            return prev_results[result_token] # found result on branch matched by letter
                        #elif looking_for == None:
                        #    return prev_results[result_token] # return first match from first branch
        
        # Went through all branches so return any results found
        if not expect_set or len(results_array) == 0:
            return None
        return results_array


    def find_mystery_param(self,rep,rep_id,step_id,param_name,link_later=False,verbose=False):
        '''Finds parameters that are outputs of previous steps.'''
        # NOT EXPECTED TO OVERRIDE, though only needed for dnaseLaunch.py so far

        # Mystery params come from the non-file outputs of other steps (even from other reps!).
        # If a step definition contains an "output_values" parameter then:
        # - find_all_app_non_file_outputs will look for already run step values and place them in rep['params']
        # - create_or_extend_workflow will create links in rep['prevStepResults'] for steps not yet run 
        # If this step contains "param_links", then each param_spec should say what rep to look for: 'self', 'sister'
        #verbose=True
        
        # Not easy, that is for sure! "param_links": { "target_size": {"name":"reads_filtered", "rep":"sister"} } 
        if "param_links" not in rep["steps"][step_id]:
            return None
        param_links = rep["steps"][step_id]["param_links"]
        if param_name not in param_links:
            return None
        param_spec = param_links[param_name]
        mystery_rep_id = None
        if param_spec["rep"] == "sister": # Hellena's 'brother sistra'?  Catch the reference anyone?
            # How to find "brother sistra"?
            if "sister" not in rep:
                print "ERROR: Looking for mystery parameter but can't find sister of rep '"+rep_id+"'."
                sys.exit(0)
            mystery_rep_id = rep["sister"]
        elif param_spec["rep"] == "self":
            mystery_rep_id = rep_id
        else: # or "parent" or "aunt" ??
            print "ERROR: Looking for mystery parameter but don't (yet) support rep '"+param_spec["rep"]+"'."
            sys.exit(0)
            
        mystery_rep = self.psv["reps"][mystery_rep_id]
        if verbose:
            print "DEBUG: Found rep "+rep['rep_tech']+"'s "+param_spec["rep"]+" as mystery rep '"+mystery_rep['rep_tech']+"'"
        mystery_value = None
        if "prevStepResults" in mystery_rep and param_spec["name"] in mystery_rep["prevStepResults"]:
            mystery_value = mystery_rep["prevStepResults"][param_spec["name"]]  # not yet run result (dx link)
            if verbose:
                print "DEBUG: Found mystery param '"+param_name+"' as '"+param_spec["name"]+\
                      "' from mystery rep '"+mystery_rep['rep_tech']+"' previous results, with value: " + str(mystery_value)
            return mystery_value
        elif "params" in mystery_rep and param_spec["name"] in mystery_rep["params"]:
            mystery_value = mystery_rep["params"][param_spec["name"]] # prior run result (actual string value)
            if verbose:
                print "DEBUG: Found mystery param '"+param_name+"' as '"+param_spec["name"]+\
                      "' from mystery rep '"+mystery_rep['rep_tech']+"' with value: " + str(mystery_value)
            return mystery_value
            
        if not link_later or param_spec["rep"] == "self":  # Linking to self must be to prior step in rep
            return None
            
        # Avoid failure to link BOTH mystery 'sister' params because one is always prior to the other
        # Solution is to put in a known token '@link_later@', then have a final pass after building whole workflow
        # This link_later strategy should only be used if the mystery_step that produces the value can be found.
        mystery_out = None
        for mystery_step_key in mystery_rep["steps"].keys():
            mystery_step = mystery_rep["steps"][mystery_step_key]
            if "output_values" in mystery_step and param_spec["name"] in mystery_step["output_values"]:
                mystery_out = mystery_step["output_values"][param_spec["name"]]
                if verbose:
                    print "DEBUG: Found mystery rep '"+mystery_rep_id+"' and step "+mystery_step_key+"' has '"+ \
                                                                        param_name+" as '"+mystery_out+"'"
                break
        if mystery_out != None:
            mystery_value = '@link_later@'
            if verbose:
                print "DEBUG: Found mystery param '"+param_name+"' as '"+param_spec["name"]+\
                      "' from rep '"+param_spec["rep"]+"' with value: " + str(mystery_value)
            return mystery_value
        
        return None
        

    def wf_find_file_input(self,rep,step_id,file_token,expect_set,verbose=False):
        '''When buidling a workflow, finds a dx defined step file input (which may be for a set) 
           from prior step results, tributary results, or prior files found.'''
        # NOT EXPECTED TO OVERRIDE
        #verbose=(file_token.startswith('bam_biorep'))
        
        # Now try to find prior results or existing files to fill in this input
        file_input = []
        expect_count = 1
        if 'tributaries' in rep:
            # First: check in previous branch (e.g. rep level pipeline) is input to the combining step
            if expect_set and \
               (file_token.endswith('_A') or file_token.endswith('_B') or file_token.endswith('_ABC')):
                expect_count = len(rep['tributaries'])
            prev_input = self.find_results_from_prev_branch(rep, file_token, expect_set)
            if prev_input != None:
                if isinstance(prev_input, list):
                    file_input.extend( prev_input )
                else:
                    file_input.append( prev_input )
            
            if verbose:
                print "DEBUG: Found %d file_inputs from prev_branch" % len(file_input)
            
        if len(file_input) != expect_count and file_token in rep['prevStepResults']:
            if isinstance(rep['prevStepResults'][file_token], list):
                file_input.extend( rep['prevStepResults'][file_token] )
            else:
                file_input.append( rep['prevStepResults'][file_token] )
            if verbose:
                print "DEBUG: Now %d file_inputs after prevStepResults" % len(file_input)

        if len(file_input) == 0:
            return None
        if not expect_set:
            return file_input[0]
        return file_input
            

    def wf_find_param_input(self,rep,rep_id,step_id,param,app_param,verbose=False):
        '''When buidling a workflow, finds a dx defined step param input (sets not yet supported)
           from rep, psv or in rare cases another step''' 
        # NOT EXPECTED TO OVERRIDE
        
        # Now try to find prior results or existing files to fill in this input
        param_input = None
        if param in rep:
            param_input = rep[param]
        #elif "params" in rep and param in rep["params"]: # watch it... this may have cross-fertilizing params
        #    param_input = rep["params"][param]
        elif param in self.psv:
            param_input = self.psv[param]
        else:
            # Rare case of non-file outputs of one step being inputs of another
            param_input = self.find_mystery_param(rep, rep_id, step_id, app_param,link_later=True)
        if param_input == None and not self.template:
            print "ERROR: step '"+step_id+"' unable to locate '"+param+ \
                                              "' in pipeline specific variables (psv)."
            sys.exit(1)
            
        return param_input


    def create_or_extend_workflow(self,rep, rep_id, wf=None,app_proj_id=None):
        '''
        This function will populate or extend a workflow with the steps in rep['stepsToDo'] and return 
        the workflow unlaunched.   It relies on the rep object containing *almost* all information to run
        the PIPELINE_BRANCH associated with the rep.  The rep should contain rep['priors'] defining existing files;
        rep['steps'] objects defining applet, input and output requirements; etc.  A few common parameters will
        be pulled from the psv (pipeline specific variables) object.
        '''
        # NOT EXPECTED TO OVERRIDE

        # Make sure results folder exists first.
        if not self.test:
            if not self.project_has_folder(self.project, rep['resultsFolder']):
                self.project.new_folder(rep['resultsFolder'],parents=True)
        self.build_applets_if_necessary()
                
        if len(rep['stepsToDo']) < 1:
            return None
        if app_proj_id == None:
            app_proj_id = self.proj_id

        # create a workflow object
        if self.multi_rep:
            wf_folder = self.psv['resultsFolder']
            wf_name = self.psv['name']
        else:
            wf_name = rep['name']
            wf_folder = rep['resultsFolder']
        if not self.test and wf == None:
            wf = dxpy.new_dxworkflow(title=wf_name,name=wf_name,folder=wf_folder,
                                           project=self.proj_id,description=self.psv['description'])

        # NOTE: prevStepResults dict contains links to result files to be generated by previous steps
        if 'prevStepResults' not in rep:
            rep['prevStepResults'] = {}
        steps = rep['steps']
        for step in rep['stepsToDo']:
            step_link_later = False
            app_name = steps[step]['app']
            app = self.find_applet_by_name(app_name, app_proj_id)
            inp_defs = app.describe().get('inputSpec') or []
            app_inputs = {}

            # file inputs
            for file_token in steps[step]['inputs'].keys():
                app_inp = steps[step]['inputs'][file_token]
                # Need inp_def to see if set is expected
                for inp_def in inp_defs:
                    if inp_def["name"] == app_inp:
                        break
                if inp_def["name"] != app_inp:
                    print "ERROR: Pipeline definition for applet '"+app_name+"' file_token '"+file_token+"'"
                    print "value of '"+app_inp+"' not found in DX definition."
                    sys.exit(1)
                expect_set = (inp_def["class"] == 'array:file')
                app_inputs[ app_inp ] = self.wf_find_file_input(rep,step,file_token,expect_set)
                if app_inputs[ app_inp ] == None:
                    del app_inputs[ app_inp ] # Because templating will have empty initial input files, this must be done.
                
            # Non-file app param inputs
            if 'params' in steps[step]:
                for param in steps[step]['params'].keys():
                    app_param = steps[step]['params'][param]
                    app_inputs[ app_param ] = self.wf_find_param_input(rep,rep_id,step,param,app_param)
                    # Since a mystery parameter could be from a sister rep, it may require linking later
                    if isinstance(app_inputs[ app_param ],str) and app_inputs[ app_param ] == '@link_later@':
                        step_link_later = True  # Need to save this step's stage for later
                        # Need a final pass of workflow once all branches are built
                        if self.link_later == None:
                            self.link_later = []
                        later_link = { 'rep': rep, 'rep_id': rep_id, 'step_id': step, 'param': app_param }
                        self.link_later.append(later_link)

            # Now we are ready to add wf stage
            if not self.test:
                if self.template:
                    stage_id = wf.add_stage(app, stage_input=app_inputs) # Templates should not fill in a result folder!
                else:
                    stage_id = wf.add_stage(app, stage_input=app_inputs, folder=rep['resultsFolder'])
                if step_link_later:  # At least one parameter will requre linking later so save stage_id
                    steps[step]['stage_id'] = stage_id
                    
            # outputs, which we will need to link to
            for file_token in steps[step]['results'].keys():
                app_out = steps[step]['results'][file_token]
                if self.test:
                    rep['prevStepResults'][ file_token ] = 'fake-for-testing'
                else:
                    rep['prevStepResults'][ file_token ] = dxpy.dxlink({'stage': stage_id, \
                                                                'outputField': app_out })
            # need to add output_values to previous results also (to support param_link requests).
            if "output_values" in steps[step]:
                for out_token in steps[step]["output_values"].keys():
                    app_out = steps[step]["output_values"][out_token]
                    if self.test:
                        rep['prevStepResults'][ out_token ] = 'fake-for-testing'
                    else:
                        rep['prevStepResults'][ out_token ] = dxpy.dxlink({'stage': stage_id, \
                                                                     'outputField': app_out })

        return wf


    def workflow_final_pass(self, wf, app_proj_id=None,verbose=False):
        '''
        Final pass of workflow just in case any '@link_later@' parameters need to be filled in.
        '''
        # NOT EXPECTED TO OVERRIDE
        if self.link_later == None:
            return
        
        # Since non-file outputs from one step might be input params to another step and since
        # those outputs might be in 'sister' branches of the workflow, this final pass may be
        # needed to fill in the links that were not available when the wf was built branch by branch.
        # TODO: This logic could also be used for file results->inputs, but is not currently.
        #verbose=True

        print "Final pass through workflow..."
        
        for later_link in self.link_later:
            rep       = later_link['rep']
            rep_id    = later_link['rep_id']
            step_id   = later_link['step_id']
            app_param = later_link['param']
            if verbose:
                print "DEBUG: looking for '"+rep['rep_tech']+"' step '"+step_id+"' param '"+app_param+"'"
            mystery_param = self.find_mystery_param(rep,rep_id,step_id, app_param)
            if mystery_param == None:
                print "ERROR: Link later for '"+rep['rep_tech']+"' step '"+step_id+"' param '"+app_param+"' failed to find link."
                sys.exit(1)
            # Could also assert that mystery_param is a dx link
            steps = rep['steps']
            step = steps[step_id]
            stage_id = step['stage_id'] # should not be a key error!
            stage = wf.get_stage(stage_id)
            if verbose:
                print "Stage:"
                print stage
            app_inputs = stage['input']
            if verbose:
                print "DEBUG: found '"+app_inputs[app_param]+"' which is being replaced with '"+mystery_param+"'"
            app_inputs[app_param] = mystery_param
            wf.update_stage(stage_id, stage_input=app_inputs)
                
        return


    def report_run_plans(self,run=None):
        '''Report the plans before executing them.'''
        # NOT EXPECTED TO OVERRIDE

        if run == None:
            run = self.psv
        print "Running: "+run['title']
        if 'subTitle' in run:
            print "         "+run['subTitle']

        if not self.no_refs:
            print "- Reference files:"
            # Should be enough to show the ref priors from 'a' rep for single-rep, multi or combined
            for token in self.psv['ref_files']:
                print "  " + self.file_path_from_fid( self.psv['reps']['a']['priors'][token], True)
        
        print "- Results written to: " + self.psv['project'] + ":" +run['resultsFolder']

        print "- Steps to run:"
        to_run_count = 0
        if self.multi_rep:
            for ltr in sorted( self.psv['reps'].keys() ):
                rep = self.psv['reps'][ltr]
                for step in rep['path']:
                    if step in rep['stepsToDo']:
                        print "  * "+rep['rep_tech'] +': '+ rep['steps'][step]['app']+" will be run"
                    else:
                        if not step.find('concat') == 0:
                            print "    "+rep['rep_tech'] +': '+ rep['steps'][step]['app']+" has already been run"
                to_run_count += len(rep['stepsToDo'])
        else:
            to_run_count += len(run['stepsToDo'])
            for step in run['path']:
                if step in run['stepsToDo']:
                    print "  * "+run['steps'][step]['app']+" will be run"
                else:
                    if not step.find('concat') == 0:
                        print "    "+run['steps'][step]['app']+" has already been run"
        if to_run_count == 0:
            print "* All expected results are in the results folder, so there is nothing to do."
            print "  If this experiment/replicate needs to be rerun, then use the --force flag to "
            print "  rerun all steps; or remove suspect results from the folder before relaunching."
            sys.exit(0)


    def workflow_report_and_build(self,run,proj_id=None,verbose=False):
        '''Builds the multi-rep/combined-rep workflow from parts, reporting plans and returning wf.'''
        # NOT EXPECTED TO OVERRIDE
        # Report the plans
        if self.multi_rep:
            self.report_run_plans()
        else:
            self.report_run_plans(run)

        # Build the workflow...
        wf = None
        if self.multi_rep:
            for rep_id in sorted( self.psv['reps'].keys() ):
                rep = self.psv['reps'][rep_id]
                if len(rep['stepsToDo']) > 0: # Going to get noisy with multiple reps
                    if self.test:
                        print "Testing workflow assembly for "+rep['rep_tech']+"..."
                    else:
                        print "Assembling workflow for "+rep['rep_tech']+"..."
                    wf = self.create_or_extend_workflow(rep, rep_id, wf=wf)
        if not self.multi_rep:
            if len(run['stepsToDo']) > 0:
                if self.test:
                    print "Testing workflow assembly for "+run['rep_tech']+"..."
                else:
                    print "Assembling workflow for "+run['rep_tech']+"..."
                wf = self.create_or_extend_workflow(run, None, wf=wf)
        
        # Need a final pass over the whole workflow to patch in params set to '@link_later@'
        if wf != None:  # Only non-test case
            self.workflow_final_pass(wf)
            
        return wf


    def launch_pad(self,wf,run):
        '''Launches or just advertises preassembled workflow.'''
        # NOT EXPECTED TO OVERRIDE
        if wf == None:
            print "ERROR: failure to assemble workflow!"
            sys.exit(1)

        print "Template workflow '" + wf.name + "' has been assembled in "+run['resultsFolder'] + "."

    def run(self):
        '''Runs launch from start to finish using command line arguments.'''
        # NOT EXPECTED TO OVERRIDE

        args = self.get_args()

        print "Retrieving pipeline specifics..."
        self.psv = self.pipeline_specific_vars(args)
        print "Running in project ["+self.proj_name+"]..."
        
        print "Building apps dictionaries..."
        # The point of using 'branches' instead of 'self.PIPELINE_BRANCHES' is that the app dicts *may* be discovered from dx
        branches = {}
        file_globs = {}
        for branch_id in self.PIPELINE_BRANCH_ORDER:
            steps_and_globs = self.assemble_steps_and_globs(branch_id,self.PIPELINE_BRANCHES[branch_id])
            branch = {}
            branch["steps"] = steps_and_globs[0]
            file_globs.update( steps_and_globs[1] )  # Should only need one set of file globs for all pipeline branches/steps
            branches[branch_id] = branch

        # finding fastqs, inputs and prior results in a stadardized way
        for branch_id in self.PIPELINE_BRANCH_ORDER:
            self.find_inputs_and_priors(branch_id,file_globs)

        # finding pipeline specific reference files in a stadardized way
        self.find_all_ref_files()

        # deternine steps to run in a stadardized way
        self.determine_steps_needed(args.force)

        # Preperation is done. Now build up multi-rep workflow
        if self.multi_rep or self.combined_reps:
            run = self.psv
        else:            
            run = self.psv['reps']['a']
        wf = self.workflow_report_and_build(run,self.proj_id)

        # Exit if test only
        if self.test:
            print "TEST ONLY - exiting."
            sys.exit(0)

        # Roll out to pad and possibly launch
        self.launch_pad(wf,run)
                
        print "(success)"

if __name__ == '__main__':
    '''Run from the command line.'''
    # EXPECTED to crash and burn on base class!
    launch = Launch()
    launch.run()

