import os
import sys
import dxpy
import requests
import json
import urlparse

import logging

REFERENCE_FILES = {} ## Dict to cache known Reference Files
FILES = {} ## Dict to cache files
APPLETS = {} ## Dict to cache known applets

KEYFILE = 'keypairs.json'  ## see processkey() Note this file must be in gitignore!
DEFAULT_SERVER = 'https://www.encodeproject.org'
S3_SERVER='s3://encode-files/'

logger = logging.getLogger("DXENCODE")  # not sure this goes here.

def processkey(key):
    ''' check encodedD access keys; assuming the format:
    {
        "default":
                {"server":"https://www.encodeproject.org", "key":"rand_user_name", "secret":"rand_password"},
        "test":
                {"server":"https://test.encodedcc.org", "key":"rand_user_name", "secret":"rand_password"},
        "www":
                {"server":"https://www.encodeproject.org", "key":"rand_user_name", "secret":"rand_password"}
    }
    '''
    if key:
        keysf = open(KEYFILE,'r')
        keys_json_string = keysf.read()
        keysf.close()
        keys = json.loads(keys_json_string)
        key_dict = keys[key]
    else:
        key_dict = {}
    AUTHID = key_dict.get('key')
    AUTHPW = key_dict.get('secret')
    if key:
        SERVER = key_dict.get('server')
    else:
        SERVER = 'https://www.encodeproject.org/'

    if not SERVER.endswith("/"):
        SERVER += "/"

    return (AUTHID,AUTHPW,SERVER)
    ## TODO possibly this should return a dict

def encoded_get(url, AUTHID=None, AUTHPW=None):
    ''' executes GET on Encoded server without without authz '''
    ##TODO possibly add try/except looking for non 4xx?
    HEADERS = {'content-type': 'application/json'}
    if AUTHID and AUTHPW:
        response = requests.get(url, auth=(AUTHID,AUTHPW), headers=HEADERS)
    else:
        response = requests.get(url, headers=HEADERS)
    return response


def project_has_folder(project, folder):
    ''' Checks for a folder in a given DX project '''
    ##TODO Deprecate for find_or_create?
    folders = project.list_folder()['folders']

    return folder in folders


def file_path_from_fid(fid,projectToo=False):
    '''Returns full dx path to file from a file id.'''
    try:
        dxlink = FILES[fid]
    except:
        logger.error("File %s not cached, trying id" % fid)
        dxlink = fid

    fileDict = dxpy.describe(dxlink) # FILES contain dxLinks
    path = fileDict['folder'] + '/' + fileDict['name']
    if projectToo:
        projDict = dxpy.describe(fileDict['project'])
        path = projDict['name'] + ':' + path
    return path


def get_project(projectName, level=None):
    '''Returns the DXProject by name or errors out if not found.'''
    try:
        project = dxpy.find_one_project(name=projectName, name_mode='exact',
                                        level=level, return_handler=False)
    except:
        print "Could not find 1 and only 1 project named '"+projectName+"'."
        sys.exit(1)

    return dxpy.DXProject(project['id'])

def find_or_create_folder(project, sub_folder, root_folder='/'):
    ''' Finds or creates a sub_folder in the specified parent (root) folder'''
    folder = root_folder+sub_folder
    logger.debug("Creating %s (%s)" % (folder, root_folder))
    if folder in project.list_folder(root_folder)['folders']:
        return folder
    else:
        return project.new_folder(folder)

def get_bucket(SERVER, AUTHID, AUTHPW, f_obj):
    ''' returns aws s3 bucket and file name from encodeD file object (f_obj)'''

    #make the URL that will get redirected - get it from the file object's href property
    encode_url = urlparse.urljoin(SERVER,f_obj.get('href'))
    logger.debug(encode_url)

    #stream=True avoids actually downloading the file, but it evaluates the redirection
    r = requests.get(encode_url, auth=(AUTHID,AUTHPW), headers={'content-type': 'application/json'}, allow_redirects=True, stream=True)
    try:
        r.raise_for_status
    except:
        logger.error('%s href does not resolve' %(f_obj.get('accession')))
        sys.exit()
    logger.debug(r)

    #this is the actual S3 https URL after redirection
    s3_url = r.url
    logger.debug(s3_url)

    #release the connection
    r.close()

    #split up the url into components
    o = urlparse.urlparse(s3_url)

    #pull out the filename
    filename = os.path.basename(o.path)

    #hack together the s3 cp url (with the s3 method instead of https)
    return filename, S3_SERVER.rstrip('/') + o.path

def move_files(fids, folder, projectId):
    '''Moves files to supplied folder.  Expected to be in the same project.'''
    for fid in fids:
        try:
            dxlink = FILES[fid]
        except:
            logger.error("File %s not in cache, trying id" % fid)
            dxlink = fid
        fileDict = dxpy.describe(dxlink) # FILES contain dxLinks
        if fileDict['project'] != projectId:
            print "ERROR: Failed to move '" + fileDict['name'] + "' as it is not in '" + \
                                                                                projectId + "'."
            sys.exit(1)
    proj = dxpy.DXProject(projectId)
    if not project_has_folder(proj, folder):
        proj.new_folder(folder,parents=True)
    proj.move(folder,fids)


def find_target_file_set(fileSet,targetFolder,project=None):
    '''Looks for a set of files in a target destination, returning the set.'''
    proj = project
    path = targetFolder
    if targetFolder.find(':') != -1:
        proj, path = targetFolder.split(':')
    if proj.find('project-') == 0:
        projId = proj
    else:
        projId = get_project(proj, level='VIEW').get_id()
    targetFids = []
    for oneFile in fileSet:
        parts = oneFile.rsplit('/',1)
        parts.reverse()
        fid = find_file(path + "/" + parts[0],projId)
        if fid != None:
            targetFids += [ fid ]
    return targetFids


def find_and_copy_read_files(priors, readSet, testOnly, readToken, resultsFolder, projectId=None):
    '''Looks for read files and copies them to results folder if necessary.'''
    readFiles = []
    if len(readSet) > 0:
        readFiles = find_file_set(readSet,projectId)
        if not testOnly:
            priorReads = find_target_file_set(readSet,resultsFolder,projectId)
            if len(priorReads) == len(readFiles):
                readFiles = priorReads
            else: # any files to be moved, move all
                if readToken in priors:
                    del priors[readToken] # make sure to regenerate the combined file if necessary
                readFiles = copy_files(readFiles, projectId, resultsFolder)
        if len(readFiles) > 1:        # If more than 1 file, the 'set' will need the 'concat' step.
            priors[readToken+'_set'] = readFiles
        else:                         # otherwise, the 1 file is same as 'concat' result.
            priors[readToken] = readFiles[0]
    return readFiles


def find_file_set(fileSet,projectId=None):
    '''Find all files in a set, and prints error(s) and exits if any are missing.'''
    fids = []
    if len(fileSet) > 0:
        for oneFile in fileSet:
            fid = find_file(oneFile,projectId,verbose=True)
            if fid != None:
                fids += [ fid ]
        if len(fids) != len(fileSet):
            print "ERROR: expecting " + str(len(fileSet)) + " but only found " + str(len(fids)) + "."
            sys.exit(1) # verbose already gave an error message(s)

    return fids


def copy_files(fids, projectId, folder, overwrite=False):
    '''Copies array of dx file dicts to project:/folder, returning new array of dx file dicts.'''
    newFids = []
    for fid in fids:
        fileDict = dxpy.describe(FILES[fid]) # FILES contain dxLinks
        if fileDict['project'] == projectId:
            # cannot copy into the same project!!!
            # so just leave in place and pretend that we did!
            #proj = dxpy.DXProject(projectId)
            #proj.move(folder,[fid])
            newFids += [ fid ]
            continue

        # Check to see if file already exists.
        alreadyThere = find_file(folder+'/'+fileDict['name'],projectId)
        if alreadyThere is None or overwrite:
            # remove what is alreadyThere?
            #if alreadyThere is not None:
            #    proj = dxpy.DXProject(projectId)
            #    proj.remove_objects([alreadyThere])
            dxFile = dxpy.get_handler(FILES[fid])
            newLink = dxpy.dxlink(dxFile.clone(projectId, folder))
        else:
            newLink = FILES(alreadyThere)
        if newLink == None:
            print "ERROR: Failed in copy of '" + fileDict['project'] + ":" + fileDict['name'] + \
                    "' to '" + projectId + ":" + folder + "'."
            sys.exit(1)
        newDict = dxpy.describe(newLink)
        FILES[newDict['id']] = newLink
        newFids += [ newDict['id'] ]

    return newFids

def resolve_project(project_name, level=None):
    ''' Convert project name into DXProject object '''
    try:
        project = dxpy.find_one_project(name=project_name, name_mode='exact',
                                        level=level, return_handler=False)
    except:
        print 'Could not find 1 and only 1 project named %s; ' % format(project_name)
        sys.exit(1)

    return dxpy.DXProject(project['id'])


def get_file_link(fid, project=None):
    ''' returns a dxlink from cache or directly'''
    if not FILES.get(fid):
        FILES[fid] = dxpy.dxlink(fid,project=None)

    return FILES[fid]

def find_file(filePath,project=None,verbose=False,multiple=False):
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
        projId = get_project(proj, level='VIEW').get_id()
    mode = 'exact'
    if filePath.find('*') or filePath.find('?'):
        mode = 'glob'
    fileDicts = list(dxpy.find_data_objects(classname='file', folder=path, name=fileName,
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
            FILES[fileDict['id']] = dxpy.dxlink(fileDict)
            fids += [ fileDict['id'] ]
        return fids
    else:
        #print "- FOUND '" + proj + ":" + filePath + "'."
        FILES[fileDicts[0]['id']] = dxpy.dxlink(fileDicts[0])
        return fileDicts[0]['id']


def find_reference_file_by_name(reference_name, project_name):
    '''Looks up a reference file by name in the project that holds common tools. From Joe Dale's code.'''
    project = dxpy.find_one_project(name=project_name, name_mode='exact', return_handler=False)
    cached = '*'
    if (reference_name, project['id']) not in REFERENCE_FILES:
        found = dxpy.find_one_data_object(classname="file", name=reference_name,
                                          project=project['id'],
                                          recurse=True,
                                          zero_ok=False, more_ok=False, return_handler=True)
        REFERENCE_FILES[(reference_name, project['id'])] = found
        cached = ''

    logger.debug(cached + "Resolved %s to %s" % (reference_name, REFERENCE_FILES[(reference_name, project['id'])].get_id()))
    return dxpy.dxlink(REFERENCE_FILES[(reference_name, project['id'])])


def find_applet_by_name(applet_name, applets_project_id):
    '''Looks up an applet by name in the project that holds tools.  From Joe Dale's code.'''
    cached = '*'
    if (applet_name, applets_project_id) not in APPLETS:
        found = dxpy.find_one_data_object(classname="applet", name=applet_name,
                                          project=applets_project_id,
                                          zero_ok=False, more_ok=False, return_handler=True)
        APPLETS[(applet_name, applets_project_id)] = found
        cached = ''

    logger.debug(cached + "Resolved %s to %s" % (applet_name, APPLETS[(applet_name, applets_project_id)].get_id()))
    return APPLETS[(applet_name, applets_project_id)]


