#!/usr/bin/env python
# coding: utf-8
# test-py-parallel 0.0.1
# Generated by dx-app-wizard.
#
# Scatter-process-gather execution pattern: Your app will split its
# input into multiple pieces, each of which will be processed in
# parallel, after which they are gathered together in some final
# output.
#
# This pattern is very similar to the "parallelized" template.  What
# it does differently is that it formally breaks out the "scatter"
# phase as a separate black-box entry point in the app.  (As a side
# effect, this requires a "map" entry point to call "process" on each
# of the results from the "scatter" phase.)
#
# Note that you can also replace any entry point in this execution
# pattern with an API call to run a separate app or applet.
#
# The following is a Unicode art picture of the flow of execution.
# Each box is an entry point, and vertical lines indicate that the
# entry point connected at the top of the line calls the entry point
# connected at the bottom of the line.  The letters represent the
# different stages in which the input is transformed, e.g. the output
# of the "scatter" entry point ("array:B") is given to the "map" entry
# point as input.  The "map" entry point calls as many "process" entry
# points as there are elements in its array input and gathers the
# results in its array output.
#
#          ┌──────┐
#       A->│ main │->D (output from "postprocess")
#          └┬─┬─┬─┘
#           │ │ │
#          ┌┴──────┐
#       A->│scatter│->array:B
#          └───────┘
#             │ │
#            ┌┴──────────────┐
#   array:B->│      map      │->array:C
#            └─────────┬─┬─┬─┘
#               │      │ . .
#               │     ┌┴──────┐
#               │  B->│process│->C
#               │     └───────┘
#            ┌──┴────────┐
#   array:C->│postprocess│->D
#            └───────────┘
#
# A = original app input, split up by "scatter" into pieces of type B
# B = an input that will be provided to a "process" entry point
# C = the output of a "process" entry point
# D = app output aggregated from the outputs of the "process" entry points
#
# See https://wiki.dnanexus.com/Developer-Portal for documentation and
# tutorials on how to modify this file.
#
# DNAnexus Python Bindings (dxpy) documentation:
#   http://autodoc.dnanexus.com/bindings/python/current/

import os
import dxpy
import subprocess
import shlex
import glob
import logging
import json
import re

DEBUG = True

logger = logging.getLogger(__name__)
logger.addHandler(dxpy.DXLogHandler())
logger.propagate = False

if DEBUG:
    logger.setLevel(logging.DEBUG)
else:
    logger.setLevel(logging.INFO)


STRIP_EXTENSIONS = ['.gz', '.fq', '.fastq', '.fa', '.fasta']
ALIGN_SCRIPT = '/usr/bin/dname_align_se.sh'
QC_SCRIPT = '/usr/bin/qc_metrics.py'
VERSION_SCRIPT = '/usr/bin/tool_versions.py'
PROPERTY_SCRIPT = '/usr/bin/parse_property.py'


def strip_extensions(filename, extensions):
    basename = filename
    for extension in extensions:
        basename = basename.rpartition(extension)[0] or basename
    return basename

def merge_bams(bam_files, bam_root, use_cat, use_sort, nthreads):

    fnames = []
    for bam in bam_files:
        dxbam = dxpy.DXFile(bam)
        dxfn = dxbam.describe()['name']
        logger.info("* Downloading %s... *" % dxfn)
        dxpy.download_dxfile(bam, dxfn)
        fnames.append(dxfn)

    outfile_name = bam_root
    logger.info("* Merged alignments file will be: %s *" % outfile_name + '.bam')
    if len(fnames) == 1:
        # UNTESTED 
        rep_outfile_name = bam_root + '_bismark_biorep'
        logger.info("* Only one input file (%s), no merging required." % fnames[0])
        os.rename(fnames[0], outfile_name + '.bam')

    else:
        if use_cat:
            for fn in fnames:
                if not os.path.isfile('sofar.bam'):
                    os.rename(fn, 'sofar.bam')
                else:
                    logger.info("* Merging...")
                    # NOTE: keeps the first header
                    catout = subprocess.check_output(['samtools', 'cat', 'sofar.bam', fn, '>', 'merging.bam'])
                    logger.info(catout)
                    os.rename('merging.bam', 'sofar.bam')

            # At this point there is a 'sofar.bam' with one or more input bams

            logger.info("* Files merged into %s (via cat) *" % outfile_name + '.bam')

        else:
            # use samtools merge
            # UNTESTED
            filelist = " ".join(fnames)
            logger.info("Merging via merge %s " % filelist)
            mergeout = subprocess.check_output(['samtools', 'merge', '-nf', 'sofar.bam'] + fnames)
            # this gets renamed later
            logger.info(mergeout)

        if use_sort:
            # sorting needed due to samtools cat
            # UNTESTED
            logger.info("* Sorting merged bam...")
            sortout = subprocess.check_output(['samtools', 'sort', '-@', nthreads, '-m', '6G', '-f,' 'sofar.bam', 'sorted.bam'])
            logger.info(sortout)
            os.rename('sorted.bam', outfile_name + '.bam')
        else:
            os.rename('sofar.bam', outfile_name + '.bam')

    return outfile_name + '.bam'


def merge_reports(outfile_name, report_files, bam_root):
    out = ["### Combined Bismark map report for split fastq ### "]
    logger.info(out[0])
    report_file_names = []
    for report in report_files:
        dxreport = dxpy.DXFile(report)
        rfn = dxreport.describe()['name']
        logger.info("* Downloading %s... *" % rfn)
        dxpy.download_dxfile(report, rfn)
        report_file_names.append(rfn)
        out.append("###################################")
        out.append("### Map report for %s ###" % rfn)
        out = out + open(rfn, 'r').readlines()

    outfh = open(outfile_name+'_map_report.txt', 'w')
    outfh.write("\n".join(out))
    return (outfile_name+'_map_report.txt', report_file_names)


def merge_qc(outfile_name, report_files):

    qc_stats = ''
    if os.path.isfile(QC_SCRIPT):
        qc_stats = json.loads('{'+subprocess.check_output(['qc_metrics.py', '-n', 'bismark_map', '-f'] + [','.join(report_files)])+'}')

    logger.debug("** merge_qc: %s " % outfile_name)
    logger.info("* Collect bam stats...")
    subprocess.check_call(['samtools', 'index', outfile_name+'.bam'])
    subprocess.check_call(['samtools', 'flagstat', outfile_name+'.bam', '>', outfile_name + '_flagstat.txt'])
    subprocess.check_call(['samtools', 'stats', outfile_name+'.bam', '>', outfile_name + '_samstats.txt'])
    logger.info(subprocess.check_output(['head', '-3', outfile_name + '_samstats.txt']))
    logger.info(subprocess.check_ouptput(['grep', '^SN', outfile_name + '_samstats.txt', '|', 'cut', '-f', '2-', '>', outfile_name + ' _samstats_summary.txt']))

    logger.info("* Prepare metadata...")
    reads = 0
    read_len = 0
    if os.path.isfile(QC_SCRIPT):
        meta = subprocess.check_output(['qc_metrics.py', '-n', 'samtools_flagstats', '-f', outfile_name+'_flagstat.txt'])
        qc_stats.extend(json.loads('{'+meta+'}'))
        reads = subprocess.check_output(['qc_metrics.py', '-n', 'samtools_flagstats', '-f', outfile_name+'_flagstat.txt', '-k', 'total'])
        meta = subprocess.check_output(['qc_metrics.py', '-n', 'samtools_stats', '-d', ':', '-f', outfile_name+'_samstats.txt'])
        qc_stats.extend(json.loads('{'+meta+'}'))
        read_len = subprocess.check_output(['qc_metrics.py', '-n', 'samtools_stats',  '-d', ':', '-f', outfile_name+'_samstats_summary.txt', '-k', 'average length'])

    logger.info(json.dumps(qc_stats))
    # All qc to one file per target file:
    qc_file = outfile_name + '_qc.txt'
    fh = open(qc_file, 'w')
    fh.write("===== samtools flagstat =====\n")

    subprocess.check_call(['cat', outfile_name + '_flagstat.txt', '>>', qc_file])
    fh.write("===== samtools stats =====\n")
    subprocess.check_call(['cat', outfile_name + '_samstats.txt', '>>', qc_file])

    fh.close()

    return (qc_file, reads, read_len, json.dumps(qc_stats))


@dxpy.entry_point("postprocess")
def postprocess(bam_files, report_files, bam_root, nthreads=8, use_cat=False, use_sort=False):
    # This is the "gather" phase which aggregates and performs any
    # additional computation after the "map" (and therefore after all
    # the "process") jobs are done.

    if DEBUG:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)

    logger.debug("** In Postprocess - refactoed dme-merge-bams - *")

    if os.path.isfile(VERSION_SCRIPT):
        versions = subprocess.check_call(['tool_versions.py', '--dxjson', 'dnanexus-executable.json'])


    merged_bam = merge_bams(bam_files, bam_root, use_cat, use_sort, nthreads)

    (merged_report, report_file_names) = merge_reports(bam_root, report_files, bam_root)

    (merged_qc, nreads, read_length, metadata) = merge_qc(bam_root, report_file_names)

    props = {
        'SW': versions,
        'reads': nreads,
        'read_length': read_length
    }
    output = {
        "bam_techrep": dxpy.dxlink(dxpy.upload_local_file(merged_bam, details=metadata, properties=props)),
        "bam_techrep_qc": dxpy.dxlink(dxpy.upload_local_file(merged_qc), detais=metadata, properties={'SW': versions}),
        "map_techrep": dxpy.dxlink(dxpy.upload_local_file(merged_report), detais=metadata, properties={'SW': versions}),
        "reads": nreads,
        "metadata": metadata
    }
    return output

@dxpy.entry_point("process")
def process(scattered_input, dme_ix, ncpus, reads_root):
    # Fill in code here to process the input and create output.

    if DEBUG:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)

    dme_ix = dxpy.DXFile(dme_ix)

    # The following line(s) download your file inputs to the local file system
    # using variable names for the filenames.

    dxpy.download_dxfile(dme_ix.get_id(), "index.tgz")
    fq = dxpy.DXFile(scattered_input)
    name = fq.describe()['name']
    dxpy.download_dxfile(fq.get_id(), name)
    bam_root = name + '_techrep'

    logger.info("* === Calling DNAnexus and ENCODE independent script... ===")
    logger.debug("** DIR: %s" % os.listdir('./'))
    logger.debug(subprocess.check_output(['head', name]))
    if os.path.isfile(ALIGN_SCRIPT):
        logger.debug("** Executable %s exists" % ALIGN_SCRIPT)
    else:
        logger.debug("** Executable %s DOES NOT exist" % ALIGN_SCRIPT)
        exit(1)
    logger.debug('** command line: %s index.tgz %s %s %s' % (ALIGN_SCRIPT, name, ncpus, bam_root))
    map_out = subprocess.check_output([ALIGN_SCRIPT, 'index.tgz', name, str(ncpus), bam_root, 'no_stats'])
    logger.info("* === Returned from dname_align_se  ===")

    # As always, you can choose not to return output if the
    # "postprocess" stage does not require any input, e.g. rows have
    # been added to a GTable that has been created in advance.  Just
    # make sure that the "postprocess" job does not run until all
    # "process" jobs have finished by making it wait for "map" to
    # finish using the depends_on argument (this is already done for
    # you in the invocation of the "postprocess" job in "main").

    logger.debug("** DIR: %s" % os.listdir('./'))
    logger.debug("** OUTPUT DIR: %s" % os.listdir('output/'))

    os.rename(bam_root+'_bismark.bam', bam_root+'.bam')
    return {
        "bam_file": dxpy.dxlink(dxpy.upload_local_file(bam_root+'.bam')),
        "report_file": dxpy.dxlink(dxpy.upload_local_file(bam_root+'_bismark_map_report.txt'))
    }


@dxpy.entry_point("map")
def map_entry_point(array_of_scattered_input, process_input):
    # The following calls "process" for each of the items in
    # *array_of_scattered_input*, using as input the item in the
    # array, as well as the rest of the fields in *process_input*.
    if DEBUG:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)

    logger.debug("** in map entry point with %s *" % process_input)
    process_jobs = []
    for item in array_of_scattered_input:
        logger.debug("** scattering: %s *" % item)
        process_input["scattered_input"] = item
        process_jobs.append(dxpy.new_dxjob(fn_input=process_input, fn_name="process"))

    logger.info("* %s scatter jobs started *" % len(array_of_scattered_input))
    bams = []
    reports = []
    for subjob in process_jobs:
        bams.append(subjob.get_output_ref('bam_file'))
        reports.append(subjob.get_output_ref('report_file'))

    return {
        "bam_files": bams,
        "report_files": reports,
    }


@dxpy.entry_point("scatter")
def scatter(orig_reads, split_size):
    # Fill in code here to do whatever is necessary to scatter the
    # input.
    if DEBUG:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)

    splitsize = split_size * 1000000 * 4
    # each FQ read is 4 lines
    os.mkdir('splits')

    for f in orig_reads:
        reads_filename = dxpy.describe(f)['name']
        reads_basename = strip_extensions(reads_filename, STRIP_EXTENSIONS)
        dxpy.download_dxfile(dxpy.DXFile(f).get_id(), reads_filename)

        reads_root_name = simplify_name() or reads_basename

        logger.info('* RUNNING /bin/zcat %s | /usr/bin/split -l %d -d - %s ' % (reads_filename, splitsize, 'splits/' + reads_root_name))
        split_out = subprocess.check_output('/bin/zcat %s | /usr/bin/split -l %d -d - %s ' % (reads_filename, splitsize, 'splits/' + reads_root_name), shell=True)

    logger.info(split_out)
    splits = os.listdir('splits')
    logger.info("* Return from scatter: %s *" % splits)

    # SHould we gzip here?
    return {
        "array_of_scattered_input": [ 
            dxpy.dxlink(dxpy.upload_local_file('splits/' + split_file)) for split_file in splits]
        }


def simplify_name():
    # Try to simplify the names

    rep_root = ''
    if os.path.isfile(PROPERTY_SCRIPT):
        logger.debug(" ".join(["Simplify Name:", 'parse_property.py', '--job', os.environ['DX_JOB_ID'],'--root_name', '--quiet']))
        rep_root = subprocess.check_output(['parse_property.py', '--job', os.environ['DX_JOB_ID'],'--root_name', '--quiet'])
        rep_root.strip()
        if not re.match('\w', rep_root):
            return ''
        logger.debug("** Simplified Name: %s", rep_root)
    return rep_root


@dxpy.entry_point("main")
def main(reads, dme_ix, ncpus, splitsize):

    # The following line(s) initialize your data object inputs on the platform
    # into dxpy.DXDataObject instances that you can start using immediately.

    #dx_reads = [dxpy.DXFile(item) for item in reads]

    # The following line(s) download your file inputs to the local file system
    # using variable names for the filenames.



    # We first create the "scatter" job which will scatter some input
    # (replace with your own input as necessary).
    logger.info("* Start Scatter with %d files %sM read splits *" % (len(reads), splitsize))

    scatter_job = dxpy.new_dxjob(fn_input={
                                 'orig_reads': reads,
                                 'split_size': splitsize,
                                 },
                                 fn_name="scatter")

    # We will want to call "process" on each output of "scatter", so
    # we call the "map" entry point to do so.  We can also provide
    # here additional input that we want each "process" entry point to
    # receive, e.g. a GTable ID to which the "process" function should
    # add rows of data.

    reads_root = simplify_name() or strip_extensions(dxpy.describe(reads[0])['name'], STRIP_EXTENSIONS)

    map_input = {
        "array_of_scattered_input": scatter_job.get_output_ref("array_of_scattered_input"),
        "process_input": {
            "reads_root": reads_root,
            "ncpus": ncpus,
            "dme_ix": dme_ix
            }
        }
    logger.info("* Start Map with: %s *" % map_input)
    map_job = dxpy.new_dxjob(fn_input=map_input, fn_name="map")

    # Finally, we want the "postprocess" job to run after "map" is
    # done calling "process" on each of its inputs.  Note that a job
    # is marked as "done" only after all of its child jobs are also
    # marked "done".
    logger.info("* Waiting for map job to finish...")
    postprocess_input = {
        "bam_files": map_job.get_output_ref("bam_files"),
        "report_files": map_job.get_output_ref("report_files"),
        "bam_root": reads_root + '_techrep'
        }
    logger.info("* Start Post process with: %s *" % postprocess_input)
    postprocess_job = dxpy.new_dxjob(fn_input=postprocess_input,
                                     fn_name="postprocess",
                                     depends_on=[map_job])

    # If you would like to include any of the output fields from the
    # postprocess_job as the output of your app, you should return it
    # here using a job-based object reference.
    #
    # return { "app_output_field": postprocess_job.get_output_ref("final_output"), ...}
    #
    # Tip: you can include in your output at this point any open
    # objects (such as gtables) which will be closed by a job that
    # finishes later.  The system will check to make sure that the
    # output object is closed and will attempt to clone it out as
    # output into the parent container only after all subjobs have
    # finished.

    output = {}
    output["bam_techrep"] = dxpy.dxlink(postprocess_job.get_output_ref("bam_techrep"))
    output["bam_techrep_qc"] = dxpy.dxlink(postprocess_job.get_output_ref("bam_techrep_qc"))
    output["map_techrep"] = dxpy.dxlink(postprocess_job.get_output_ref("map_techrep"))
    output["reads"] = postprocess_job.get_output_ref("reads")
    output["metadata"] = postprocess_job.get_output_ref("metadata")

    return output

dxpy.run()
