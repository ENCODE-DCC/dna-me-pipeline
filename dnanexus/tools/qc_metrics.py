#!/usr/bin/env python2.7
# qc_metrics.py v1.1 Creates a json string of qc_metrics for a given applet.
#                    Write request to stdout and verbose info to stderr.  This allows easy use in dx app scripts.

# imports needed for Settings class:
import os, sys, string, argparse, json

# For a given metric name, expect the following parsing:
EXPECTED_PARSING = {
    "vertical":   {"type": "vertical",   "lines": "", "columns": "", "delimit": None},
    "horizontal": {"type": "horizontal", "lines": "", "columns": "", "delimit": None},
    "singleton":  {"type": "singleton", "delimit": None},
    "edwBamStats":            {"type": "vertical",   "lines": "", "columns": "", "delimit": None},
    "edwComparePeaks":        {"type": "edwComparePeaks"},
    "fastqStatsAndSubsample": {"type": "fastqstats"},
    "samtools_flagstats":     {"type": "flagstats"},
    "samtools_stats":         {"type": "samstats"},
# Implicit:
#    "bismark_extract":        {"type": "bismark_extract"},
#    "bismark_map":            {"type": "bismark_map"},
#    "bismark_lambda":         {"type": "bismark_lambda"},
}

def strip_comments(line,ws_too=False):
    """
    Strips comments from a line (and opptionally leading/trailing whitespace).
    """
    bam = -1
    ix = 0
    while True:
        bam = line[ix:].find('#',bam + 1)
        if bam == -1:
            break
        bam = ix + bam
        if bam == 0:
            return ''
        if line[ bam - 1 ] != '\\':
            line = line[ 0:bam ]
            break  
        else: #if line[ bam - 1 ] == '\\': # ignore '#' and keep looking
            ix = bam + 1
            #line = line[ 0:bam - 1 ] + line[ bam: ]
            
    if ws_too:
        line = line.strip()
    return line 

def string_or_number(a_string):
    try:
        return int(a_string)
    except:
        try:
            return float(a_string)
        except:
            return a_string 

def readline_may_continue(fh):
    """
    Another readLine, but this one supports the '\' continuation character
    """
    line = ''
    while True:
        rawLine = fh.readline()
        if rawLine == '':
            return None
    
        line = line + rawLine.strip()
        if len(line) > 0 and line[ len(line) - 1 ] == '\\':
            #line = line[ 0:len(line) - 1 ]
            line = line[ :-1 ]
            continue
            
        break
        
    return line

def expand_seq(seq_string,one_to_zero=False,verbose=False):
    '''Exapands a string sequence in from of 1,3-6,9,18:21 into 1,3,4,5,6,9,18,19,20,21'''
    seq = []
    if seq_string != None and seq_string != '':    
        seq_parts = seq_string.split(',')
        for val in seq_parts:
            val = val.strip()
            fromto = val.split('-')
            if len(fromto) == 1:
                fromto = val.split(':')
            if len(fromto) == 2:
                if one_to_zero:
                    beg = int(fromto[0]) - 1
                    end = int(fromto[1])
                else:
                    beg = int(fromto[0])
                    end = int(fromto[1]) + 1
                #for n in range(beg,end)
                seq.extend(range(beg,end))
            else:
                if one_to_zero:
                    seq.append(int(val) - 1)
                else:
                    seq.append(int(val))

    if verbose:
        print "seq: "
        print seq
    return seq


def parse_pair(line,columns='',delimit=None,verbose=False):
    '''
    Reads a single line extracting the key-value pair.  Returns a tuple.
    '''

    if line == '':
        return None

    key = ''
    val = ''
    # columns could be '1-3,4' meaning key:1-3 and val:4
    if columns != None and columns != '':
        parts = line.split(delimit)
        col_parts = columns.split(',')
        if len(col_parts) > 0:
            only_cols = expand_seq(col_parts[0],one_to_zero=True,verbose=verbose)
            for col in only_cols:
                if len(parts) > col:
                    if len(key) > 0:
                        key = key + ' '
                    key = key + parts[col].strip()
        if len(col_parts) > 1:
            only_cols = expand_seq(col_parts[1],one_to_zero=True,verbose=verbose)
            for col in only_cols:
                if len(parts) > col:
                    if len(val) > 0:
                        val = val + ' '
                    val = val + parts[col].strip()
    else:
        parts = line.split(delimit,1)
        key = parts[0].strip()
        if len(parts) > 1:
            val = parts[1].strip()
    
    return (key, val)
    
def parse_line(line,columns='',delimit=None,verbose=False):
    '''
    Reads a single line extracting columns and returning the list.
    '''
    if line == '':
        return []

    cols = [] 
    parts = line.split(delimit)
    # columns could be '1,2-3,4-6,9' meaning 2-3 and 4-6 should be merged into single column!
    if columns != None and columns != '':
        col_parts = columns.split(',')
        for col_part in col_parts:
            only_cols = expand_seq(col_part,one_to_zero=True,verbose=verbose)
            col = only_cols[0]
            if len(parts) > col:
                val = parts[col].strip()
                for col in only_cols[1:]:
                    if len(parts) > col:
                        val = val + ' ' + parts[col].strip()
                cols.append(val)
    else:
        for part in parts:
            cols.append(part.strip())
    
    return cols
    
def read_vertical(filePath,lines='',columns='',delimit=None,verbose=False):
    '''
    Reads a file which should contain nothing but 'key value' pairs, one per line. 
    '''
    pairs = {}
    only_lines = expand_seq(lines,verbose)
    line_no = 0
    fh = open(filePath, 'r')
    while True:
        line = readline_may_continue( fh )
        if line == None:
            break
        line_no += 1
        if len(only_lines) > 0 and line_no not in only_lines:
            continue
        if verbose:
            print "["+line+"]"
        line = strip_comments(line,True)
        if line == '':
            continue
        if (line.startswith('#') or line == ''):
            continue
        key, val = parse_pair(line,columns,delimit,verbose)
        pairs[key] = string_or_number(val)
    fh.close()
    
    return pairs
    
def read_horizontal(filePath,lines='',columns='',delimit=None,verbose=False):
    '''
    Reads a file which should contain only two lines: 
    tab separated header line and values line, which are converted to pairs 
    '''
    pairs = {}
    keys = None
    values = None

    only_lines = expand_seq(lines,verbose)
    only_cols = expand_seq(lines,verbose)
    
    fh = open(filePath, 'r')
    line_no = 0
    while True:
        line = readline_may_continue( fh )
        if line == None:
            break
        line_no += 1
        line_no += 1
        if len(only_lines) > 0 and line_no not in only_lines:
            continue
        if verbose:
            print "["+line+"]"
        line = strip_comments(line,True)
        if line == '':
            continue
        if (line.startswith('#') or line == ''):
            continue
        if keys == None:
            keys = parse_line(line,columns,delimit,verbose)
            if verbose:
                print keys
            continue
            
        if values == None:
            values = parse_line(line,columns,delimit,verbose)
            if verbose:
                print values
            break
    fh.close()
    
    for ix, key in enumerate(keys):
        if len(values) > ix:
            pairs[key] = string_or_number(values[ix])
        else:
            pairs[key] = ''

    return pairs
           
def read_singleton(filePath,key,delimit=None,verbose=False):
    '''
    Generic case of single value file. 
    '''
    # TODO support selecting by line and columns!
    pairs = {}

    fh = open(filePath, 'r')
    line = readline_may_continue( fh )
    if line != None:
        if verbose:
            print "["+line+"]"
        line = strip_comments(line,True)
        if line != '':
            values = parse_line(line,delimit=delimit,verbose=verbose)
            pairs[key] = string_or_number(values[0])
        else:
            pairs[key] = ''
    fh.close()
    return pairs

### Now special case routines                
def read_bismark_map(filePath,verbose=False):
    '''
    SPECIAL CASE for bismark map reports. 
    '''
    pairs = {}

    fh = open(filePath, 'r')
    now_lambda = False
    while True:
        line = readline_may_continue( fh )
        if line == None:
            break
        if verbose:
            print "["+line+"]"
        line = strip_comments(line,True)
        if line == '':
            continue
       
        if line == "===== bismark lambda =====":
            now_lambda = True
            continue
            
        parts = line.replace('\t',' ').split(':')
        if len(parts) != 2:
            continue
        # Sequences analysed in total:	556192132
        # Number of alignments with a unique best hit from the different alignments:	442897911
        # Mapping efficiency:	79.6%
        # Sequences with no alignments under any condition:	49874999
        # Sequences did not map uniquely:	63419222
        # Sequences which were discarded because genomic sequence could not be extracted:	71
        # Number of alignments to (merely theoretical) complementary strands being rejected in total:	0
        # Total number of C's analysed:	8120095875
        # Total methylated C's in CpG context:	242388214
        # Total methylated C's in CHG context:	8128428
        # Total methylated C's in CHH context:	28485904
        ## Total methylated C's in Unknown context:
        # Total unmethylated C's in CpG context:	77451246
        # Total unmethylated C's in CHG context:	1698056505
        # Total unmethylated C's in CHH context:	6065585578
        ## Total unmethylated C's in Unknown context
        # C methylated in CpG context:	75.8%
        # C methylated in CHG context:	0.5%
        # C methylated in CHH context:	0.5%
        if parts[0].startswith("Total ") \
        or parts[0].startswith("Sequences ") \
        or parts[0].startswith("Mapping ") \
        or parts[0].startswith("Number of alignments ") \
        or parts[0].startswith("C methylated in C"):
            var = parts[0]
            if now_lambda:
                var = "lambda " + var
            pairs[var] = string_or_number(parts[1].replace(' ',''))

        # CT/CT:	222976469	((converted) top strand)
        # CT/GA:	219921442	((converted) bottom strand)
        # GA/CT:	0	(complementary to (converted) top strand)
        # GA/GA:	0	(complementary to (converted) bottom strand)
        if parts[0] in ["CT/CT","CT/GA","GA/CT","GA/GA"]:
            subparts = parts[1].split()
            var = parts[0] + ' ' + ' '.join(subparts[1:])
            if now_lambda:
                var = "lambda " + var
            pairs[var] = string_or_number(subparts[0].replace(' ','')) 

    fh.close()
    return pairs
    
def read_bismark_combine(filePaths,verbose=False):
    '''
    SPECIAL CASE for combining multiple bismark map reports. 
    '''
    metrics = {}
    files = filePaths.split(',')
    for oneFile in files:
        cur_metrics = read_bismark_map(oneFile,verbose)
        if not metrics:
            metrics = cur_metrics
            continue
            
        # most things are just summed:
        for key in cur_metrics.keys():
            if key.startswith('lambda '):
                var = key[7:]
            else:
                var = key
            if not var.startswith("C methylated in ") and var != "Mapping efficiency":
                if key in metrics.keys():
                    metrics[key] += cur_metrics[key]
                else:
                    metrics[key] = cur_metrics[key]
        
    # Only need to calc percents after all files have been summed and there was more than one file.
    if len(files) > 1 and metrics:          
        # For 3 C contexts the percentages are calculated from the summed numbers
        for context in ["CpG","CHG","CHH","Unknown"]:
            key_me = "Total methylated C's in "+context+" context"
            key_un = "Total unmethylated C's in "+context+" context"
            if key_me in metrics.keys() and key_un in metrics.keys() and metrics[key_me] + metrics[key_un] > 0:
                pcnt = (metrics[key_me] * 100.0) / (metrics[key_me] + metrics[key_un])
                metrics["C methylated in "+context+" context"] = str(round(pcnt,1)) + '%'
            else: #if verbose:
                if key_me not in metrics.keys():
                    print >> sys.stderr, "Couldn't find '"+key_me+"' in metrics." 
                elif key_un not in metrics.keys():
                    print >> sys.stderr, "Couldn't find '"+key_un+"' in metrics."
                else:
                    print >> sys.stderr, "metrics ["+key_me+"] + metrics ["+key_un+"] is zero."
            key_me = "lambda " + key_me
            key_un = "lambda " + key_un
            if key_me in metrics.keys() and key_un in metrics.keys() and metrics[key_me] + metrics[key_un] > 0:
                pcnt = (metrics[key_me] * 100.0) / (metrics[key_me] + metrics[key_un])
                metrics["lambda C methylated in "+context+" context"] = str(round(pcnt,1)) + '%'
            else: #if verbose:
                if key_me not in metrics.keys():
                    print >> sys.stderr, "Couldn't find '"+key_me+"' in metrics." 
                elif key_un not in metrics.keys():
                    print >> sys.stderr, "Couldn't find '"+key_un+"' in metrics."
                else:
                    print >> sys.stderr, "metrics ["+key_me+"] + metrics ["+key_un+"] is zero."
          
        # Mapping efficiency is one off      
        key_tot = "Sequences analysed in total"
        key_map = "Number of alignments with a unique best hit from the different alignments"
        if key_tot in metrics.keys() and key_map in metrics.keys() and metrics[key_tot] > 0:
            pcnt = (metrics[key_map] * 100.0) / (metrics[key_tot])
            metrics["Mapping efficiency"] = str(round(pcnt,1)) + '%'
        else: #if verbose:
            if key_tot not in metrics.keys():
                print >> sys.stderr, "Couldn't find '"+key_tot+"' in metrics." 
            elif key_map not in metrics.keys():
                print >> sys.stderr, "Couldn't find '"+key_map+"' in metrics."
            else:
                print >> sys.stderr, "metrics ["+key_tot+"] is zero."
        key_tot = "lambda " + key_tot
        key_map = "lambda " + key_map
        if key_tot in metrics.keys() and key_map in metrics.keys() and metrics[key_tot] > 0:
            pcnt = (metrics[key_map] * 100.0) / (metrics[key_tot])
            metrics["lambda Mapping efficiency"] = str(round(pcnt,1)) + '%'
        else: #if verbose:
            if key_tot not in metrics.keys():
                print >> sys.stderr, "Couldn't find '"+key_tot+"' in metrics." 
            elif key_map not in metrics.keys():
                print >> sys.stderr, "Couldn't find '"+key_map+"' in metrics."
            else:
                print >> sys.stderr, "metrics ["+key_tot+"] is zero."
        
    return metrics


def read_bismark_split(filePath,verbose=False):
    '''
    SPECIAL CASE for bismark splitting_reports. 
    '''
    pairs = {}


    fh = open(filePath, 'r')
    while True:
        line = readline_may_continue( fh )
        if line == None:
            break
        if verbose:
            print "["+line+"]"
        line = strip_comments(line,True)
        if line == '':
            continue
       
        parts = line.replace('\t',' ').split(':')
        
        #Processed 9260092 lines in total
        if len(parts) == 1 and parts[0].startswith("Processed ") and parts[0].endswith(" lines in total"):
            pairs["Bismark result lines processed"] = string_or_number(line.split()[1])
            continue
        if len(parts) != 2:
            continue
        #Bismark Extractor Version: v0.14.4
        #Bismark result file: single-end (SAM format)
        #Output specified: comprehensive
        #Total number of methylation call strings processed: 9260092
        #Total number of C's analysed:	184118189
        #Total methylated C's in CpG context:	6330842
        #Total methylated C's in CHG context:	1106648
        #Total methylated C's in CHH context:	3837779
        #Total C to T conversions in CpG context:	3361736
        #Total C to T conversions in CHG context:	37848548
        #Total C to T conversions in CHH context:	131632636
        #C methylated in CpG context:	65.3%
        #C methylated in CHG context:	2.8%
        #C methylated in CHH context:	2.8%
        if parts[0].startswith("Total ") \
        or parts[0].startswith("C methylated in ") \
        or parts[0].startswith("Bismark Extractor ") \
        or parts[0].startswith("Bismark result ") \
        or parts[0].startswith("Output "):
            var = parts[0]
            pairs[var] = string_or_number(parts[1].replace(' ',''))

    fh.close()
    return pairs
    
def read_edwComparePeaks(filePath,verbose=False):
    '''
    SPECIAL CASE for edwComparePeaks. 
    '''  
    pairs = read_vertical(filePath,verbose=verbose)
    # Fix one value
    values = pairs["unionSize"].split()
    pairs["unionSize"] = string_or_number(values[1])
    return pairs
    
def read_fastqstats(filePath,verbose=False):
    '''
    SPECIAL CASE for fastqStatsAndSubsample. 
    '''  
    pairs = read_vertical(filePath,verbose=verbose)
    # Fix some values:
    for key in ["aAtPos","cAtPos","gAtPos","tAtPos","nAtPos","qualPos"]:
        values = pairs[key].split(',')
        numbers = []
        for val in values:
            if val != "":
                numbers.append(string_or_number(val))
        pairs[key] = numbers
    return pairs
    
def read_flagstats(filePath,verbose=False):
    '''
    SPECIAL CASE for samtools flagstats. 
    '''
    pairs = {}

    fh = open(filePath, 'r')
    while True:
        line = readline_may_continue( fh )
        if line == None:
            break
        if verbose:
            print "["+line+"]"
        line = strip_comments(line,True)
        if line == '':
            continue
        # 2826233 + 0 in total (QC-passed reads + QC-failed reads)
        if line.find("QC-passed reads") > 0:
        # 2826233 + 0 in total (QC-passed reads + QC-failed reads)
            parts = line.split()
            pairs["total"] = string_or_number(parts[0]) 
            pairs["total_qc_failed"] = string_or_number(parts[2]) 
        # 0 + 0 duplicates
        elif line.find("duplicates") > 0:
            parts = line.split()
            pairs["duplicates"] = string_or_number(parts[0]) 
            pairs["duplicates_qc_failed"] = string_or_number(parts[2]) 
        # 2826233 + 0 mapped (100.00%:-nan%)
        elif "mapped" not in pairs and line.find("mapped") > 0:
            parts = line.split()
            pairs["mapped"] = string_or_number(parts[0]) 
            pairs["mapped_qc_failed"] = string_or_number(parts[2])
            val = parts[4][1:].split(':')[0] 
            pairs["mapped_pct"] = string_or_number(val) 
        # 2142 + 0 paired in sequencing
        elif line.find("paired in sequencing") > 0:
            parts = line.split()
            if int(parts[0]) <= 0: # Not paired-end, so nothing more needed
                break
            pairs["paired"] = string_or_number(parts[0]) 
            pairs["paired_qc_failed"] = string_or_number(parts[2]) 
        # 107149 + 0 read1
        elif line.find("read1") > 0:
            parts = line.split()
            pairs["read1"] = string_or_number(parts[0]) 
            pairs["read1_qc_failed"] = string_or_number(parts[2]) 
        # 107149 + 0 read2
        elif line.find("read2") > 0:
            parts = line.split()
            pairs["read2"] = string_or_number(parts[0]) 
            pairs["read2_qc_failed"] = string_or_number(parts[2]) 
        # 2046 + 0 properly paired (95.48%:-nan%)
        elif line.find("properly paired") > 0:
            parts = line.split()
            pairs["paired_properly"] = string_or_number(parts[0]) 
            pairs["paired_properly_qc_failed"] = string_or_number(parts[2]) 
            val = parts[5][1:].split(':')[0] 
            pairs["paired_properly_pct"] = string_or_number(val) 
        # 0 + 0      singletons (0.00%:-nan%)
        elif line.find("singletons") > 0:
            parts = line.split()
            pairs["singletons"] = string_or_number(parts[0]) 
            pairs["singletons_qc_failed"] = string_or_number(parts[2]) 
            val = parts[4][1:].split(':')[0] 
            pairs["singletons_pct"] = string_or_number(val) 
        # 2046212 + 0 with itself and mate mapped
        elif line.find("with itself and mate mapped") > 0:
            parts = line.split()
            pairs["with_itself"] = string_or_number(parts[0]) 
            pairs["with_itself_qc_failed"] = string_or_number(parts[2]) 
        # 0 + 0 with mate mapped to a different chr (mapQ>=5)
        elif line.find("with mate mapped to a different chr") > 0:
            parts = line.split()
            pairs["diff_chroms"] = string_or_number(parts[0]) 
            pairs["diff_chroms_qc_failed"] = string_or_number(parts[2])
            break

    fh.close()
    return pairs
    
def read_samstats(filePath,verbose=False):
    '''
    SPECIAL CASE of samtools stats 
    '''
    pairs = read_vertical(filePath,delimit=':',verbose=verbose)
    val = pairs['reads MQ0']
    if isinstance(val,str):
        val = val.split('\t')
        pairs['reads MQ0'] = string_or_number(val[0])
    return pairs


def main():
    parser = argparse.ArgumentParser(description =  "Creates a json string of qc_metrics for a given applet. " + \
                                                    "Returns string to stdout and formatted json to stderr.")
    parser.add_argument('-n','--name', required=True,
                        help="Name of metrics in file.")
    parser.add_argument('-f', '--file',
                        help='File containing QC metrics.',
                        required=True)
    parser.add_argument('-t', '--type',
                        help='Type of parsing to be done, if not already understood by name.',
                        choices=['pairs', 'horizontal', 'vertical'],
                        default=None,
                        required=False)
    parser.add_argument('-l', '--lines',
                        help='Only include 1-based numbered lines (e.g. "1,2,5").',
                        default='',
                        required=False)
    parser.add_argument('-c', '--columns',
                        help='Only include 1-based numbered columns (e.g. "1,2,5").',
                        default='',
                        required=False)
    parser.add_argument('-k', '--key',
                        help='Prints just the value for this key.',
                        default=None,
                        required=False)
    parser.add_argument('--keypair',
                        help='Prints the key: value pair for this key.',
                        default=None,
                        required=False)
    parser.add_argument('-d', '--delimit',
                        help='Delimiter to use.',
                        default=None,
                        required=False)
    parser.add_argument('-j', '--json', action="store_true", required=False, default=False, 
                        help="Prints pretty json.")
    parser.add_argument('-v', '--verbose', action="store_true", required=False, default=False, 
                        help="Make some noise.")

    args = parser.parse_args(sys.argv[1:])
    if len(sys.argv) < 3:
        parser.print_usage()
        return
        
    metrics = {}
    
    if args.name in EXPECTED_PARSING:
        parsing = EXPECTED_PARSING[args.name]
    elif args.type != None and args.type in EXPECTED_PARSING: 
        parsing = EXPECTED_PARSING[args.type]
    else: 
        parsing = {"type": args.name }
        
    if args.lines != '':
        parsing["lines"] = args.lines
    if args.columns != '':
        parsing["columns"] = args.columns
    if args.delimit != None:
        parsing["delimit"] = args.delimit
    
    # Read and parse the file into metrics dict
    if parsing["type"] == 'vertical':
        metrics = read_vertical(args.file,parsing["lines"],parsing["columns"],parsing["delimit"],args.verbose)
    elif parsing["type"] == 'horizontal':
        metrics = read_horizontal(args.file,parsing["lines"],parsing["columns"],parsing["delimit"],args.verbose)
    elif parsing["type"] == 'singleton':
        metrics = read_singleton(args.file,args.key,parsing["delimit"],args.verbose)
    
    # Now special case routines                
    elif parsing["type"] == "edwComparePeaks":
        metrics = read_edwComparePeaks(args.file,args.verbose)
    elif parsing["type"] == "bismark_map" or parsing["type"] == "bismark_ref" or parsing["type"] == "bismark_lambda":
        metrics = read_bismark_combine(args.file,args.verbose)
    elif parsing["type"] == "bismark_extract":
        metrics = read_bismark_split(args.file,args.verbose)
    elif parsing["type"] == "fastqstats":
        metrics = read_fastqstats(args.file,args.verbose)
    elif parsing["type"] == 'flagstats':
        metrics = read_flagstats(args.file,args.verbose)
    elif parsing["type"] == 'samstats':
        metrics = read_samstats(args.file,args.verbose)
    else:
        sys.stderr.write('Unknown metric request\n')
        parser.print_usage()
        return

    # Print out the metrics
    if args.key != None and parsing["type"] != 'singleton':
        if args.key in metrics:
            print json.dumps(metrics[args.key])
            sys.stderr.write(json.dumps(metrics[args.key],indent=4) + '\n')
        else:
            print ''   
            sys.stderr.write('(not found)\n')
    elif args.keypair != None:
        if args.keypair in metrics:
            print '"' + args.keypair + '": ' + json.dumps(metrics[args.keypair])
            sys.stderr.write('"' + args.keypair + '": ' + json.dumps(metrics[args.keypair],indent=4) + '\n')
        else:
            print '"' + args.keypair + '": '
            sys.stderr.write('"' + args.keypair + '": \n')
    else: 
        print '"' + args.name + '": ' + json.dumps(metrics,sort_keys=True)
        sys.stderr.write('"' + args.name + '": ' + json.dumps(metrics,indent=4,sort_keys=True) + '\n')
    
if __name__ == '__main__':
    main()

