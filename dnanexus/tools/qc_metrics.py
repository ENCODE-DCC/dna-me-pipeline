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
    "hotspot":                {"type": "hotspot"},
    "samtools_flagstats":     {"type": "flagstats"},
    "samtools_stats":         {"type": "samstats"},
    "phantompeaktools_spp":   {"type": "spp"},
    "pbc":                    {"type": "pbc"}
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
    # columns could be '1-3,4' meaning key:2-3 and val:4
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
    Generic case of single vale file. 
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
    
def read_hotspot(filePath,verbose=False):
    '''
    SPECIAL CASE of read_horizontal that is customized for 'hotspot' output. 
    '''
    pairs = {}
    keys = None
    values = None

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
        #  total tags  hotspot tags    SPOT
        #    2255195       1083552   0.4804
        if keys == None:
            keys = parse_line(line,columns="1:2,3:4,5",delimit=None,verbose=verbose)
            if verbose:
                print keys
            continue
            
        if values == None:
            values = parse_line(line,columns='',delimit=None,verbose=verbose)
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

def read_spp(filePath,verbose=False):
    '''
    SPECIAL CASE for phantompeakqualtools spp. 
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
        
        # strandshift(min): -500 
        # strandshift(step): 5 
        # strandshift(max) 1500 
        # exclusion(min): -500 
        # exclusion(max): -1 
        # FDR threshold: 0.01 
        key, val = parse_pair(line,delimit=':',verbose=verbose)
        if key in ["strandshift(min)","strandshift(step)","exclusion(min)","exclusion(max)","FDR threshold"]:
            pairs[key] = string_or_number(val)
        # strandshift(max) 1500 
        elif line.find("strandshift(max)") == 0: # no colon?
            key, val = parse_pair(line,verbose=verbose)
            pairs[key] = string_or_number(val)
        # done. read 2286836 fragments
        elif line.find("done. read") > 0:
            parts = line.split()
            pairs["read fragments"] = string_or_number(parts[2]) 
        # ChIP data read length 36 
        elif line.find("ChIP data read length") > -1:
            key, val = parse_pair(line,columns='1:4,5',delimit=None,verbose=verbose)
            pairs[key] = string_or_number(val)
       # Minimum cross-correlation value 0.03001068 
        # Minimum cross-correlation shift 1500 
        # Window half size 265 
        # Phantom peak location 40 
        # Phantom peak Correlation 0.08697874 
        elif line.find("Minimum cross-correlation") == 0 \
          or line.find("Window half size") == 0 \
          or line.find("Phantom peak") == 0:
            key, val = parse_pair(line,columns='1:3,4',delimit=None,verbose=verbose)
            pairs[key] = string_or_number(val)
        # Top 3 cross-correlation values 0.086978740217396,0.0864436517524779 
        elif line.find("Top 3 cross-correlation values") == 0:
            key, val = parse_pair(line,columns='1:4,5',delimit=None,verbose=verbose)
            numbers = []
            for num in val.split(','):
                numbers.append(string_or_number(num))
            pairs[key] = numbers 
        # Top 3 estimates for fragment length 40,55 
        elif line.find("Top 3 estimates for fragment length") == 0:
            key, val = parse_pair(line,columns='1:6,7',delimit=None,verbose=verbose)
            numbers = []
            for num in val.split(','):
                numbers.append(string_or_number(num))
            pairs[key] = numbers 
        # Normalized Strand cross-correlation coefficient (NSC) 2.898259 
        # Relative Strand cross-correlation Coefficient (RSC) 1 
        elif line.find("Strand cross-correlation") > 0:
            key, val = parse_pair(line,columns='1:5,6',delimit=None,verbose=verbose)
            pairs[key] = string_or_number(val)
        # Phantom Peak Quality Tag 1 
        elif line.find("Phantom Peak Quality Tag") == 0:
            key, val = parse_pair(line,columns='1:4,5',delimit=None,verbose=verbose)
            pairs[key] = string_or_number(val)

    fh.close()
    return pairs

def read_pbc(filePath,verbose=False):
    '''
    SPECIAL CASE for pbc. 
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
        
        # 2286836	2219898	2176268	37175	0.970729	0.980346	58.541170
        parts = line.split()
        if len(parts) > 0:
            pairs["reads"] = string_or_number(parts[0])
        if len(parts) > 1:
            pairs["locations"] = string_or_number(parts[1])
        if len(parts) > 2:
            pairs["mapped by 1 read"] = string_or_number(parts[2])
        if len(parts) > 3:
            pairs["mapped by 2 reads"] = string_or_number(parts[3])
        if len(parts) > 4:
            pairs["locations per read"] = string_or_number(parts[4])
        if len(parts) > 5:
            pairs["proportion of 1 read locations"] = string_or_number(parts[5])
        if len(parts) > 6:
            pairs["ratio: 1 read over 2 read locations"] = string_or_number(parts[6])
        break
    fh.close()
    return pairs
 
            
def main():
    parser = argparse.ArgumentParser(description =  "Creates a json string of qc_metrics for a given applet. " + \
                                                    "Returns string to stdout and formatted json to stderr.")
    parser.add_argument('-n','--name', required=True,
                        help="Name of metrics in file.")
    parser.add_argument('-f', '--file',
                        help='File containing QC metrics.',
                        required=True)
    #parser.add_argument('-t', '--type',
    #                    help='Type of parsing to be done.',
    #                    choices=['pairs', 'horizontal'],
    #                    required=True)
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
    else: 
        parsing = EXPECTED_PARSING["vertical"]
        
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
    elif parsing["type"] == "fastqstats":
        metrics = read_fastqstats(args.file,args.verbose)
    elif parsing["type"] == 'flagstats':
        metrics = read_flagstats(args.file,args.verbose)
    elif parsing["type"] == 'hotspot':
        metrics = read_hotspot(args.file,args.verbose)
    elif parsing["type"] == 'samstats':
        metrics = read_samstats(args.file,args.verbose)
    elif parsing["type"] == 'spp':
        metrics = read_spp(args.file,args.verbose)
    elif parsing["type"] == 'pbc':
        metrics = read_pbc(args.file,args.verbose)

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
        print '"' + args.name + '": ' + json.dumps(metrics)
        sys.stderr.write('"' + args.name + '": ' + json.dumps(metrics,indent=4) + '\n')
    
if __name__ == '__main__':
    main()

