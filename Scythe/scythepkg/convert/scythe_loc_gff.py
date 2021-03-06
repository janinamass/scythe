#!/usr/bin/env python
import re
import sys
import getopt

"""
scythe_loc_gff.py
~~~
Convert gff3 format to loc format.
"""


def usage():
    """Print usage."""
    print("""
    #####################################
    #  scythe_loc_gff.py  -f FILE.gff3  #
    #####################################

    -f, --file=gff3_FILE
    -o, --output=FILE        output file [default: gff3_FILE.loc]
    [ -s, --sequences ]          fasta file with sequences,
                             in case "longest" sequences are
                             not labelled in the gff

    -h, --help               prints this
    -H, --HELP               show help on format
    """)
    sys.exit(2)


def formatHelp():
    print("""
    #------------ loc output format ---------------------------#

    LOCUS0\tTRANSCRIPT0_0\tTRANSCRIPT0_1\t...\tTRANSCRIPT0_n
    LOCUS1\tTRANSCRIPT1_0\t...\tTRANSCRIPT1_m
    .
    .
    .
    LOCUSk\tTRANSCRIPTk_0\t...TRANSCRIPTk_l

    #----------------------------------------------------------#


    #-------- gff version 3 input format ----------------------#
    example (tab is shown as \\t):

    ##gff-version 3
    L1i\\texample\\tgene\\t1000\\t9000\\t.\\t+\\t.\\tID=L1;Name=L1;Note=example
    L1.a\\texample\\tmRNA\\t1000\\t9000\\t.\\t+\\t.\\tID=L1.a;Parent=L1;Name=L1.a


    Note that this script only relies on the "ID" and "Parent" tags,
    "Name" will be ignored. If "longest" is specified (eg phytozome)
    and =="1", the transcript will be placed on the first position
    for its gene.
    #----------------------------------------------------------#

    """)


def checkGff3(infile):
    with open(infile, 'r') as f:
        l = f.readline().strip()
        if l == "##gff-version 3":
            return True
        else:
            return False


def read_gff_attrib(infile, attrib):
    with open(infile, 'r') as f:
        for l in f:
            if l.strip() == "":
                pass
            elif l.startswith("#"):
                pass
            else:
                tmp = l.split("\t")[2].lower().strip()
                if tmp == attrib.lower():
                    ln = l.split("\t")[8].strip()
                    yield (ln)


def read_gene_mrna(infile):
    genes = {}
    longest = set()
    p_id = r"""(ID=)(.*?)(;|$)"""
    p_parent = r"""(Parent=)(.*?)(;|$)"""
    o_longest = r"""(longest=)(.*?)(;|$)"""
    for g in read_gff_attrib(infile, "gene"):
        m = re.findall(p_id, g)
        if m:
            mid = m[0][1]
            genes[mid] = []
    for g in read_gff_attrib(infile, "mrna"):
        m = re.findall(p_id, g)
        if m:
            mid = m[0][1]
        m = re.findall(p_parent, g)
        if m:
            # phytozome tags their representive transcript as 'longest'
            # add 'longest' as first element
            if m[0][1] in genes:
                l = re.findall(o_longest, g)
                if l:
                    if l[0][1] == "1":
                        longest.add(mid)
                        genes[m[0][1]].insert(0, mid)
                    else:
                        genes[m[0][1]].append(mid)
                else:
                    genes[m[0][1]].append(mid)
            else:
                raise Warning("{}: parent is not in known genes".format(mid))

    return genes.copy()


def writeLoc(genesDct, outfile):
    f = lambda x: "\t".join(genesDct[x])
    tmp = [g + "\t" + f(g) for g in sorted(genesDct.keys())]
    res = "\n".join(tmp)
    with open(outfile, 'w') as out:
        out.write(res)


def read_gff2loc(infile, outfile):
    infile = open(infile, 'r')
    outfile = open(outfile, 'w')
    loci = {}
    longest = {}
    rawstr = r"""(Name=)(.*);pacid.*(longest=)(.*);(Parent=)(.*)"""
    for ln in infile:
        s = ln
        m = re.findall(rawstr, s)
        if len(m) > 0:
            name = m[0][1]
            isLongest = m[0][3]
            parent = m[0][5]
            if isLongest == str(1):
                if parent in longest:
                    print("#Warning " + parent + " has more than one default model\nCheck your gff -> ",
                          longest[parent], name)
                longest[parent] = name  # longest will be printed to 2nd col
            elif isLongest == str(0):
                if parent in loci:
                    loci[parent].append(name)
                else:
                    loci[parent] = [name]

    s_def = sorted(longest.keys())
    for k_def in s_def:
        try:
            outfile.write(k_def + "\t" + longest[k_def] + "\t" + "\t".join(loci[k_def]) + "\n")
        except KeyError:
            outfile.write(k_def + "\t" + longest[k_def] + "\n")
        if k_def in loci:
            del loci[k_def]
    s = sorted(loci.keys())
    for k in s:
        try:
            outfile.write(k + "\t" + longest[k] + "\t" + "\t".join(loci[k]) + "\n")
        except KeyError:
            print("#Warning " + k + " has no default model\n")
            outfile.write(k + "\t" + "\t".join(loci[k]) + "\n")
    return loci


def getSequenceLengths(fasta=None):
    faDct = {}
    if not fasta:
        return None
    else:
        for header, sequence in FastaParser().read_fasta(fasta," ",0):
            faDct[header] = len(sequence)
    return faDct.copy()



class FastaParser(object):
    def read_fasta(self, fasta, delim = None, asID = 0):
        """read from fasta fasta file 'fasta'
        and split sequence id at 'delim' (if set)\n
        example:\n
        >idpart1|idpart2\n
        ATGTGA\n
        and 'delim="|"' returns ("idpart1", "ATGTGA")
        """
        name = ""
        fasta = open(fasta, "r")
        while True:
            line = name or fasta.readline()
            if not line:
                break
            seq = []
            while True:
                name = fasta.readline()
                name = name.rstrip()
                if not name or name.startswith(">"):
                    break
                else:
                    seq.append(name)
            joinedSeq = "".join(seq)
            line = line[1:]
            if delim:
                #vprint("delim set to "+delim)
                line = line.split(delim)[asID]
            yield (line.rstrip(), joinedSeq.rstrip())
        fasta.close()



def main():
    outfile = None
    infile = None
    fasta = None
    fastaDct = None
    try:
        opts, args = getopt.gnu_getopt(sys.argv[1:], "f:s:hHo:", ["file=","sequences=", "help", "HELP", "output="])
    except getopt.GetoptError as err:
        print(str(err))
        usage()
    for o, a in opts:
        if o in ("-f", "--file"):
            infile = a
        elif o in ("-h", "--help"):
            usage()
        elif o in ("-H", "--HELP"):
            formatHelp()
        elif o in ("-o", "--output"):
            outfile = a
        elif o in ("-s", "--sequences"):
            fasta = a
        else:
            assert False, "unhandled option"

    if infile is None:
        usage()
    if outfile is None:
        outfile = infile + ".loc"
    if not checkGff3(infile):
        raise Warning("Gff3 line (##gff-version 3) missing at beginning of file. This might not work.")
    if fasta:
        fastaDct = getSequenceLengths(fasta)

    genes = read_gene_mrna(infile)
    if fastaDct:
        f = lambda x: [fastaDct[y] for y in x]
        for k in genes:
            l = genes[k]
            try:
                mxindex = f(l).index(max(f(l)))
                mxel = l[mxindex]
                if len(genes[k]) >1:
                    genes[k].remove(mxel)
                    genes[k].insert(0, mxel)
            except KeyError:
                sys.stderr.write("{} not in fasta\n".format(k))



    writeLoc(genes, outfile)

if __name__ == "__main__":
    main()
