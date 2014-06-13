#!/usr/env/ python
import re
import sys
import getopt

def usage():
    print ("""
    ################################
    #  Scythe_proteinOrtho2grp.py  #
    ################################

    general options:
    -f, --file=INFILE        (output of OrthoMCL)
    -l, --loc=LOCFILE        concatenated .loc file or comma separated list of .loc files
    -o, --output=OUTFILE     output file [default: INFILE.grp]
    -p, --prefix=PREFIX      eg 'gr_' of 'gr_0:'
    -h, --help               prints this

    other:
    -d, --delimiter=SEP      in case your ids have an additional prefix,
                             eg set -d"|" and -a1 if your file has identifiers like SPXY|genemodelID
    -i, --asID=NUM


    ------------
    .grp format: GroupID\tgeneIDiSp1\tgeneIDjSp2\t...geneIDkSpn
    """)
    sys.exit(2)

def read_omcl2grp(infile, locfile, outfile, prefix, delim, asID):
    mod2loc = {}
    locfiles = locfile.split(",")
    for locfile in locfiles:
        locfile = locfile.strip()
        locfile = open(locfile, 'r')
        for ln in locfile:
            ln = ln.rstrip()
            g, gm = ln.split("\t")[0],ln.split("\t")[1:]
            for m in gm:
                mod2loc[m]=g
        locfile.close()

    infile = open(infile, 'r')
    outfile = open(outfile, 'w')
    cnt=0
    omit = 0
    for ln in infile:
        ln = ln.strip()
        grp, ln = ln.split(":")[0],ln.split(":")[1]
        grp = grp.split("(")[0]
        grp = grp.split(prefix)[1]
        ln = ln.strip().split(" ")
        spec = [tmp.split("(")[1] for tmp in ln]
        gm = [tmp.split("(")[0] for tmp in ln]
        tmp = set()
        for s in spec:
            if s not in tmp:
                tmp.add(s)
            else:
                ln=""
                gm = []
                omit +=1
                continue #paralogs
        allid = []
        for g in gm:
            if delim and asID:
                g =  g.split(delim)[asID]
            allid.append(mod2loc[g])
        if len(allid)>0:
            cnt +=1
            outfile.write(grp+"\t"+"\t".join(allid)+"\n")

    infile.close()
    outfile.close()
    print("# ",str(cnt)," orthogroups formatted ("+str(omit)+" omitted).")
def main():
    ###################################
    outfile = None
    infile = None
    prefix = None
    delim = None
    asID = None
    ###################################
    try:
        opts, args = getopt.gnu_getopt(sys.argv[1:], "f:ho:p:l:d:i:", ["file=","help", "output=", "prefix=","loc=","delim=","asID="])
    except getopt.GetoptError as err:
        print (str(err))
        usage()
    for o, a in opts:
        if o in ("-f", "--file"):
            infile=a
        elif o in ("-l", "--loc"):
            locfile=a
        elif o in ("-h", "--help"):
            usage()
        elif o in ("-o", "--output"):
            outfile = a
        elif o in ("-d", "--delim"):
            delim = a
        elif o in ("-p", "--prefix"):
            prefix = a
        elif o in ("-i", "--asID"):
            asID = int(a)

        else:
            assert False, "unhandled option"
    ########################################
    if infile is None:
        usage()
    if outfile is None:
        outfile = infile+".grp"
    if prefix is None:
        usage()
    if asID and not delim or delim and not asID:
        usage()
    ########################################
    print("#  Output: ", outfile)
    read_omcl2grp(infile, locfile, outfile, prefix, delim, asID)

if __name__ == "__main__":
    main()
