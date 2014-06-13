#!/usr/bin/env python
import re
import sys
import getopt

def usage():
    print ("""
    ################################
    #  Scythe_proteinOrtho2grp.py  #
    ################################

    -f, --file=INFILE        (output of ProteinOrtho, pairs, unambiguous)
    -l, --loc=LOCFILE        concatenated .loc file or comma separated list of .loc files
    -o, --output=OUTFILE     output file [default: INFILE.grp]
    -h, --help               prints this

    ------------
    .grp format: GroupID\tgeneIDiSp1\tgeneIDjSp2\t...geneIDkSpn
    """)
    sys.exit(2)

def read_proteinortho2grp(infile, locfile, outfile):
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
    grp= 0
    outfile = open(outfile, 'w')
    infile = open(infile, 'r')
    for ln in infile:
        if ln.startswith("#"):
            continue
        ln = ln.rstrip()
        ids = ln.split("\t")[3:]
        allid = []
        for i in ids:
            if i == "*":
                continue
            allid.append(mod2loc[i])
        outfile.write(str(grp)+"\t"+"\t".join(allid)+"\n")
        grp+=1
    infile.close()
    print("# ",str(grp)," orthogroups formatted.")

def main():
    ###################################
    outfile = None
    infile = None
    locfile = None
    ###################################
    try:
        opts, args = getopt.gnu_getopt(sys.argv[1:], "f:ho:l:", ["file=","help", "output=","loc="])
    except getopt.GetoptError as err:
        print (str(err))
        usage()
    for o, a in opts:
        if o in ("-f", "--file"):
            infile=a
        elif o in ("-h", "--help"):
            usage()
        elif o in ("-o", "--output"):
            outfile = a
        elif o in ("-l","--loc"):
            locfile = a
            print(a)
        else:
            assert False, "unhandled option"
    ########################################
    if infile is None:
        usage()
    if outfile is None:
        outfile = infile+".grp"
    if locfile is None:
        usage()
    ########################################
    print("#  Output: ", outfile)
    read_proteinortho2grp(infile, locfile, outfile)

if __name__ == "__main__":
    main()
