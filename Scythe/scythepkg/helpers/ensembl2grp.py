import sys, getopt

def usage():
    print ("""
    ######################################
    #  Scythe_ensembleTsv2grp.py (v0.1)  #
    ######################################
    -f, --files=list of ensembl tsv files (eg sA.tsv,sB.tsv,sC.tsv)
    -o, --output=FILE         output file
    -h, --help                prints this
    #----------------------------------#

    """)
    sys.exit(2)
def main():
    ###################################
    outfile = None
    infiles = None
    ###################################
    try:
        opts, args = getopt.gnu_getopt(sys.argv[1:], "f:ho:l:", ["files=","help", "output=","lengths_file="])
    except getopt.GetoptError as err:
        print (str(err))
        usage()
    for o, a in opts:
        if o in ("-f", "--files"):
            infiles = a
        elif o in ("-h", "--help"):
            usage()
        elif o in ("-o", "--output"):
            outfile = a
        else:
            assert False, "unhandled option"

    if not infiles:
        usage()
    if not outfile:
        outfile = "out.grp"
    infiles = infiles.split(",")
    readTsvFiles(infiles, outfile)
########################################

def readTsvFiles(listoftsv, outfile):
    print(outfile)
    if listoftsv is None:
        return(-1)
    ortho=dict()#
    seen = set()
    done = set()
    res = ""
    for t in listoftsv:
        infile = open(t, "r")
        for l in infile:
            l = l.strip()
            l = l.split("\t")
            seen.add(l[0])
            seen.add(l[1])

            if l[0] not in ortho:
                ortho[l[0]] = [l[1]]
            else:
                ortho[l[0]].append(l[1])
            if l[1] not in ortho:
                ortho[l[1]] = [l[0]]
            else:
                ortho[l[1]].append(l[0])

    cntr = 0
    for s in seen:
        notyetdone = [o for o in ortho[s] if o not in done]
        if s not in done and notyetdone:
            res+=str(cntr)+"\t"+s+"\t"+"\t".join(notyetdone)
            res+="\n"
            cntr+=1
            done.add(s)
            for d in ortho[s]:
                done.add(s)

    out=open(outfile,"w")
    out.write(res)
if __name__ == "__main__":
    main()
