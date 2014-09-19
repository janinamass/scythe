#!/usr/bin/env python

import re,sys, getopt

def usage():
    print ("""
    usage:  scythe_mergeSubsets.py -g groups.grp -o new.grp

    options:
    -g, --grp=FILE.grp
    -o, --output=OUTFILE.grp    output file [default: FILE.allspec.grp]
    [-r, --rename               discard old orthogroup ids and start numbering from 0]
    -h, --help                  prints this
    [-n, --numspec=N    min number of species ]
    ------------
    .grp format: GroupID\tgeneIDiSp1\tgeneIDjSp2\t...geneIDkSpn
    """)
    sys.exit(2)


class Orthogroup(object):
    def __init__(self,list):
        self._list = list[1:]
        self._isSubset = None
        self._size = len(self._list)
        self._id = list[0]
        self._ignore = False
    def expand(self, other):
        tmpset = set()
        tmplist = []
        dct = {}
        purge= False
        for i in self._list:
            tmplist.append(i)
            if i[0:6] in dct:
                if i != dct[i[0:6]]:
                    purge=True
            else:
                dct[i[0:6]]=i
        for j in other._list:
            if j[0:6] in dct:
                if j != dct[j[0:6]]:
                    purge=True
            else:
                dct[j[0:6]]=j
            tmplist.append(j)
        tmpset = set(tmplist)
        newlist=list(tmpset)
        newid= str(self._id)+"_"+str(other._id)

        initlist=[newid]
        for i in newlist:
            initlist.append(i)
        new = Orthogroup(initlist)
        if purge:
            new._isSubset=True
            new._ignore=True
        return(new)
    def toString(self):
        s = "\t".join(self._list)
        return(s)

def mergeSubsets(grp, outfile, rename=False):
    write = True
    if outfile==None:
        write = False
    knownGenes= set()
    kglist = []
    knownOrthogroups = []
    godict = {}
    proc = set()
    try:
        g = open(grp,'r')
    except IOError as e:
        print(e)
    for l in g:
            l = l.rstrip()
            tmp = l.split("\t")
            newOrtho= Orthogroup(tmp)
            for t in tmp[1:]:
                if t in knownGenes:
                    tmplist = []
                    for o in godict[t]:
                        new = newOrtho.expand(o)
                        new._isSubset=False
                        o._isSubset = True
                        newOrtho._isSubset = True
                        tmplist.append(new)
                    for l in tmplist:
                        godict[t].append(l)
                else:
                    knownGenes.add(t)
                    godict[t] = [newOrtho]
    g.close()
    try:
        g = open(grp,'r')
    except IOError as e:
        print(e)
    processed =[]
    for l in g:
        l = l.rstrip()
        tmp = l.split("\t")
        for t in tmp[1:]:

            nmm = [l for l in godict[t] if  l._isSubset !=True and l._ignore==False]
            out=None
            if nmm !=[]:
                out = nmm[-1]
                if out._id not in proc:
                    processed.append(out._id+"\t"+out.toString())

                proc.add(out._id)
    g.close()
    if write:
        try:
            o = open(outfile,'w')
        except IOError as e:
            print(e)

        seen = set()

        cnt = 0
        ok = 0
        redundant = 0
        isOK = True
        for p in processed:
            isOK = True
            tmp = p.split("\t")
            for t in tmp[1:]:
                if t in seen:
                    isOK = False
                seen.add(t)
            if isOK:
                ok+=1
                if not rename:
                    o.write(str(p))
                if rename:
                    o.write(str(cnt)+"\t"+"\t".join(p.split("\t")[1:]))
                o.write("\n")
                cnt+=1
            else:
                redundant+=1
        o.close()
    return(processed)

def filterGroups(grp, outfile, numspec=None,rename=False):
    write = True
    if outfile==None:
        write = False
    try:
        g = open(grp,'r')
    except IOError as e:
        print(e)
        usage()

    if numspec:
        max = numspec
    if not numspec:
        max = 0
        for l in g:
            l = l.strip()
            tmp = len(l.split("\t")[1:])
            if tmp > max:
                max = tmp
        g.close()
    numspec=max

    if write:
        try:
            g = open(grp,'r')
        except IOError as e:
            print(e)
            usage()
        try:
            out = open(outfile,'w')
        except IOError as e:
            print(e)
            usage()

        cnt=0
        for l in g:
            ln = l.rstrip()
            tmp = len(ln.split("\t")[1:])
            if tmp == numspec:
                if not rename:
                    out.write(l)
                else:
                    out.write(str(cnt)+"\t"+"\t".join(l.split("\t")[1:]))
                    cnt+=1
        out.close()
    return(numspec)
def main():
    out = None
    outNoSubsets = None
    rename=False
    numspec = None
    grp = None
    ###################################
    try:
        opts, args = getopt.gnu_getopt(sys.argv[1:], "g:n:o:hr", ["grp=","numspec=","output=", "help","rename"])
    except getopt.GetoptError as err:
        print (str(err))
        usage()
    for o, a in opts:
        if o in ("-g", "--grp"):
            grp= a
        elif o in ("-h", "--help"):
            usage()
        elif o in ("-n", "--numspec"):
            numspec = int(a)
        elif o in ("-o", "--output"):
            out = a
            tmp=a
        elif o in ("-r", "--rename"):
            rename = True
        else:
            assert False, "unhandled option"
    ########################################
    if grp is None:
        usage()
    if out is None:
        tmp="tmp."
        if grp.endswith(".grp"):
            tmp = grp[:-3]
    out = tmp+"shared_by_all.grp"
    outNoSubsets = tmp+"full.grp"
    ########################################
    print("# Output: ", out)
    print("# Output: ", outNoSubsets)

    if not numspec:
        numspec= filterGroups(grp,None, None, rename)
    print(numspec)
    mergeSubsets(grp,outNoSubsets, rename)
    filterGroups(outNoSubsets, out, numspec, True)



if __name__=="__main__":
    main()

