#!/bin/env python3

import sys
import getopt
import  logging

logger = logging.getLogger("clean_loc")
logger.setLevel(logging.DEBUG)
# create file handler which logs even debug messages
fh = logging.FileHandler('clean_loc.log')
fh.setLevel(logging.DEBUG)
# create console handler with a higher log level
ch = logging.StreamHandler()
ch.setLevel(logging.ERROR)
# create formatter and add it to the handlers
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
fh.setFormatter(formatter)
ch.setFormatter(formatter)
# add the handlers to the logger
logger.addHandler(fh)
logger.addHandler(ch)

def usage():
    print ("""

    options:
    -l, --loc=FILE.loc
    -o, --output=OUTFILE.loc output file [default: out.loc]
    -h, --help         prints this
    ------------
    """)
    sys.exit(2)

def main():
    out = None
    loc = None
    ###################################
    try:
        opts, args = getopt.gnu_getopt(sys.argv[1:], "l:o:h", ["loc=","output=", "help"])
    except getopt.GetoptError as err:
        print (str(err))
        usage()
    for o, a in opts:
        if o in ("-l", "--loc"):
            loc= a
        elif o in ("-h", "--help"):
            usage()
        elif o in ("-o", "--out"):
            out  = a
        elif o in ("-o", "--output"):
            out = a
        else:
            assert False, "unhandled option"
    ########################################
    if loc is None:
        usage()
    if out is None:
        out = "out.loc"
    logger.debug("loc: {}, out: {}".format(loc,out))
    with open(loc, 'r') as infile:
        res = ""
        for l in infile:
            l = l.strip().split("\t")
            k = [l[0]]
            v = l[1:]
            logger.debug("k{} v{}".format(k, v))
            unique_v = set(v)
            if res:
                res +="\n"
            res+="\t".join(k+list(unique_v))
    with open(out, 'w') as outfile:
        outfile.write(res)

if __name__ == "__main__":
    main()