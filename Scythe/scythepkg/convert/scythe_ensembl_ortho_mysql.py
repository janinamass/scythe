#! /usr/bin/env python
import mysql.connector
import sys
import getopt

def q(s):
    return ('"'+s+'"')

def getGenome_Db_Ids(specnames, release):
    cmd = 'use ensembl_compara_'+str(release)+";"
    cnx = mysql.connector.connect(user='anonymous', host='ensembldb.ensembl.org' )
    curA = cnx.cursor(buffered=True)
    curA.execute(cmd)
    specnames = [q(s) for s in specnames]
    res = dict()
    for s in specnames:
        cmd = 'SELECT genome_db_id FROM genome_db WHERE name='+s+";"
        curA.execute(cmd)
        res[s] = curA.fetchall()[0][0]
    return(res)

def getSpecies_Set_Ids(genomedb_ids, release):
    cmd = 'use ensembl_compara_'+str(release)+";"
    cnx = mysql.connector.connect(user='anonymous', host='ensembldb.ensembl.org' )
    curA = cnx.cursor(buffered=True)
    curA.execute(cmd)
    genomedb_ids=[str(s) for s in genomedb_ids]
    s = "genome_db_id  IN ("+ ",".join(genomedb_ids)+") "
    cmd = "SELECT species_set_id,genome_db_id, count(species_set_id )"
    cmd += "from species_set  WHERE ("
    cmd += s
    cmd += ") group by species_set_id having count(species_set_id) > 1  ;"
    curA.execute(cmd)
    res = curA.fetchall()
    res = [r[0] for r in res]
    return (res)

def getMethodLinkSpecies_Set_Ids(speciesSet_ids, release):
    cmd = 'use ensembl_compara_'+str(release)+";"
    cnx = mysql.connector.connect(user='anonymous', host='ensembldb.ensembl.org' )
    curA = cnx.cursor(buffered=True)
    curA.execute(cmd)
    speciesSet_ids=[str(s) for s in speciesSet_ids]
    s = "species_set_id  IN ( "+ ",".join(speciesSet_ids)+") "
    cmd = "SELECT method_link_species_set_id, species_set_id"
    cmd += " FROM method_link_species_set  WHERE "
    cmd += s
    cmd += " ;"
    curA.execute(cmd)
    res = curA.fetchall()
    curA.close()
    res = [r[0] for r in res]
    return (res)

def getHomologyId(method_link_species_set_ids, release=71):
    cmd = 'use ensembl_compara_'+str(release)+";"
    cnx = mysql.connector.connect(user='anonymous', host='ensembldb.ensembl.org' )
    curA = cnx.cursor(buffered=True)
    curA.execute(cmd)
    method_link_species_set_ids=[str(s) for s in method_link_species_set_ids]
    #s = "method_link_species_set_id ="+" OR method_link_species_set_id=".join(method_link_species_set_ids)
    s = "method_link_species_set_id IN ("+",".join( method_link_species_set_ids)+") "
    cmd = "SELECT homology_id, method_link_species_set_id"
    cmd += " FROM homology WHERE ("
    cmd += s
    cmd += ");"
    curA.execute(cmd)
    res = curA.fetchall()
    curA.close()
    return (res)

def getHomologyMemberId(homology_ids, release):
    cmd = 'use ensembl_compara_'+str(release)+";"
    cnx = mysql.connector.connect(user='anonymous', host='ensembldb.ensembl.org' )
    curA = cnx.cursor(buffered=True)
    curA.execute(cmd)
    homology_ids=[str(s) for s in homology_ids]
    s = "homology_id ="+" OR homology_id=".join(homology_ids)
    cmd = "SELECT homology_member.member_id, member.member_id,member.stable_id, homology_id "
    cmd += " FROM homology_member LEFT JOIN  member ON "
    cmd+= "(homology_member.member_id = member.member_id) WHERE ("
    cmd += s
    cmd += ");"
    curA.execute(cmd)
    res = curA.fetchall()
    curA.close()
    return (res)

def fetch1to1orthologs(method_link_species_set_ids,release):
    """Return (homology_id, genome_db_id, stable_id)"""
    method_link_species_set_ids=[str(s) for s in method_link_species_set_ids]

    cmd = 'use ensembl_compara_'+str(release)+";"
    cnx = mysql.connector.connect(user='anonymous', host='ensembldb.ensembl.org' )
    curA = cnx.cursor(buffered=True)
    curA.execute(cmd)
    s = "method_link_species_set.method_link_species_set_id IN ("+ ",".join(method_link_species_set_ids)+") "
    cmd = " SELECT DISTINCT homology.homology_id,  member.genome_db_id,  member.stable_id from"
    cmd +=" homology, homology_member, member, species_set, method_link_species_set "
    cmd +=" WHERE ( "
    cmd +=s
    cmd +='AND homology.description="ortholog_one2one"'
    cmd +="AND "+s
    cmd +="AND species_set.species_set_id = method_link_species_set.species_set_id "
    cmd +="AND homology.method_link_species_set_id = method_link_species_set.method_link_species_set_id "
    cmd +="AND homology_member.member_id=member.member_id "
    cmd +="AND homology_member.homology_id=homology.homology_id "
    cmd +=")"

    curA.execute(cmd)
    res = curA.fetchall()
    curA.close()
    return (res)



def makeTable(genomeDB_ids, orthoSpec2stableIds, orthoIds, outfile):
    global SPECID
    global OUTFILES
    res = dict()
    specset = dict()
    for o in orthoIds:
        for i in genomeDB_ids:
            if (o,i) in orthoSpec2stableIds:
                if not o in res:
                    res[o]=[i]
                else:
                    res[o].append(i)
    for r in res:
        #print(res[r][0],res[r][1])
        if not ((res[r][0],res[r][1])) in specset:
                 specset[(res[r][0],res[r][1])]=[(orthoSpec2stableIds[(r,res[r][0])],orthoSpec2stableIds[(r,res[r][1] )]) ]

        else:
             specset[(res[r][0],res[r][1])].append((orthoSpec2stableIds[(r,res[r][0])],orthoSpec2stableIds[(r,res[r][1])]))

    #print(specset)
    for s in specset:
        st = SPECID[s[0]]+"__"+SPECID[s[1]]
        print(st)
        outfile = st
        out = open(outfile,"w")
        print("writing to {}".format(outfile))
        OUTFILES.append(outfile)
        tuplist = specset[s]
        for tup in tuplist:
            tup = str(tup[0])+"\t"+str(tup[1])
            out.write(tup)
            out.write("\n")
        out.close()


def fetchOrthoFromMySQL(specieslist = [], release=None, outfile=None):
    if not release :
        sys.exit("No release specified.")
    if len(specieslist) < 2:
        sys.exit("Need more than one species.")
    print(specieslist)
    print("...")
    name2genomeDB_ids = getGenome_Db_Ids(specieslist, release)
    #print(name2genomeDB_ids)
    for k,v in name2genomeDB_ids.items():
        SPECID[k] = v
        SPECID[v] = k[1:-1]
    orthoIds = []
    orthoSpec2stableIds = dict()
    members = dict()

    genomeDB_ids = name2genomeDB_ids.values()
    #print(genomeDB_ids)
    speciesSet_ids = getSpecies_Set_Ids(genomeDB_ids,release)
    #print(speciesSet_ids,"specsetids")
    methodLinkSpeciesSet_ids = getMethodLinkSpecies_Set_Ids(speciesSet_ids, release)
    #print(methodLinkSpeciesSet_ids)
    res = fetch1to1orthologs(methodLinkSpeciesSet_ids, release=release)
    print("...")
    for r in res:
        if not (r[0],r[1]) in orthoSpec2stableIds:
            #r[0] is shared between orthologous genes from different species
            #r[1] is the genome_id
            #r[2] is the stable_id of the gene
            if not r[0] in orthoIds:
                orthoIds.append(r[0])
            orthoSpec2stableIds[(r[0],r[1])]=r[2]#[0]

    print("writing table\n")
    makeTable(genomeDB_ids, orthoSpec2stableIds, orthoIds, outfile)
    return(OUTFILES)

def usage():
    """Usage"""
    print("""

    usage: scythe_ensembl_ortho_mysql.py -s species_1,species_2 -r INT

    options:
    -s, --species=STR   comma-separated list of species (eg 'homo_sapiens,gorilla_gorilla')
    -r, --release=NUM   ensembl version (eg '75')
    -h, --help          prints this
    """)
    sys.exit(2)

SPECID = {}
OUTFILES = []

def main():
    speciesList = None
    release = None

    try:
        opts, args = getopt.gnu_getopt(sys.argv[1:], "s:r:h",["species=","release=","help"])
    except getopt.GetoptError as err:
        sys.stderr.write(str(err))
        usage()
    for o, a in opts:
        if o in ("-s", "--species"):
            speciesList = a.split(",")
        elif o in ("-r", "--release"):
            release = int(a)
        elif o in ("-h", "--help"):
            usage()
        else:
            assert False, "unhandled option"
    if not speciesList:
        usage()
    if not release:
        usage()

    fetchOrthoFromMySQL(specieslist = speciesList, release = release )

if __name__=="__main__":
    main()

