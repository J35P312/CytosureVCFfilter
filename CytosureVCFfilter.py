import argparse

parser = argparse.ArgumentParser("""Filter SV vcf files for usage with Marcel's cytosure2vcf tools""",add_help=False)
parser.add_argument('--vcf'        ,required=True,type=str, help="minimum variant size, smaller variants will be removed")
parser.add_argument('--size'        , default=5000,type=int, help="minimum variant size, smaller variants will be removed")
parser.add_argument('--tiddit'        , default=8,type=int, help="TIDDIT support, remove variants which are called only by tiddit if their support is too low(default = 9)")
parser.add_argument('--nator'        , default=0.3,type=float, help="minimum CNVnator normalised RD deviation")
parser.add_argument('--frequency'        , default=0.01,type=float, help="Maximum frequency(default = 0.01)")
parser.add_argument('--frequency_tag'        , default="FRQ",type=str, help="the freqyency tag of the info field")
args, unknown = parser.parse_known_args()


for line in open(args.vcf):
    ok=False
    if line[0] == "#":
        print line.strip()
        continue
    content=line.strip().split("\t")

    #test sv size, skip if too small
    if "END=" in content[7]:
        end=int(content[7].split("END=")[-1].split(";")[0].split("\t")[0])
        svlen=abs(int(content[1])-end)
        if svlen < args.size:
            continue
    if not "END" in content[7] and "WINA=" in content[7] and "WINB=" in content[7] and ":" in content[4]:
        chromosomeA=content[0]
        posA=int(content[1])
        chromosomeB=content[4].split("N[")[-1].split(":")[0]
        posB=int(content[4].split(":")[-1].split("[")[0])
        if chromosomeA == chromosomeB and abs(posA-posB) < args.size:
            continue

    #test tiddit support, flag as low quallity if too low 
    if "WINA=" in content[7] and "WINB=" in content[7]:
        support=0
        if ";LTE" in content[7]:
            support+=int(content[7].split(";LTE=")[-1].split(";")[0].split("\t")[0])            
        if" ;SR" in content[7]:
            support+=int(content[7].split(";SR=")[-1].split(";")[0].split("\t")[0])   
        if support >= args.tiddit:
            ok=True

    #test cnvnator support, discard if no tiddit support and too low cnvnator support; flag as ok if CNVnator support is ok, but tiddit support is too low
    if "natorRD=" in content[7]:
        normRD=abs(float(content[7].split(";natorRD=")[-1].split(";")[0])-1)
        if normRD < args.nator and not ok:
            continue
        else:
            ok=True

    if args.frequency_tag in content[7]:
        frequency=float(content[7].split(";{}=".format(args.frequency_tag) )[-1].split(";")[0].split("\t")[0])
        if frequency > args.frequency:
            continue    

    #variants which are rare enough and called by both tiddit and cnvnator are ok regardless of support    
    if "natorRD=" in content[7] and "WINA=" in content[7] and "WINB=" in content[7]:
        ok=True

    #tiddit and cnvnator filters will be skipped if the variant is not called using those callers
    if not ok and not "natorRD=" in content[7] and not "WINB=" in content[7]:
        ok=True
    #remove second bnd of a variant
    if "WINA=" in content[7] and "WINB=" in content[7] and not "_1\tN\t" in line:
        continue

    if ok:
        print line.strip()
