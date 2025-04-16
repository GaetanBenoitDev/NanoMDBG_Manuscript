
import os, sys, argparse



def main(argv):
    
    parser = argparse.ArgumentParser()

    parser.add_argument("gfa", help="input gfa file")
    parser.add_argument("fasta", help="output fasta file")
    
    args = parser.parse_args()

    gfa_to_fasta(args.gfa, args.fasta)

def gfa_to_fasta(gfa_filename, fasta_filename):

    nbUnitigs = 0
    #unitig_sizes = {}

    fasta_file = open(fasta_filename, "w")

    for line in open(gfa_filename):
        line = line.strip()
        fields = line.split("\t")

        if fields[0] == "S":
            unitig_name = fields[1]
            unitig_sequence = fields[2]
            fasta_file.write(">" + unitig_name + "\n")
            fasta_file.write(unitig_sequence + "\n")
            nbUnitigs += 1
            #unitig_sizes[unitig_name] = len(unitig_sequence)

    print("nb unitigs: " + str(nbUnitigs))

    fasta_file.close()
    #return unitig_sizes

if __name__ == "__main__":
    main(sys.argv[1:])  
