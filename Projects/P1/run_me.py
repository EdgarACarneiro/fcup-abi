from bioseq import BioSeq


if __name__ == '__main__':
    print("Hello and welcome to the demonstrarion of the module 'BioSeq'!\n")
    print("The program will stop at checkpoints and you, as the User, will have to press Enter to continue. Ok?\n")
    input("Press Enter to continue...")

    print("Ok, nice!\n Lets create a dna sequence with the sequence 'AGCTAGCTAGCTACATCAAACGATCGTCGCTAGCTAGCTAC', ok?\n")
    input("Press Enter to continue...")

    dna = BioSeq.create_bio_seq("AGCTAGCTAGCTACATCAAACGATCGTCGCTAGCTAGCTAC")
    dna.pretty_print()

    print("Now, lets load a genetic_code dictionary!\n")
    dna.read_genetic_code('tests/files/genetic_code.txt')
    input("Press Enter to continue...")
    dna.pretty_print()

    input("Press Enter to continue...")
    print("Ok, look, we can also transcript our dna into a RNA!\n")
    print(str(dna.transcription()))

    input("Press Enter to continue...")
    print("We do have a genetic code dictionary and a dna, so we will translate it now!\n")
    input("Press Enter to continue...")
    print(str(dna.translate()))

    input("Press Enter to continue...")
    print("We can also see the codon usage of our sequence for the 'R' aminoacid!\n")
    input("Press Enter to continue...")
    print(dna.codon_usage("R"))

    input("Press Enter to continue...")
    print("Lets see all the obtainable reading frames!\n")
    input("Press Enter to continue...")
    for rf in dna.reading_frames():
        print(str(rf))

    input("Press Enter to continue...")
    print("We can now try other dna sequence! What about this one? 'ATGAAATTATGAATGAGCCTCAGCTGAAGCATCGCGCATCAGACTACGCTCAGACTCAGACTCAGCATTATAGTGAATGTTAATAAATAAAATAA' \n")
    input("Press Enter to continue...")
    print("Lets see all its opeaning reading frames!\n")
    dna = BioSeq.create_bio_seq("ATGAAATTATGAATGAGCCTCAGCTGAAGCATCGCGCATCAGACTACGCTCAGACTCAGACTCAGCATTATAGTGAATGTTAATAAATAAAATAA")
    dna.read_genetic_code('tests/files/genetic_code.txt')

    for orfs in dna.all_orfs():
        print(orfs.get_seq())

    input("Press Enter to continue...")

    print("\n-------------------------\n")

    print("Hey, lets try an RNA now! Can it be this one? 'UCGA'?\n")
    rna = BioSeq.create_bio_seq("UCGA", "rna")
    
    input("Press Enter to continue...")
    print("Lets see the frequency of its symbols!\n")
    for (k, v) in rna.freq_symbols().items():
        print(k + ': ' + str(v))
    input("Press Enter to continue...")

    print("Lets see its gc percent!\n")
    print(rna.gc_percent())
    input("Press Enter to continue...")

    print("And what about its reverse complement?\n")
    print("Reverse Complement: " + rna.reverse_complement())
    input("Press Enter to continue...")

    print("We can also save and load sequences! Let's do it with the RNA to the 'test/files/test_save_load.csv' file!\n")
    input("Press Enter to continue...")
    rna.save('tests/files/test_save_load.csv')
    print("The rna is now saved, open the csv!\n")
    input("And now we print what we load from that file!\n")
    BioSeq.load('tests/files/test_save_load.csv').pretty_print()

    input("Press Enter to continue...")

    print("\n-------------------------\n")

    print("Finally lets see the contents of the FASTA file in the 'tests/files/sequence.fasta'!\n")
    input("Press Enter to continue...")
    fasta = BioSeq.read_fasta_file('tests/files/sequence.fasta')
    for (k, v) in fasta.items():
        print('> ' + k + '\n' + v + '\n\n')