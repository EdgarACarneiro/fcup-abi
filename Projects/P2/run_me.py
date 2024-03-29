from bioseq import BioSeq
from seq_align import read_submat_file,\
                        subst_matrix,\
                        pretty_matrix,\
                        global_align_multiple_solutions,\
                        recover_global_align_multiple_solutions,\
                        local_align_multiple_solutions,\
                        recover_local_align_multiple_solutions,\
                        compare_pairwise_global_align,\
                        compare_pairwise_local_align,\
                        compare_pairwise_num_global_align,\
                        compare_pairwise_num_local_align

def wait_input():
    input("Press Enter to continue...\n")

if __name__ == '__main__':

    print("Hello and welcome to the demonstrarion of the module 'SeqAlign'!\n")
    print("The program will stop at checkpoints and you, as the User, will have to press Enter to continue. Ok?\n")
    wait_input()

    print("Ok, nice!\nLets start with two sequences: 'GATTACA' & 'GCATGCT', ok?\n")
    wait_input()

    seq1 = 'GATTACA'
    seq2 = 'GCATGCT'

    print("Now we can either create a substitution matrix or read one from a file!\n")
    print("Lets create one first with the alphabet 'ACTG' and 1 for matches and -1 for mismatches.\n")
    wait_input()

    sm_dna = subst_matrix("ACGT", 1, -1)
    print(sm_dna)
    print("We will call this substitution matrix 'sm_dna', ok?\n")
    wait_input()

    print("Now, lets load one from a file.\nCan we use the one in the file 'tests/files/blosum62.mat'? Is it okay for you?\n")
    wait_input()

    sm_blosum = read_submat_file('tests/files/blosum62.mat')
    print(sm_blosum)
    print("We will call this substitution matrix 'sm_blosum', ok?\n")
    wait_input()

    print("Since we have a substitution matrix and the sequences we can get some global and local alignments!\n\n")
    print("::: GLOBAL ALIGNMENT with multiple solutions :::\n")
    print("\nLets make the global alignment with our sequences, the 'sm_dna' and a gap of -1, shall we?\n")
    wait_input()

    ga_score, ga_trace = global_align_multiple_solutions(seq1, seq2, sm_dna, -1)
    print("Score matrix obtained:\n")
    pretty_matrix(ga_score, " " + seq1, " " + seq2)
    wait_input()

    print("\nTrace matrix obtained:\n")
    pretty_matrix(ga_trace, " " + seq1, " " + seq2)
    wait_input()

    print("Using the computed matrixes we can recover the multiple optimal global alignments.\nLets do it?\n")
    wait_input()

    rga = recover_global_align_multiple_solutions(ga_trace, seq1, seq2)
    print("The obtained global alignments were:\n")
    for align in rga:
        print('Seq1: ' + str(align[0]) +  '\nSeq2: ' + str(align[1]) + '\n')
    wait_input()

    print("::: LOCAL ALIGNMENT with multiple solutions :::\n")
    print("\nLets make the local alignment with our sequences, the 'sm_dna' and a gap of -1, shall we?\n")
    wait_input()

    ga_score, ga_trace, max_score = local_align_multiple_solutions(seq1, seq2, sm_dna, -1)
    print("Score matrix obtained:\n")
    pretty_matrix(ga_score, " " + seq1, " " + seq2)
    wait_input()

    print("\nTrace matrix obtained:\n")
    pretty_matrix(ga_trace, " " + seq1, " " + seq2)
    wait_input()

    print("\nMaximum value in the score matrix: " + str(max_score) + '\n')
    wait_input()

    print("Using the computed matrixes we can recover the multiple optimal local alignments.\nLets do it?\n")
    wait_input()

    rga = recover_local_align_multiple_solutions(ga_score, ga_trace, seq1, seq2)
    print("The obtained local alignments were:\n")
    for align in rga:
        print('Seq1: ' + str(align[0]) +  '\nSeq2: ' + str(align[1]) + '\n')
    wait_input()

    print("\n-------------------------\n")

    print("Now shall we make it a little bit more complex?\n")
    wait_input()

    print("So, lets load and use the proteins present in the 'tests/files/protein_sequences.fas' file!\n\
        We need to use the read fasta functionality from the bioseq library to load them.\n")
    wait_input()

    print("The protein sequences are:\n")
    seqs = BioSeq.read_fasta_file('tests/files/protein_sequences.fas')
    for key, value in seqs.items():
        print("> " + str(key) + ": " + str(value) + "\n")
    wait_input()

    print("::: GLOBAL ALIGNMENT with multiple solutions :::\n")
    print("\nLets make the global alignment with two of our protein sequences: sp|B0C882: & sp|A1TQI0, the 'sm_dna' and a gap of -3, shall we?\n")
    p_seq1 = seqs["sp|B0C882"]
    p_seq2 = seqs["sp|A1TQI0"]
    wait_input()

    ga_score, ga_trace = global_align_multiple_solutions(p_seq1, p_seq2, sm_blosum, -3)
    print("Score matrix obtained:\n")
    pretty_matrix(ga_score, " " + p_seq1, " " + p_seq2)
    wait_input()

    print("\nTrace matrix obtained:\n")
    pretty_matrix(ga_trace, " " + p_seq1, " " + p_seq2)
    wait_input()

    print("Using the computed matrixes we can recover the multiple optimal global alignments.\nLets do it?\n")
    wait_input()

    rga = recover_global_align_multiple_solutions(ga_trace, p_seq1, p_seq2)
    print("The obtained global alignments were:\n")
    for align in rga:
        print('Seq1: ' + str(align[0]) +  '\nSeq2: ' + str(align[1]) + '\n')
    wait_input()

    print("::: LOCAL ALIGNMENT with multiple solutions :::\n")
    print("\nLets make the local alignment with two other protein sequences: sp|A0KL54 & sp|B7JC18, the 'sm_dna' and a gap of -1, shall we?\n")
    p_seq1 = seqs["sp|A0KL54"]
    p_seq2 = seqs["sp|B7JC18"]
    wait_input()

    ga_score, ga_trace, max_score = local_align_multiple_solutions(p_seq1, p_seq2, sm_blosum, -3)
    print("Score matrix obtained:\n")
    pretty_matrix(ga_score, " " + p_seq1, " " + p_seq2)
    wait_input()

    print("\nTrace matrix obtained:\n")
    pretty_matrix(ga_trace, " " + p_seq1, " " + p_seq2)
    wait_input()

    print("\nMaximum value in the score matrix: " + str(max_score) + '\n')
    wait_input()

    print("Using the computed matrixes we can recover the multiple optimal local alignments.\nLets do it?\n")
    wait_input()

    rga = recover_local_align_multiple_solutions(ga_score, ga_trace, p_seq1, p_seq2)
    print("The obtained local alignments were:\n")
    for align in rga:
        print('Seq1: ' + str(align[0]) +  '\nSeq2: ' + str(align[1]) + '\n')
    wait_input()

    print("\n-------------------------\n")

    print("We can also obtain the pairwise global and local alignment of the the sequences belonging to a list:\n")
    wait_input()

    print("Lets try it out with the protein sequences we had!\nFirst, the comparation of scores of the pairwise global alignments:\n")
    wait_input()

    seqs = list(seqs.values())
    print("This might take a little bit...")
    compare_pairwise_global_align(seqs, sm_blosum, -3)
    wait_input()

    print("And now, the comparation of scores of the pairwise local alignments:\n")
    wait_input()

    print("This might take a little bit...")
    compare_pairwise_local_align(seqs, sm_blosum, -3)
    wait_input()

    print("\n::: EXTRA :::\n")

    print("Finally, we can also obtain the pairwise global and local alignment of the the sequences belonging to a list:\n")
    wait_input()

    print("Lets try it out with the protein sequences we had!\nFirst, the comparation of the number of pairwise global alignments:\n")
    wait_input()

    print("This might take a little bit...")
    compare_pairwise_num_global_align(seqs, sm_blosum, -3)
    wait_input()

    print("And now, the comparation of the number of pairwise local alignments:\n")
    wait_input()

    print("This might take a little bit...")
    compare_pairwise_num_local_align(seqs, sm_blosum, -3)
    wait_input()

    print("And that is it!! This are this module functionalities, hope you enjoyed it!")
    wait_input()