from src import Pipeline, UPGMA, MyGraph, MyBlast, MultipleAlignment, SubstMatrix, BioSeq, Seq, Dna, Rna, Protein


def wait_input():
    input("Press Enter to continue...\n")

def take_while():
    print("This might take a little bit...")

if __name__ == '__main__':

    print("Hello and welcome to the demonstration of the module 'Pipeline'!\n")
    print("The program will stop at checkpoints and you, as the User, will have to press Enter to continue. Ok?\n")
    wait_input()

    print("Ok, nice!\nFirst, lets perform a step-by-step demonstration of the Pipeline and then execute a full automatic Pipeline run, ok?\n")
    wait_input()

    print("Our Pipeline module constructor accepts to strings, representing the filenames of the fasta files.\n" +
          "We will use the fasta file 'src/tests/files/source.fasta' for the query sequence and the fasta file 'src/tests/files/seqdump.txt' for the database.\n")
    wait_input()

    p = Pipeline('src/tests/files/source.fasta',
                 'src/tests/files/seqdump.txt', 10)

    print("As you can see, our query sequence is:\n\t%s - %s" %
          (p.query_id, p.query_seq.get_seq()))
    print("Noticed that it output a Protein? That is because our Pipeline is capable of infering the Sequence type!\n")
    wait_input()

    print("Now lets look at our database!\nDATABASE:\n")
    for el in p.database:
        print("\t%s - %s" % (el[0], el[1].get_seq()))
    wait_input()

    print("Also, we must nor forget to set the the alignment settings that we will use, and by that I mean the Substitution Matrix and the gap penalty value.\n")
    print("One can either dinamically create a substitution matrix or load one from a file.\nIn this demo lets use the famous 'blosum62' substitution matrix and a gap penalty of -8. Its that ok for you?\n")
    wait_input()

    p.change_alignment_settings(SubstMatrix.read_submat_file(
        'src/test/files/blosum62.mat'), -8)

    print("We can now see the Alignment configuration that was used:\nSubstitution Matrix:\n")
    print('\t%s' % str(p.align_config[0]))
    print("Gap Penalty:\n\t%s" % str(p.align_config[1]))
    print("The Alignment Configuration can be changed using the Pipeline::change_alignment_settings() function.\n")
    wait_input()

    print("Ok, now that we are setup, shall we start executing our Pipeline?")
    wait_input()

    print("First, we will create a database copy were only the sequences of different specie from the query sequence are featured.\n")
    wait_input()
    print("\n\t:::Step 1 - Create copy database without similar specie sequences:::\n")
    db_copy = [(_id, seq)
               for _id, seq in p.database
               if Pipeline.get_specie_from_id(_id) !=
               p.get_specie()]

    print("This leaves us with the following database copy:\n\tSPECIE - SEQUENCE")
    for al in db_copy:
        print(p.get_specie_from_seq(al), al.get_seq())
    wait_input()

    print("Second, lets run the BLAST method so we can get the 10 best alignments between the query sequence and the remaining sequences of different species!\n")
    wait_input()
    print("\n\t:::Step 2 - Running Blast and getting top alignments:::\n")
    blast = MyBlast()
    for _, seq in db_copy:
        blast.add_sequence_database(str(seq))
    take_while()
    top_alignments = blast.best_alignments(str(p.query_seq), p.TOP)

    print("Top %d Alignments obtained from Blast" % p.TOP)
    for al in top_alignments:
        print(al)
    wait_input()

    print("Now, lets map the obtained top alignments into the respective sequences and join them with the query sequence so we can after perform some operations over that set of sequences.\n")
    wait_input()

    db_seqs = [val for _, val in db_copy]
    best_seqs = [p.query_seq] + list(map(
        lambda align: db_seqs[align[4]], top_alignments))
    for seq in best_seqs:
        print(p.get_specie_from_seq(seq), seq.get_seq())
    wait_input()

    print("For the next step, the run the Multiple Sequence Alignment (MSA) Algorithm over our set of sequences.\n")
    print("\n\t:::Step 3 - Running MSA with the respective top alignments:::\n")
    # Multiple Sequence Alignment
    msa = MultipleAlignment(best_seqs, p.align_config).align_consensus()
    take_while()

    # Printing the Multiple Sequence Alignemnt
    print("Multiple Sequence Alignment Result:\n")
    msa.pretty_print()

    print("\nAs the result of the MSA we obtain the evolution of the consensus sequence! Interesting isn't it?")
    wait_input()

    print("For the fourth step, lets obtain the Phylogenetic tree resultant of our set of sequences, using the UPGMA method!")
    wait_input()

    print("\n\t:::Step 4 - Obtaining the Ultrametric Tree from the top alignments:::\n")
    upgma = UPGMA(best_seqs, p.align_config)
    take_while()

    # Printing UPGMA distance matrix
    print("Here is the Distance Matrix for our set of sequences computed by the UPGMA method:")
    print("Distances Matrix obtained by the UPGMA method:")
    upgma.dists_mat.print_mat()

    # Producing the Ultrametric Tree
    tree = upgma.run()
    take_while()

    # Printing the tree with a mapping for the species
    print("Phylogenetics Tree created:\n")
    tree.print_tree({
        idx: p.get_specie_from_seq(seq)
        for idx, seq in enumerate(best_seqs)
    })
    wait_input()

    print("That was awesome!\nNow, for the final step, lets assemble that same information in a graph, using the distance matrix created by the UPGMA algorithm.")
    wait_input()

    print("First, lets see a tree with a cut value of 10!")

    print("\n\t:::Step 5 - Creating Graph using UPGMA distance matrix and cut value of %d:::\n" % p.cut)
    # Creating Graph from distance matrix with the given cut value
    g = MyGraph.create_from_num_matrix(upgma.dists_mat, p.cut)
    take_while()
    g.print_graph_and_metrics()
    wait_input()

    print("Now lets see how the graph we get from cut 15 looks like!")
    wait_input()
    p.cut = 15
    g = MyGraph.create_from_num_matrix(upgma.dists_mat, p.cut)
    take_while()
    g.print_graph_and_metrics()
    wait_input()

    print("And what about cut 5?")
    wait_input()
    p.cut = 5
    g = MyGraph.create_from_num_matrix(upgma.dists_mat, p.cut)
    take_while()
    g.print_graph_and_metrics()
    wait_input()

    print("And that is it for our step-by-step demonstration!\n")
    

    print("\n-------------------------\n")

    print("To sum this presentation lets see the pipeline running autonomously, with the same parameters as the previous demonstration!\n")
    wait_input()

    print("First we create and setup the Pipeline:")
    print("             #fasta_query_filename          #fasta_database_filename       #cut #subst_matrix_filename        #gap penalty")
    print("p = Pipeline('src/tests/files/source.fasta', 'src/tests/files/seqdump.txt', 10, 'src/test/files/blosum62.mat', -8)")
    p = Pipeline('src/tests/files/source.fasta', 'src/tests/files/seqdump.txt', 10, 'src/test/files/blosum62.mat', -8)
    wait_input()

    print("And now we run our Pipeline with: p.execute().\nHere's its output:\n")
    p.execute()
    wait_input()

    print("And that is it!! This are this module functionalities, hope you enjoyed it!")
    wait_input()
