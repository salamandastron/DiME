
This program is an implementation of the method described in the paper "Statistical significance of communities in network", written by Andrea Lancichinetti, Filippo Radicchi and José J. Ramasco. Each feedback is very welcome. If you  find a bug or have problems, or want to give advises, please contact us:

andrea.lancichinetti@isi.it
radicchi@isi.it
jramasco@isi.it


Turin, 9 Oct 2009
---------------------------------------------------------------





-------------------- How to compile -----------------------------------
In order to compile, type:

make

PLEASE NOTE THAT YOU NEED TO INSTALL THE GSL LIBRARY.
http://www.gnu.org/software/gsl/

-------------------- How to run the program ---------------------------

To run the program, type:

./compare [OPTIONS]



MANDATORY OPTIONS:

-f <filename>
                This option is used to specify the file where the network data are stored:
                filename should include a list of edges and their weights (optionally). Only integers number are allowed. This is the format:
                source_node target_node (weight)
                Indeed, since the program works only on undirected graphs, it does not care about which node is the source and which one is the target. Repetitions and self-loops are neglected
-c <filename2>
                filename2 is the community file, i.e. it must contain the groups:
                each row should contain the nodes belonging to the same community separated by a tab or a blank space.


OTHER OPTIONS:

-t <number> :
                set the threshold equal to number
-cscore :
                the program uses the original definition of the c-score: this is slower than the default definition which gives comparable results
-nobcore :
                the program doesn't search for the b-q-core. this can be used if you want a fast computation of the c-cores
-list :
                when using this option the community file is expected to be a list of "integer integer" which are considered node-membership
-v :
                this flags increments the verbosity


-------------------- Examples ---------------------------
        
./compare -f network.dat -c community.dat -t 0.1
./compare -f network.dat -c community_list.dat -list -v



Have a look at the files network.dat and community.dat (community_list.dat) to see an example.


-------------------- Output ---------------------------

The program will compute the worst node in the cluster, the c-score, the b-score, the c-core and the b-core with a given threshold (default value is 5%).
It will also produce five output files.


[network file].table
		This file contains a summary of the results of the program. For each group, you can find six numbers:
		the group label, the number of nodes in the group, the c-score, the b-score, the number of nodes in the c-core and the number of nodes in the b-core
c_core.dat
		For each group, you can find the nodes in the c-core.
c_erased.dat
		For each group, you can find the nodes which are NOT in the c-core.
b_core.dat, b_erased.dat
		These two files are the anologous for the b-score.


-------------------- Seed for the random number generator ---------------------------
 
In the file time_seed.dat you can edit the seed which generates the random numbers. After reading, the program will increase this number by 1 (this is done to generate different networks running the program again and again). If the file is erased, it will be produced by the program again.


