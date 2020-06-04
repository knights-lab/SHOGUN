/usr/bin/time -v shogun assign_taxonomy --database /home/bhillmann/petard/data/rep94 --aligner bowtie2 --no-capitalist --input /home/bhillmann/petard/data/shogun_mem/BIOM/89437/alignment.bowtie2.100k.sam --output scratch/profile.tsv
        Command being timed: "shogun assign_taxonomy --database /home/bhillmann/petard/data/rep94 --aligner bowtie2 --no-capitalist --input /home/bhillmann/petard/data/shogun_mem/BIOM/89437/alignment.bowtie2.100k.sam --output scratch/profile.tsv"
        User time (seconds): 0.67
        System time (seconds): 0.04
        Percent of CPU this job got: 88%
        Elapsed (wall clock) time (h:mm:ss or m:ss): 0:00.81
        Average shared text size (kbytes): 0
        Average unshared data size (kbytes): 0
        Average stack size (kbytes): 0
        Average total size (kbytes): 0
        Maximum resident set size (kbytes): 84964
        Average resident set size (kbytes): 0
        Major (requiring I/O) page faults: 122
        Minor (reclaiming a frame) page faults: 16865
        Voluntary context switches: 424
        Involuntary context switches: 6
        Swaps: 0
        File system inputs: 46464
        File system outputs: 80
        Socket messages sent: 0
        Socket messages received: 0
        Signals delivered: 0
        Page size (bytes): 4096
        Exit status: 0
        
(shogun_v107)
20.05.14 4:08:27 : bhillmann@x86_64-conda_cos6-linux-gnu:/home/bhillmann/petard/data/shogun_mem/BIOM/89437 (master)
/usr/bin/time -v shogun assign_taxonomy --database /home/bhillmann/petard/data/rep94 --aligner bowtie2 --no-capitalist --input /home/bhillmann/petard/data/shogun_mem/BIOM/89437/alignment.bowtie2.sam --output scratch/profile.tsv
        Command being timed: "shogun assign_taxonomy --database /home/bhillmann/petard/data/rep94 --aligner bowtie2 --no-capitalist --input /home/bhillmann/petard/data/shogun_mem/BIOM/89437/alignment.bowtie2.sam --output scratch/profile.tsv"
        User time (seconds): 216.58
        System time (seconds): 8.98
        Percent of CPU this job got: 99%
        Elapsed (wall clock) time (h:mm:ss or m:ss): 3:45.87
        Average shared text size (kbytes): 0
        Average unshared data size (kbytes): 0
        Average stack size (kbytes): 0
        Average total size (kbytes): 0
        Maximum resident set size (kbytes): 7,097,516
        Average resident set size (kbytes): 0
        Major (requiring I/O) page faults: 0
        Minor (reclaiming a frame) page faults: 2424534
        Voluntary context switches: 3
        Involuntary context switches: 2108
        Swaps: 0
        File system inputs: 0
        File system outputs: 424
        Socket messages sent: 0
        Socket messages received: 0
        Signals delivered: 0
        Page size (bytes): 4096
        Exit status: 0        

(shogun)
20.05.14 4:14:27 : bhillmann@x86_64-conda_cos6-linux-gnu:/home/bhillmann/petard/data/shogun_mem/BIOM/89437 (master)
/usr/bin/time -v shogun assign_taxonomy --database /home/bhillmann/petard/data/rep94 --aligner bowtie2 --no-capitalist --input /home/bhillmann/petard/data/shogun_mem/BIOM/89437/alignment.bowtie2.sam --output scratch/profile.tsv
        Command being timed: "shogun assign_taxonomy --database /home/bhillmann/petard/data/rep94 --aligner bowtie2 --no-capitalist --input /home/bhillmann/petard/data/shogun_mem/BIOM/89437/alignment.bowtie2.sam --output scratch/profile.tsv"
        User time (seconds): 217.41
        System time (seconds): 8.15
        Percent of CPU this job got: 99%
        Elapsed (wall clock) time (h:mm:ss or m:ss): 3:45.83
        Average shared text size (kbytes): 0
        Average unshared data size (kbytes): 0
        Average stack size (kbytes): 0
        Average total size (kbytes): 0
        Maximum resident set size (kbytes): 7097384
        Average resident set size (kbytes): 0
        Major (requiring I/O) page faults: 0
        Minor (reclaiming a frame) page faults: 2380710
        Voluntary context switches: 1
        Involuntary context switches: 560
        Swaps: 0
        File system inputs: 0
        File system outputs: 424
        Socket messages sent: 0
        Socket messages received: 0
        Signals delivered: 0
        Page size (bytes): 4096
        Exit status: 0

# Sorting the file
/usr/bin/time -v sort -t $'\t' -k 1,1 -n BIOM/89437/alignment.bowtie2.sam > sort.sam
	Command being timed: "sort -t 	 -k 1,1 -n BIOM/89437/alignment.bowtie2.sam"
	User time (seconds): 731.02
	System time (seconds): 42.83
	Percent of CPU this job got: 487%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 2:38.81
	Average shared text size (kbytes): 0
	Average unshared data size (kbytes): 0
	Average stack size (kbytes): 0
	Average total size (kbytes): 0
	Maximum resident set size (kbytes): 43131592
	Average resident set size (kbytes): 0
	Major (requiring I/O) page faults: 0
	Minor (reclaiming a frame) page faults: 10782498
	Voluntary context switches: 7424
	Involuntary context switches: 53039
	Swaps: 0
	File system inputs: 66757336
	File system outputs: 66757320
	Socket messages sent: 0
	Socket messages received: 0
	Signals delivered: 0
	Page size (bytes): 4096
	Exit status: 0

/usr/bin/time -v shogun assign_taxonomy --database /home/bhillmann/petard/data/rep94 --aligner bowtie2 --no-capitalist --input sort.sam --output profile.tsv
	Command being timed: "shogun assign_taxonomy --database /home/bhillmann/petard/data/rep94 --aligner bowtie2 --no-capitalist --input sort.sam --output profile.tsv"
	User time (seconds): 404.95
	System time (seconds): 6.63
	Percent of CPU this job got: 99%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 6:52.06
	Average shared text size (kbytes): 0
	Average unshared data size (kbytes): 0
	Average stack size (kbytes): 0
	Average total size (kbytes): 0
	Maximum resident set size (kbytes): 83584
	Average resident set size (kbytes): 0
	Major (requiring I/O) page faults: 0
	Minor (reclaiming a frame) page faults: 21026
	Voluntary context switches: 17
	Involuntary context switches: 3029
	Swaps: 0
	File system inputs: 0
	File system outputs: 600
	Socket messages sent: 0
	Socket messages received: 0
	Signals delivered: 0
	Page size (bytes): 4096
	Exit status: 0

# after implementing LCA first pass
/usr/bin/time -v shogun assign_taxonomy --database /home/bhillmann/petard/data/rep94 --aligner bowtie2 --no-capitalist --input sort.sam --output profile.tsv
	Command being timed: "shogun assign_taxonomy --database /home/bhillmann/petard/data/rep94 --aligner bowtie2 --no-capitalist --input sort.sam --output profile.tsv"
	User time (seconds): 185.54
	System time (seconds): 6.06
	Percent of CPU this job got: 99%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 3:12.61
	Average shared text size (kbytes): 0
	Average unshared data size (kbytes): 0
	Average stack size (kbytes): 0
	Average total size (kbytes): 0
	Maximum resident set size (kbytes): 90840
	Average resident set size (kbytes): 0
	Major (requiring I/O) page faults: 0
	Minor (reclaiming a frame) page faults: 22889
	Voluntary context switches: 17
	Involuntary context switches: 11228
	Swaps: 0
	File system inputs: 0
	File system outputs: 624
	Socket messages sent: 0
	Socket messages received: 0
	Signals delivered: 0
	Page size (bytes): 4096
	Exit status: 0
