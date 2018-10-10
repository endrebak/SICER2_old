
import pytest

from SICER2.src.reads_to_bins import add_reads_to_dict, files_to_bin_counts


args = {"bin_size": 200, "fragment_size": 150, "drop_duplicates": True}


def test_files_to_bins():

    d = files_to_bin_counts(["tests/test.bed", "tests/test.bed"], args)
    print(d)

    print(d["chr7", "+"][0])
    assert list(d["chr7", "+"][0]) == [20246600, 91135000]


def test_reads_to_dict():

    d = add_reads_to_dict("tests/test.bed")
    print(d)

    assert list(d["chr7", "+"]) == [20246668, 20246668, 91135110]

	# Command being timed: "bin/SICER2 -t /mnt/scratch/projects/epipp/data/bam/Exp1_0h_H3K27me3.bed /mnt/scratch/projects/epipp/data/bam/Exp1_3h_H3K27me3.bed /mnt/scratch/projects/epipp/data/bam/Exp1_6h_H3K27me3.bed /mnt/scratch/projects/epipp/data/bam/Exp1_9h_H3K27me3.bed /mnt/scratch/projects/epipp/data/bam/Exp1_12h_H3K27me3.bed /mnt/scratch/projects/epipp/data/bam/Exp1_15h_H3K27me3.bed /mnt/scratch/projects/epipp/data/bam/Exp1_18h_H3K27me3.bed /mnt/scratch/projects/epipp/data/bam/Exp1_21h_H3K27me3.bed /mnt/scratch/projects/epipp/data/bam/Exp1_24h_H3K27me3.bed /mnt/scratch/projects/epipp/data/bam/Exp1_Unsynchronized_H3K27me3.bed /mnt/scratch/projects/epipp/data/bam/Exp2_0h_H3K27me3.bed /mnt/scratch/projects/epipp/data/bam/Exp2_3h_H3K27me3.bed /mnt/scratch/projects/epipp/data/bam/Exp2_6h_H3K27me3.bed /mnt/scratch/projects/epipp/data/bam/Exp2_9h_H3K27me3.bed /mnt/scratch/projects/epipp/data/bam/Exp2_12h_H3K27me3.bed /mnt/scratch/projects/epipp/data/bam/Exp2_15h_H3K27me3.bed /mnt/scratch/projects/epipp/data/bam/Exp2_18h_H3K27me3.bed /mnt/scratch/projects/epipp/data/bam/Exp2_21h_H3K27me3.bed /mnt/scratch/projects/epipp/data/bam/Exp2_24h_H3K27me3.bed /mnt/scratch/projects/epipp/data/bam/Exp2_Unsynchronized_H3K27me3.bed -c /mnt/scratch/projects/epipp/data/bam/Exp1_0h_Input.bed /mnt/scratch/projects/epipp/data/bam/Exp1_3h_Input.bed /mnt/scratch/projects/epipp/data/bam/Exp1_6h_Input.bed /mnt/scratch/projects/epipp/data/bam/Exp1_9h_Input.bed /mnt/scratch/projects/epipp/data/bam/Exp1_12h_Input.bed /mnt/scratch/projects/epipp/data/bam/Exp1_15h_Input.bed /mnt/scratch/projects/epipp/data/bam/Exp1_18h_Input.bed /mnt/scratch/projects/epipp/data/bam/Exp1_21h_Input.bed /mnt/scratch/projects/epipp/data/bam/Exp1_24h_Input.bed /mnt/scratch/projects/epipp/data/bam/Exp1_Unsynchronized_Input.bed /mnt/scratch/projects/epipp/data/bam/Exp2_0h_Input.bed /mnt/scratch/projects/epipp/data/bam/Exp2_3h_Input.bed /mnt/scratch/projects/epipp/data/bam/Exp2_6h_Input.bed /mnt/scratch/projects/epipp/data/bam/Exp2_9h_Input.bed /mnt/scratch/projects/epipp/data/bam/Exp2_12h_Input.bed /mnt/scratch/projects/epipp/data/bam/Exp2_15h_Input.bed /mnt/scratch/projects/epipp/data/bam/Exp2_18h_Input.bed /mnt/scratch/projects/epipp/data/bam/Exp2_21h_Input.bed /mnt/scratch/projects/epipp/data/bam/Exp2_24h_Input.bed /mnt/scratch/projects/epipp/data/bam/Exp2_Unsynchronized_Input.bed"
	# User time (seconds): 573.16
	# System time (seconds): 32.99
	# Percent of CPU this job got: 653%
	# Elapsed (wall clock) time (h:mm:ss or m:ss): 1:32.77
	# Average shared text size (kbytes): 0
	# Average unshared data size (kbytes): 0
	# Average stack size (kbytes): 0
	# Average total size (kbytes): 0
	# Maximum resident set size (kbytes): 1705780
	# Average resident set size (kbytes): 0
	# Major (requiring I/O) page faults: 0
	# Minor (reclaiming a frame) page faults: 188591
	# Voluntary context switches: 1156
	# Involuntary context switches: 128254
	# Swaps: 0
	# File system inputs: 0
	# File system outputs: 0
	# Socket messages sent: 0
	# Socket messages received: 0
	# Signals delivered: 0
	# Page size (bytes): 4096
	# Exit status: 0
