from collections import defaultdict
import sys
import numpy as np

from libc.stdio cimport FILE, fopen, fwrite, fscanf, fclose, fseek, SEEK_END, ftell, stdout, stderr, EOF

from cython.operator import dereference
from libcpp.algorithm cimport sort as stdsort
from libcpp.map cimport map as cppmap
from libcpp.algorithm cimport unique
from libcpp.vector cimport vector
from libcpp.string cimport string

cdef extern from "<algorithm>" namespace "std" nogil:
    OutputIter merge[InputIter1, InputIter2, OutputIter] (InputIter1 first1, InputIter1 last1,
                                                          InputIter2 first2, InputIter2 last2,
                                                          OutputIter result)

cimport cython

import numpy as np

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.initializedcheck(False)
cpdef count_reads_per_bin(tags):

    # print("tags", tags)
    cdef:
        long[::1] bins
        long[::1] counts
        Vector v
        int i
        int last = -1
        int current
        int count = 0
        int nparsed
        int last_two_equal

    bins_counts = dict()
    for k, v in tags.items():
        nparsed = 0
        count = 0

        # print(k)
        bin_arr = np.ones(len(v), dtype=int) * -1
        bins = bin_arr
        count_arr = np.ones(len(v), dtype=int)
        counts = count_arr

        # print(v)
        if len(v) >= 1:
            last = v.wrapped_vector[0]
            # bins[0] = last
            # counts[0] = 1

        for i in range(0, len(v)):
            current = v.wrapped_vector[i]
            # print("current", current)
            # print("last", last)
            # print("count", count)

            if current != last:
                # print("adding", last, count)
                bins[nparsed] = last
                counts[nparsed] = count
                last = current
                count = 1
                nparsed += 1
            else:
                count += 1

        # print("v", v)
        # print(bin_arr)
        # print(count_arr)
        last_two_equal = v.wrapped_vector[len(v) - 2] == v.wrapped_vector[len(v) - 1]
        if last_two_equal:
            bins[nparsed] = v.wrapped_vector[len(v) - 1]
            counts[nparsed] = count
            nparsed += 1
        else:
            # print("else")
            # print("adding", v.wrapped_vector[len(v) - 1], 1)
            bins[nparsed] = v.wrapped_vector[len(v) - 1]
            counts[nparsed] = 1
            nparsed += 1

        # print(bin_arr)
        # print(count_arr)
        bins_counts[k] = (bin_arr[:nparsed], count_arr[:nparsed])

    return bins_counts


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.initializedcheck(False)
cdef class Vector:

    cdef vector[int] wrapped_vector

    cdef push_back(self, int num):
        self.wrapped_vector.push_back(num)

    def sort(self):
        stdsort(self.wrapped_vector.begin(), self.wrapped_vector.end())

    def unique(self):
        self.wrapped_vector.erase(unique(self.wrapped_vector.begin(), self.wrapped_vector.end()), self.wrapped_vector.end())

    def __str__(self):
        return "[" + ", ".join([str(i) for i in self.wrapped_vector]) + "]"

    def __repr__(self):
        return str(self)

    def __len__(self):
        return self.wrapped_vector.size()

    def __iter__(self):
        # slow, only implemented to ease testing
        return (v for v in self.wrapped_vector)

    def merge(self, Vector other):

        cdef vector[int] o = vector[int](len(self) + len(other))
        merge(self.wrapped_vector.begin(), self.wrapped_vector.end(),
              other.wrapped_vector.begin(), other.wrapped_vector.end(),
              o.begin())

        cdef Vector output = Vector()
        output.wrapped_vector = o

        return output


@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.initializedcheck(False)
cpdef files_to_bin_counts(files, args):

    cdef:
        int bin_size = args["bin_size"]
        int half_fragment_size = args["fragment_size"] / 2
        Vector v
        Vector v2
        cdef long[::1] bin_arr


    sum_tags = defaultdict(list)
    sys.stderr.write("Parsing file:\n")
    sys.stderr.flush()
    for f in files:

        sys.stderr.write("  " + f + "\n")
        sys.stderr.flush()

        tags = add_reads_to_dict(f)
        # print("tags", tags)

        for v in tags.values():
            v.sort()

        if args["drop_duplicates"]:
            for k, v in tags.items():
                v.unique()

        for (chromosome, strand), v in tags.items():

            if strand == "+":
                for i in range(len(v)):
                    v.wrapped_vector[i] = v.wrapped_vector[i] + half_fragment_size
            else:
                for i in range(len(v)):
                    v.wrapped_vector[i] = v.wrapped_vector[i] - half_fragment_size

            for i in range(len(v)):
                v.wrapped_vector[i] = v.wrapped_vector[i] - (v.wrapped_vector[i] % bin_size)

            if chromosome not in sum_tags:
                sum_tags[chromosome] = v
            else:
                v2 = sum_tags[chromosome]
                sum_tags[chromosome] = v.merge(v2)

    sys.stderr.write("\nCounting the number of reads in each bin\n\n")
    sys.stderr.flush()

    bins_counts = count_reads_per_bin(sum_tags)

    return bins_counts



cpdef add_reads_to_dict(f):

    genome = dict()
    cdef Vector v

    for line in open(f):
        chromosome, left, right, _, _, strand = line.split()

        if strand == "+":
            five_end = int(left)
        else:
            five_end = int(right)

        if (chromosome, strand) in genome:
            v = genome[chromosome, strand]
            v.wrapped_vector.push_back(five_end)
        else:
            v = Vector()
            v.wrapped_vector.push_back(five_end)
            genome[chromosome, strand] = v

    return genome



######## this was actually slower!

# cpdef add_reads_to_dict(f):

#     genome = dict()
#     cdef:
#         Vector v
#         FILE *f_handle
#         char chromosome [10]
#         char strand [2]
#         int left
#         int right

#     fp = fopen(f.encode(), "r")

#     while (
#             fscanf(fp, "%s\t%d\t%d\t%*s\t%*d\t%s\n", chromosome, &left, &right, strand) != EOF
#     ):


#         # print("----")
#         # print("chromosome is", chromosome)
#         # print("strand is", strand)
#         # print("left is", left)
#         # print("right is", right)
#         if strand == b"+":
#             five_end = left
#         else:
#             five_end = right

#         # print("five end is ", five_end)

#         if (chromosome, strand) in genome:
#             v = genome[chromosome, strand]
#             v.wrapped_vector.push_back(five_end)
#         else:
#             v = Vector()
#             v.wrapped_vector.push_back(five_end)
#             genome[chromosome, strand] = v

#     return genome
