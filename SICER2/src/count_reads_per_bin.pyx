import numpy as np

cpdef count_reads_per_bin(tags):

    cdef:
        cdef long[::1] bins
        cdef long[::1] counts
        cdef long[::1] v
        cdef int i
        cdef int last
        cdef int current
        cdef int count = 0
        cdef int nparsed

    bins_counts = dict()
    for k, v in tags.items():
        nparsed = 0

        bin_arr = np.zeros(len(v), dtype=int)
        bins = bin_arr
        count_arr = np.zeros(len(v), dtype=int)
        counts = count_arr

        if len(v) >= 1:
            last = bins[0]

        for i in range(len(v)):
            current = v[i]

            if current != last:
                bins[nparsed] = current
                counts[nparsed] = count
                last = current
                count = 1
                nparsed += 1
            else:
                count += 1


        bins_counts[k] = (bins, counts)

    return bins_counts[:nparsed]
