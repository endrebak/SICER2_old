
from libcpp.vector cimport vector

cdef struct island:
    int start
    int end
    int chip_count
    int input_count


cdef class IslandVector:

    cdef vector[island] wrapped_vector

    cdef push_back(self, island i):
        self.wrapped_vector.push_back(i)

    def __str__(self):
        return "[" + ", ".join([str(i) for i in self.wrapped_vector]) + "]"

    def __repr__(self):
        return str(self)

    def __len__(self):
        return self.wrapped_vector.size()

    def __iter__(self):
        # slow, only implemented to ease testing
        return (v for v in self.wrapped_vector)



# def merge_nearby_bins(df, gaps_allowed, bin_size, score_threshold):

#     bin_size_minus_one = bin_size - 1

#     distance_allowed = (gaps_allowed * bin_size) + 2

#     if nrow(df) == 1:
#         return df

#     current_island = df[1, :]

#     merged_islands = similar(df, 0)

#     # just to have exact same result as SICER
#     slightly_less = score_threshold - 0.0000000001

#     for idx in 2:nrow(df):
#         dist = df[idx, :Start] - current_island[1, :End]
#         if dist <= distance_allowed:
#             current_island[:End] = df[idx, :End]
#             current_island[1, :Score] += df[idx, :Score]
#             current_island[1, :Count] += df[idx, :Count]
#             current_island[1, :InputCount] += df[idx, :InputCount]
#         else:
#             if current_island[1, :Score] > slightly_less:
#                 append!(merged_islands, current_island[:])
#             current_island = df[idx, :]

#     if current_island[1, :Score] > slightly_less:

#         append!(merged_islands, current_island[:])

#     delete!(merged_islands, 1)

def find_islands(bins_counts, int gaps_allowed, int bin_size, float score_threshold, int island_enriched_threshold):

    cdef:
        pass
        float slightly_less = score_threshold - 0.0000000001
        int i
        island _island
        IslandVector v

    chromosomes = bins_counts.keys()

    island_dict = dict()

    for chromosome in chromosomes:
        v = IslandVector()
        i = 0

        # TODO: if chromo not in chromsizes, remove
        bins, counts = bins_counts[chromosome]

        # TODO: need to iterate till you find first island > island_enriched_threshold

        for i in range(len(bins)):

            if counts[i] < island_enriched_threshold:
                continue

            _island = [bins[i], bins[i] + bin_size - 1, counts[i], 0]

            v.push_back(_island)


        # if bins.size() == 1:
        #     final_islands

        # for i in range(1, bins.size()):

        #     _bin = bins[i]

        # del bins_counts[chromosome]
