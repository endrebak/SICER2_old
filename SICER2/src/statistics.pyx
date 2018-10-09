
def compute_background_probabilities(total_chip_count, bin_size, effective_genome_fraction, gap_intervals_allowed):

    average_window_readcount = total_chip_count * (bin_size / effective_genome_fraction)
    tag_density = total_chip_count / effective_genome_fraction


    # self.tag_density = total_tags * 1.0 / genomeLength;
    # self.average = self.tag_density * windowSize;


    poisson = Poisson(average_window_readcount)

    island_enriched_threshold = compute_enriched_threshold(poisson)

    gap_contribution = compute_gap_factor(island_enriched_threshold, gap_intervals_allowed, poisson)

    boundary_contribution = compute_boundary(island_enriched_threshold, gap_intervals_allowed, poisson)

    genome_length_in_bins = effective_genome_fraction / bin_size

    score_threshold = compute_score_threshold(average_window_readcount, island_enriched_threshold,
                                              gap_contribution, boundary_contribution, genome_length_in_bins, bin_size)


    return score_threshold, island_enriched_threshold, average_window_readcount
