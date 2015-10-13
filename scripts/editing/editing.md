# Summary of editing analysis


* Download the human hg19 editing positions from the RADAR database (v2, 52,618 positions).

* Prepare a VCF based on these data, after excluding intronic positions (8,126 positions left)

* Use that VCF and GATK ASEReadCounter to assess the level of editing at each site.

* Most of these positions show no editing whatsoever, so pick the 401 positions that show at least 5% editing in at least one sample.

* At these 401 positions, compute the median editing level for each sample. This is what is shown in the scatterplot.