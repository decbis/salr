Sequence AnaLyzeR
=================

*Sequence AnaLyzeR*'s purpose is to analyze a sequence of tokens to evaluate how
random the sequence is.

The `pochisq` & related function are adapted from *Ent A Pseudorandom Number Sequence Test Program* (http://www.fourmilab.ch/random/).

Available analyzers
-------------------

character_count_tester
~~~~~~~~~~~~~~~~~~~~~~
Returns a (probability observed distribution is uniform, histogram) tuple, for each column calculated across all sequences in sequences_list.
It assumes all sequences have the same length.

Example::
    
    sequences_list = ['abc', 'bca']
    test_results = character_count_tester(sequences_list)
    for (probability, histogram) in test_results:
        print probability, ' the distribution is random'
        _chart(histogram)
