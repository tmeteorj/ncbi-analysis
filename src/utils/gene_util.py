def get_opposite_dna(dna):
    op_dna = ''
    for x in dna:
        if x == 'a': op_dna += 't'
        if x == 't': op_dna += 'a'
        if x == 'c': op_dna += 'g'
        if x == 'g': op_dna += 'c'
    return op_dna
