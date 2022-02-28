#!/usr/bin/env python3

import math

def get_percentile(sorted_list, percentile):
    rank = int(math.ceil(percentile / 100.0 * len(sorted_list)))
    if rank == 0:
        return sorted_list[0]
    return sorted_list[rank - 1]

insert_sizes = []
with open('insertsize.sam', 'rt') as sam:
    for sam_line in sam:
        try:
            sam_parts = sam_line.split('\t')
            sam_flags = int(sam_parts[1])
            if sam_flags & 2:  # read mapped in proper pair
                insert_size = int(sam_parts[8])
                if insert_size > 0:
                    insert_sizes.append(insert_size)
        except (ValueError, IndexError):
            pass

insert_sizes = sorted(insert_sizes)
#insert_size_1st = get_percentile(insert_sizes, 1.0)
insert_size_99th = get_percentile(insert_sizes, 99.0)
print(insert_size_99th)
