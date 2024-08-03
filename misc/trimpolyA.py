import re
import math
import gzip
import sys

#note that the intervals in this script are 0-based.

def find_seeds(seq:str, sub:str) -> list:
    start = 0
    seed_indexes = list()
    while True:
        start = seq.find(sub,start)
        if start >= 0:
            seed_indexes.append([start,start+4])
            start += 1
        else:
            return seed_indexes


def seed_merge(seed_indexes:list, gap_threshold:int = 2) -> list:
    merged = []
    for seed_index in seed_indexes:
        if not merged or merged[-1][1] + gap_threshold < seed_index[0]:
            merged.append(seed_index)
        else:
            merged[-1][1] = max(seed_index[1],merged[-1][1])

    return merged
    

def seed_extention(seq:str, merged:list, match_score:int = 1, match_extention_score:int = 1.1, mismatch_open_penalty:int = 1, mismatch_extention_penalty:int = 2,  penalty_buffer:int = -5,base = "A") -> list:
    '''
    this function is for rescuing As that are omitted during the find seed session.
    '''
    extended_intervals = list()
    for interval in merged:
        score = 0
        start = interval[0]
        end = interval[1]
        extended_start = start
        extended_end = end
        
        while score > penalty_buffer and start > 0:
            prev = seq[start]
            start -= 1
            current = seq[start]
            if current == base:
                if current != prev:
                    score += match_score
                    if score > 0:
                        extended_start = start
                        score = 0
                else:
                    score += match_extention_score
                    if score > 0:
                        extended_start = start
                        score = 0

            else:
                if prev != base:
                    score -= mismatch_open_penalty
                else:
                    score -= mismatch_extention_penalty
                   
        score = 0
        while score > penalty_buffer and end < len(seq):
            prev = seq[end-1]
            end += 1
            current = seq[end-1]
            if current == base:
                if current != prev:
                    score += match_score
                    if score > 0:
                        extended_end = end
                        score = 0

                else:
                    score += match_extention_score
                    if score > 0:
                        extended_end = end
                        score = 0

            else:
                if prev == base:
                    score -= mismatch_open_penalty
                else:
                    score -= mismatch_extention_penalty

        extended_intervals.append([extended_start,extended_end])

    return extended_intervals

def segment_bridge(segment_indexes:list, gap_threshold:int = 4):
    segment_indexes = seed_merge(segment_indexes,gap_threshold)
    return segment_indexes
    
def main():
    f1 = gzip.open(sys.argv[1],"rt")
    i = 1
    f = open(sys.argv[2],"w")
    for line in f1:
        if i % 4 == 2:
            line = line.strip()
            A_max = [0,0,0,0]
            T_max = [0,0,0,0]
            Aseed_indexes = segment_bridge(seed_extention(line,seed_merge(find_seeds(line,"AAAA")),base = "A"))
            Tseed_indexes = segment_bridge(seed_extention(line,seed_merge(find_seeds(line,"TTTT")),base = "T"))
            for Aseed_index in Aseed_indexes:
                start = int(Aseed_index[0])
                end = int(Aseed_index[1])
                A_num = line[start:end].count("A")
                position = int(Aseed_index[1])/len(line)
                score = A_num * (1-math.sqrt(1-(position)**2))
                if score > A_max[2]:
                    A_max = [start,end,score,A_num]

            for Tseed_index in Tseed_indexes:
                start = int(Tseed_index[0])
                end = int(Tseed_index[1])
                T_num = line[start:end].count("T")
                position = int(Tseed_index[0])/len(line)
                score = T_num * (1-math.sqrt(1-(position-1)**2))
                if score > T_max[2]:
                    T_max = [start,end,score,T_num]
            
            if T_max[2] == A_max[2] and T_max[2] == 0:
                f.write("%s\n"%line)
                i += 1
                continue

            if T_max[2] < A_max[2]:
                if A_max[3] >= 10:
                    f.write("%s\n" %line[:A_max[0]])
                else:
                    if T_max[3] < 10:
                        f.write("%s\n" %line)
                    else:
                        f.write("%s\n" %line[T_max[1]:])
            else:
                if T_max[3] >= 10:
                    f.write("%s\n" %line[T_max[1]:])
                else:
                    if A_max[3] < 10:
                        f.write("%s\n" %line)
                    else:
                        f.write("%s\n" %line[:A_max[0]])
        elif i % 4 == 1:
            line = line.replace("@",">")
            f.write(line)

        i += 1
if __name__ == "__main__":    
    main()