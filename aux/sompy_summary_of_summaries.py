import sys
assert(sys.version_info.major == 3)


# ++++++++++++++
# |   FUNC     |
# ++++++++++++++
def checkHeader(header: str) -> bool:
    expected_header_fields = ["idx","type","total.truth","total.query","tp","fp","fn","unk","ambi","recall","recall_lower","recall_upper","recall2","precision","precision_lower","precision_upper","na","ambiguous","fp.region.size","fp.rate","sompyversion","sompycmd"]
    return header.rstrip().split(',') == expected_header_fields


def header4outputfile(header: str) -> str:
    f = lambda s: s+".avg"
    l_old_header =  header.rstrip().split(',')
    l_new_header = l_old_header[:2] + list(map(f, l_old_header[2:-2])) + ["F1-score.avg"] + l_old_header[-2:]
    return ','.join(l_new_header)


def f1_score(recall: float, precision:float) -> float:
    return (2*precision*recall) / (precision + recall)


# ++++++++++++++
# |   MAIN     |
# ++++++++++++++
found_header = False
header       = ""

indel_lines_avg = []
snv_lines_avg   = []

nb_indel_lines = 0
nb_snv_lines   = 0

for line in sys.stdin:
    if line.startswith('idx', 0, 3) and not found_header:
        header = line
        assert(checkHeader(header))
        found_header = True

    elif line[0] == "0":      # indel line
        current_line_buffer = line.split(",")[2:-2]
        if not indel_lines_avg:   # if this is the first observed indel line
            indel_lines_avg = current_line_buffer
            nb_indel_lines += 1
        else:
            assert(len(indel_lines_avg)==len(current_line_buffer))
            # add values of current line to avg
            indel_lines_avg = [float(avg)*(nb_indel_lines/(nb_indel_lines+1)) + float(add)*(1/(nb_indel_lines+1)) for avg,add in zip(indel_lines_avg,current_line_buffer)]
            nb_indel_lines += 1

    elif line[0] == "1":      # snv line
        current_line_buffer = line.split(",")[2:-2]
        if not snv_lines_avg:   # if this is the first observed snv line
            snv_lines_avg = current_line_buffer
            nb_snv_lines += 1
        else:
            assert(len(snv_lines_avg)==len(current_line_buffer))
            # add values of current line to avg
            snv_lines_avg = [float(avg)*(nb_snv_lines/(nb_snv_lines+1)) + float(add)*(1/(nb_snv_lines+1)) for avg,add in zip(snv_lines_avg,current_line_buffer)]
            nb_snv_lines += 1

    else:
        continue


# ++++++++++++++
# |  OUTPUT    |
# ++++++++++++++
assert(header.rstrip().split(",")[9] == "recall")
recall_field_idx = 9
assert(header.rstrip().split(",")[13] == "precision")
precision_field_idx = 13
print(header4outputfile(header))
print("0,indels," + ','.join(str(x) for x in indel_lines_avg) + ',' + str(f1_score(indel_lines_avg[recall_field_idx], indel_lines_avg[precision_field_idx])) + ",-,-")
print("1,SNVs,"   + ','.join(str(x) for x in snv_lines_avg)   + ',' + str(f1_score(snv_lines_avg[recall_field_idx], snv_lines_avg[precision_field_idx])) + ",-,-")