import time
import re

base_def_id = None
results = []
line_count = 0

with open("def_corr.log") as f:
    for line in f:
        line_count += 1
        if line_count < 12:
            continue

        data = line.strip().split('\t')
        if len(data) != 3:
            continue

        time, def_id, def_num = data
        # str to list
        def_id_clean = re.sub(r'\s+', ' ', def_id).strip('[]')
        current_list = [list(map(int, re.findall(r'\d+', sublist))) for sublist in re.split(r'\],\s*\[', def_id_clean)]

        if base_def_id is None:
            base_def_id = current_list
            base_def_num = def_num
            results.append([time, 1])
            continue

        common_defs = sum(
            len(set(base_sublist) & set(current_sublist)) for base_sublist in base_def_id for current_sublist in
            current_list)
        corr_ratio = 1 - common_defs/int(base_def_num)
        results.append([time, corr_ratio])

with open("corr_ratio.log", 'w') as f:
    f.write('time\tcorr_ratio\n')
    for time, corr_ratio in results:
        f.write(f"{time}\t{corr_ratio}\n")



