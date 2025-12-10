from mooonpy.tools.file_utils import Path

def read_lunarlogfile(file: [Path, str], keywords: [str, list]):
    keywords = set(keywords)

    log_data = {}
    with open(file) as f:
        for line in f:
            key, _, rest = line.partition(":")
            key = key.strip()
            if key in keywords:
                value = float(rest.split()[0])
                log_data[key] = value

    return log_data