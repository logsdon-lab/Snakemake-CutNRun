from os.path import join, dirname


CONDITIONS = ["treatment", "control"]


def load_primer_list() -> dict[str, str]:
    primers = {}
    primer_list = config["primer_list"]
    with open(primer_list, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            pset, index, seq = line.strip().split(",")
            primers[index] = seq
    
    assert len(primers) != 0, f"No primers loaded from {primer_list}."
    return primers


def get_primer_pair_arg(primer_pair: list[str], primer_map: dict[str, str]) -> str:
    try:
        primer_1_idx, primer_2_idx = primer_pair
    except Exception:
        raise ValueError(f"Two primers expected. Received: {primer_pair}")

    try:
        primer_1_seq = primer_map[primer_1_idx]
        primer_2_seq = primer_map[primer_2_idx]
    except KeyError:
        raise ValueError(f"Either {primer_1_idx} or {primer_2_idx} not in primer list.")
    
    return f"-a {primer_1_seq}...{primer_2_seq}"


def get_cutadapt_arg(wc, primer_map: dict[str, str]) -> str:
    # https://cutadapt.readthedocs.io/en/stable/guide.html#linked-adapters
    # Will be linked
    sm_config = config["samples"]
    primers: list[list[str]] | list = sm_config[wc.sm]["data"][wc.condition]["primers"]
    if not primers:
        raise ValueError("Must have at least on prime pair.")
    
    args = []
    if isinstance(primers[0], list):
        for primer_pair in primers:
            primer_pair_arg = get_primer_pair_arg(primer_pair, primer_map)
            args.append(primer_pair_arg)
    else:
        args.append(get_primer_pair_arg(primers, primer_map))

    return " ".join(args) 
