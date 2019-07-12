import tabulate



def print_cutflow(output):
    """Pretty-print cutflow data to the terminal."""
    for cutflow_name in [ x for x in output.keys() if x.startswith("cutflow") and x.endswith("j")]:
        if not len(output[cutflow_name]):
            continue
        table = []
        print("----")
        print(cutflow_name)
        print("----")
        for cut, count in sorted(output[cutflow_name].items(), key=lambda x:x[1], reverse=True):
            table.append([cut, count])
        print(tabulate.tabulate(table, headers=["Cut", "Passing events"]))
