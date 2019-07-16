import tabulate

def print_cutflow(output, outfile=None):
    """Pretty-print cutflow data to the terminal."""
    for i, cutflow_name in enumerate([ x for x in output.keys() if x.startswith("cutflow") ]):
        if not len(output[cutflow_name]):
            continue
        table = []
        print("----")
        print(cutflow_name)
        print("----")
        for cut, count in sorted(output[cutflow_name].items(), key=lambda x:x[1], reverse=True):
            table.append([cut, count])
        text = tabulate.tabulate(table, headers=["Cut", "Passing events"])
        print(text)
        if outfile:
            with open(outfile, "w" if not i else "a") as f:
                f.write(rf'----\n{cutflow_name}\n----')
                f.write(text + r"\n")
