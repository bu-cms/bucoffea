import tabulate

def print_cutflow(output, outfile=None):
    """Pretty-print cutflow data to the terminal."""
    for i, cutflow_name in enumerate([ x for x in output.keys() if x.startswith("cutflow") ]):
        for dataset in output[cutflow_name].keys():
            if not len(output[cutflow_name][dataset]):
                continue
            table = []
            print("----")
            print(cutflow_name, dataset)
            print("----")
            for cut, count in sorted(output[cutflow_name][dataset].items(), key=lambda x:x[1], reverse=True):
                table.append([cut, count])
            text = tabulate.tabulate(table, headers=["Cut", "Passing events"], floatfmt=".1f")
            print(text)
            if outfile:
                with open(outfile, "w" if not i else "a") as f:
                    f.write(f'----\n{cutflow_name, dataset}\n----')
                    f.write(text + "\n")
