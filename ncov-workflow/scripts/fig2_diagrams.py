import argparse
import json

#make diagrams coloring branches of the tree used to calculate dn/ds accumulation in a time window for Fig 2

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Make diagrams for Figure 2 (dn/ds accumulation)",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--input', type=str, metavar="JSON", required=True, help="input Auspice JSON")
    parser.add_argument('--output', type=str, metavar="JSON", required=True, help="output Auspice JSON")
    parser.add_argument('--nodes-in-window', type=str, metavar="JSON", required=True, help="input JSON with names of nodes in the time window")
    args = parser.parse_args()

    with open(args.input, "r") as f:
        auspice_json = json.load(f)

    with open(args.nodes_in_window, "r") as f:
        nodes_json = json.load(f)

    nodes_in_window = nodes_json['nodes']


    # check whether each node is in time window, and add attribute that can be used for coloring
    def attach_labels(n):

        n["node_attrs"]["in_window"] = {}

        if n["name"] in nodes_in_window:
            n["node_attrs"]["in_window"]["value"] = "True"
        else:
            n["node_attrs"]["in_window"]["value"] = "False"

        if "children" in n:
            for c in n["children"]:
                attach_labels(c)

    attach_labels(auspice_json["tree"])

    # add 'colorings' option for in_window
    auspice_json['meta']['colorings'] += [{"key": "in_window", "scale": [["True", "#DB2823"], ["False", "#f6f6f6"]], "title": "in_window", "type": "categorical"}]


    with open(args.output, 'w') as f:
        json.dump(auspice_json, f, indent=2)
