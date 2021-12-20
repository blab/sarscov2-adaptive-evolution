import argparse
import json

#map emerging_lineage designation onto the lineages I want to color for figures in the manuscript

manuscript_lineage_names = {'unassigned': 'basal', 'A.23.1': 'other VOI', 'B.1.1.7 (Alpha)': 'Alpha', 'B.1.1.318': 'other VOI', 'B.1.1.519': 'other VOI', 'B.1.351 (Beta)': 'Beta', 'B.1.427/429 (Epsilon)': 'other VOI', 'B.1.525 (Eta)': 'other VOI', 'B.1.526 (Iota)': 'other VOI', 'B.1.617.1 (Kappa)': 'other VOI', 'B.1.617.2 (Delta)': 'Delta', 'B.1.619': 'other VOI', 'B.1.620': 'other VOI', 'B.1.621': 'other VOI', 'C.36.3': 'other VOI', 'C.37 (Lambda)': 'other VOI', 'P.1 (Gamma)': 'Gamma', 'P.3 (Theta)': 'other VOI'}

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Add attribute to color lineages for manuscript figures",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--input', type=str, metavar="JSON", required=True, help="input Auspice JSON")
    parser.add_argument('--output', type=str, metavar="JSON", required=True, help="output Auspice JSON")
    args = parser.parse_args()

    with open(args.input, "r") as f:
        auspice_json = json.load(f)

    # add manuscript_lineage to each node
    def attach_labels(n):

        if 'emerging_lineage' in n["node_attrs"].keys():
            n["node_attrs"]["manuscript_lineage"] = {}
            n["node_attrs"]["manuscript_lineage"]["value"] = manuscript_lineage_names[n["node_attrs"]["emerging_lineage"]["value"]]

        if "children" in n:
            for c in n["children"]:
                attach_labels(c)

    attach_labels(auspice_json["tree"])

    # add 'colorings' option for manuscript_lineage
    # just call it 'linage'
    auspice_json['meta']['colorings'] += [{"key": "manuscript_lineage", "scale": [["Alpha", "#5E1D9D"], ["Beta", "#416CCE"], ["Delta", "#89BB6B"], ["Gamma", "#E14F2A"], ["other VOI", "#DDAA3C"], ["basal", "#B5B7B8"]], "title": "Lineage", "type": "categorical"}]

    # make this the default coloring option
    auspice_json['meta']['display_defaults']['color_by'] = 'manuscript_lineage'


    with open(args.output, 'w') as f:
        json.dump(auspice_json, f, indent=2)
