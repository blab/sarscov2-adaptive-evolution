import argparse


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--tree", required=True, help="SARS-CoV-2 tree as an Auspice JSON")
    parser.add_argument("--gene", required=True, help="gene to test")
    parser.add_argument("--nonsyn-syn", required=True, help="mutation type to test")
    parser.add_argument("--iteration", required=True, help="iteration of randomization")
    parser.add_argument("--output", required=True, help="JSON with randomization results for a single iteration")

    args = parser.parse_args()

    # Do some randomizations!
    with open(args.output, "w") as output_handle:
        output_handle.write("Randomization!\n")
