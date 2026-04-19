#! /usr/bin/env python3

from anytree import AnyNode, RenderTree, PostOrderIter, PreOrderIter
from anytree.render import AsciiStyle
from anytree.search import findall
from xml.dom.minidom import parse
import sys
import re
import copy
import itertools
from collections import Counter

# Penalty formula constants
MIN_POSITION_SUPPORT = 10
GAP_PENALTY_MULTIPLIER = 15  # weight applied to sum_of_gaps/position_count in the penalty
DELTA_MODIFIER = 0.1
SUPPORT_DIV_WEIGHT = 0.2
DISTANCE_POSITION_DIVIDER = 500


def check_position_coverage(poly, data, all_parent_positions=[]):
    def return_uncovered():
        return 0, 0, 0, None

    # if a position upstream was mutated, update the branch support here (+1 for support)
    is_remutation = False

    # ignore insertions (maybe do that later...)
    if "." in poly:
        return return_uncovered()
    # haplogroup-info
    pos = re.search("[0-9]+", poly).group()
    base = re.search("[A-Zd]", poly).group()
    # pileup-info
    try:
        pile = data[pos]  # e.g. CCTC
        cov = len(pile)  # coverage
        if cov == 0:
            return return_uncovered()
        target = pile.count(base)  # on target
        perc = round((target / cov) * 100, 2)

        if poly.endswith("!"):  # a remutation
            # check, if the mutation is already represented in the branch leading here
            # e.g. 16311T in L3 --> 16311T! in several U subgroups
            if pos + base in all_parent_positions:
                return return_uncovered()
            elif any(x.startswith(pos) for x in all_parent_positions):
                # e.g. 15301A in L3'4'5'6 is unsupported, but 15301G! in N is supported
                is_remutation = True

    except KeyError:
        return return_uncovered()

    return perc, target, cov, is_remutation


def get_position_weights(xml_tree):
    positions = []
    for _xml_haplogroup in xml_tree.getElementsByTagName("haplogroup"):
        # Now parse all the positions (poly-tags) and update the dictionary
        for _xml_child in (
            _xml_haplogroup.childNodes
        ):  # child here means XML childs --> for parsing the POLYs
            if (
                _xml_child.nodeType == _xml_child.ELEMENT_NODE
                and _xml_child.tagName == "details"
            ):
                # Get data for each haplogroup-defining position
                # Get the 'poly' elements directly under the 'details' element
                for _xml_poly in _xml_child.getElementsByTagName("poly"):
                    # extract position from XML
                    # poly = e.g. 1234T
                    poly = _xml_poly.firstChild.data
                    positions.append(poly)
    counter = Counter(positions)
    weights = {pos: 1 / counter[pos] for pos in counter.keys()}
    return weights


#
# This is aweful spaghetti-code and I am sorry for this!!
#

# define input
xml_path = sys.argv[1]
pileup_path = sys.argv[2]
prefix = sys.argv[3]

# open XML file
with open(xml_path) as xml_file:
    xml_tree = parse(xml_file)

# Save uniqueness of positions
weights = get_position_weights(xml_tree)

# open pileup file
pileup_data = {}
with open(pileup_path) as pileup_file:
    for _line in pileup_file:
        _cols = _line.split("\t")
        _sequence = _cols[4]
        _quality = _cols[5]
        # this line ignores the masked bases at the end
        _good_bases = "".join([b for b, q in zip(_sequence, _quality) if q != "!"])
        # import position -> bases
        pileup_data[_cols[1]] = _good_bases.upper()

# create the anynode tree
# create dictionaries to be filled during parsing
# Nomenclature
# |
# x NODE: Position and reads
# |
# x NODE
# |
# V BRANCH (sum of nodes): Positions and reads

raw_data = {
    # Branch
    "branch_positions_covered": 0,
    "branch_positions_support": 0,
    "branch_reads_covered": 0,
    "branch_reads_support": 0,
    "unique_branch_positions_covered": 0,
    "unique_branch_positions_support": 0,
    "branch_positions": [],
    # Node
    "node_positions_covered": 0,
    "node_positions_support": 0,
    "node_reads_covered": 0,
    "node_reads_support": 0,
    "unique_node_positions_covered": 0,
    "unique_node_positions_support": 0,
    "node_positions": [],  # only positions from this haplogroup node
    "node_positions_rendered": [],  # node_positions, but including the read coverage statistics
    # More stats
    "gaps_required": 0,  # how many intermediate nodes were skipped to come here
    "sum_of_gaps": 0,
    "distance_to_root":0,
    "penalty": -1,  # for branches at the leaves, update later
    # Pre-rendered-values
    "pct_branch_position_support": 0.0,
    "pct_unique_branch_position_support": 0.0,
    "pct_branch_sequence_support": 0.0,
    "pct_node_position_support": 0.0,
    "pct_unique_node_position_support": 0.0,
    "pct_node_sequence_support": 0.0,
    "pct_branch_support_delta": 0.0,
}

# start of node-creation
node = AnyNode(id="RSRS", parent=None, data=raw_data.copy())
name_node_dict = {"RSRS": node}

# walk through the XML and fill the anynode-tree on the fly
for _xml_haplogroup in xml_tree.getElementsByTagName("haplogroup"):
    # thats the haplogroup label
    _haplogroup_name = _xml_haplogroup.getAttribute("name")

    # get the parent haplogroup name from the XML parent
    # to get the AnyTree node parent from the dict
    _xml_parent_node = _xml_haplogroup.parentNode

    if _xml_parent_node.tagName == "haplogroup":
        _xml_parent = _xml_parent_node.getAttribute("name")
    else:
        continue

    # Now get the actual parent-node (AnyTree) from the dict
    parent_node = name_node_dict[_xml_parent]

    # Add the data-dict for the current haplogroup node
    _data = copy.deepcopy(raw_data)
    # check if the parent has position-support:
    parent_support = parent_node.data["node_positions_support"] > 0

    # first update the dict
    _data.update(
        {
            "branch_positions_covered": parent_node.data[
                "branch_positions_covered"
            ],  # later: add node covered on top
            "branch_positions_support": parent_node.data[
                "branch_positions_support"
            ],  # later: add node support on top
            "unique_branch_positions_covered": parent_node.data[
                "unique_branch_positions_covered"
            ],  # later: add node covered on top
            "unique_branch_positions_support": parent_node.data[
                "unique_branch_positions_support"
            ],  # later: add node support on top
            "branch_positions": parent_node.data[
                "branch_positions"
            ].copy(),  # later: add node_positions on top
            "branch_reads_covered": parent_node.data[
                "branch_reads_covered"
            ],  # later: add node reads covered on top
            "branch_reads_support": parent_node.data[
                "branch_reads_support"
            ],  # later: add node reads support on top
            "gaps_required": 0
            if parent_support
            else parent_node.data["gaps_required"] + 1,
            "sum_of_gaps": parent_node.data["sum_of_gaps"],
            "distance_to_root":parent_node.data["distance_to_root"]+1
        }
    )

    # Now parse all the positions (poly-tags) and update the dictionary
    for _xml_child in (
        _xml_haplogroup.childNodes
    ):  # child here means XML childs --> for parsing the POLYs
        if (
            _xml_child.nodeType == _xml_child.ELEMENT_NODE
            and _xml_child.tagName == "details"
        ):
            # Get data for each haplogroup-defining position
            # Get the 'poly' elements directly under the 'details' element
            for _xml_poly in _xml_child.getElementsByTagName("poly"):
                # extract position from XML
                # poly = e.g. 1234T
                _poly = _xml_poly.firstChild.data
                _data["node_positions"].append(_poly)
                _data["branch_positions"].append(_poly)
                poly_weight = weights[_poly]

                # now parse the positions and calculate coverage
                perc, derived, cov, mutation = check_position_coverage(
                    _poly, pileup_data, parent_node.data["branch_positions"]
                )

                parsed_poly = (
                    f"{'**' if mutation else ''}{_poly} ({perc:.2f}% {derived}/{cov})"
                )

                _data["node_positions_rendered"].append(parsed_poly)

                _data["node_reads_covered"] += cov
                _data["node_reads_support"] += derived
                
                _data["branch_reads_covered"] += cov
                _data["branch_reads_support"] += derived
                
                # add the position to the count
                if cov > 0:
                    _data["node_positions_covered"] += 1
                    _data["branch_positions_covered"] += 1
                    if poly_weight == 1:
                        _data["unique_branch_positions_covered"] += 1
                        _data["unique_node_positions_covered"] += 1
                
                if perc > MIN_POSITION_SUPPORT:
                    _data["node_positions_support"] += 1
                    _data["branch_positions_support"] += 1
                    _data["branch_positions_support"] += mutation
                    if poly_weight == 1:
                        _data["unique_branch_positions_support"] += 1
                        _data["unique_node_positions_support"] += 1

    if _data["node_positions_support"] == 0:
        _data["sum_of_gaps"] += 1

    # add node to the tree
    tmp = AnyNode(id=_haplogroup_name, parent=parent_node, data=_data)
    name_node_dict.update({_haplogroup_name: tmp})


# Calculate Stats
for hap in PostOrderIter(node):
    hap.data["pct_node_sequence_support"] = (
        hap.data["node_reads_support"] / hap.data["node_reads_covered"] * 100
        if hap.data["node_reads_covered"] > 0
        else 0.0
    )
    hap.data["pct_branch_sequence_support"] = (
        hap.data["branch_reads_support"] / hap.data["branch_reads_covered"] * 100
        if hap.data["branch_reads_covered"] > 0
        else 0.0
    )
    hap.data["pct_branch_position_support"] = (
        hap.data["branch_positions_support"]
        / hap.data["branch_positions_covered"]
        * 100
        if hap.data["branch_positions_covered"] > 0
        else 0.0
    )
    hap.data["pct_unique_branch_position_support"] = (
        hap.data["unique_branch_positions_support"]
        / hap.data["unique_branch_positions_covered"]
        * 100
        if hap.data["unique_branch_positions_covered"] > 0
        else 0.0
    )
    hap.data["pct_node_position_support"] = (
        hap.data["node_positions_support"] / hap.data["node_positions_covered"] * 100
        if hap.data["node_positions_covered"] > 0
        else 0.0
    )
    hap.data["pct_unique_node_position_support"] = (
        hap.data["unique_node_positions_support"]
        / hap.data["unique_node_positions_covered"]
        * 100
        if hap.data["unique_node_positions_covered"] > 0
        else 0.0
    )
    hap.data["pct_branch_support_delta"] = abs(
        hap.data["pct_node_sequence_support"] - hap.data["pct_branch_sequence_support"]
    )

    ## Calculate Penalty
    ### Gap-penalty
    if hap.data['branch_positions_support'] > 0: # that should only not apply to the very root
        _gap_position_proportion = hap.data['sum_of_gaps'] / hap.data['branch_positions_support']
    else:
        _gap_position_proportion = 1
    
     ### Support-Delta
    
    hap.data["penalty"] = (
        _gap_position_proportion * GAP_PENALTY_MULTIPLIER
        + (100-hap.data['pct_branch_position_support']) * SUPPORT_DIV_WEIGHT
        + hap.data['pct_branch_support_delta'] * DELTA_MODIFIER
        - (hap.data['distance_to_root']*(hap.data['branch_positions_covered'] / DISTANCE_POSITION_DIVIDER ))
    )


#Test: update the penalty to penalize high penalty values upstream
for hap in PreOrderIter(node):
    if hap==node:
        continue
    hap.data['penalty'] += (hap.parent.data['penalty'] * 0.1)



def print_header(file):
    print(
        "\t".join(
            [
                "Order",
                "Parent",
                "Haplogroup",
                "PhyloTree",
                "Penalty",
                "DistanceToRoot",
                "SumOfGaps",
                "BranchPositionSupport#",
                "BranchPositionSupport%",
                "BranchPositionSupportUnique#",
                "BranchPositionSupportUnique%",
                "BranchSequenceSupport#",
                "BranchSequenceSupport%",
                "PositionSupport#",
                "PositionSupport%",
                "PositionSupportUnique#",
                "PositionSupportUnique%",
                "SequenceSupport#",
                "SequenceSupport%",
                "SupportDelta",
                "Positions",
            ]
        ),
        file=file,
    )


def print_line(row, file, n):
    d = row.node.data
    parent = row.node.parent.id if row.node.parent else "-"
    print(
        "\t".join(
            [
                str(n),
                parent,
                row.node.id,
                f"{row.pre.rstrip()} {row.node.id}",
                f"{d['penalty']:.2f}",
                f"{d['distance_to_root']}",
                f"{d['sum_of_gaps']}",
                f"{d['branch_positions_support']}/{d['branch_positions_covered']}",
                f"{d['pct_branch_position_support']:.2f}%",
                f"{d['unique_branch_positions_support']}/{d['unique_branch_positions_covered']}",
                f"{d['pct_unique_branch_position_support']:.2f}%",
                f"{d['branch_reads_support']}/{d['branch_reads_covered']}",
                f"{d['pct_branch_sequence_support']:.2f}%",
                f"{d['node_positions_support']}/{d['node_positions_covered']}",
                f"{d['pct_node_position_support']:.2f}%",
                f"{d['unique_node_positions_support']}/{d['unique_node_positions_covered']}",
                f"{d['pct_unique_node_position_support']:.2f}%",
                f"{d['node_reads_support']}/{d['node_reads_covered']}",
                f"{d['pct_node_sequence_support']:.2f}%",
                f"{d['pct_branch_support_delta']:.2f}%",
                f"{'; '.join(d['node_positions_rendered'])}",
            ]
        ),
        file=file,
    )


# Print summary stats for every haplogroup!
with open(f"{prefix}.raw.tsv", "w") as tree:
    print_header(tree)
    for n, row in enumerate(list(RenderTree(node, style=AsciiStyle)), 1):
        print_line(row, tree, n)
