#! /usr/bin/env python3

import pandas as pd
from anytree import AnyNode, RenderTree, PostOrderIter
from anytree.render import AsciiStyle
from anytree.search import find
import sys
from pathlib import Path


def get_tree(df):
    """
    quickly reconstruct the tree from the dataframe
    """
    node = AnyNode(id="RSRS", parent=None)

    name_node_dict = {"RSRS": node}

    for _, _row in df.iterrows():
        # thats the haplogroup label
        _name = _row["Haplogroup"].strip()
        if _name == "RSRS":
            continue

        _parent = _row["Parent"].strip()

        # Now get the actual parent-node (AnyTree) from the dict
        _parent_node = name_node_dict.get(_parent, node)

        # add node to the tree
        _tmp = AnyNode(
            id=_name,
            parent=_parent_node,
        )
        name_node_dict.update({_name: _tmp})
    return node


def update_table(df):
    noise_filter = df[
        (
            df["BranchPositionSupportUnique%"].apply(
                lambda x: float(x.split("%")[0]) > 50
            )
        )
        & (df["SequenceSupport%"].apply(lambda x: float(x.split("%")[0]) > 10))
    ].copy()

    return noise_filter


def extract_best(df):
    _quantile5 = df["Penalty"].quantile(0.05)
    _filter = df[df["Penalty"] < _quantile5].copy()

    # now get the penalty values
    _p = _filter["Penalty"]

    _filter.insert(
        _filter.columns.get_loc("Penalty"),
        "Support",
        (_p.max() - _p) / (_p.max() - _p.min()),
    )

    _filter["Support"] = _filter["Support"].fillna(
        1
    )  # if max and min are the same (one line remaining)

    _results = _filter[_filter.Support > 0.3]

    return _results


def get_full_path(node, haplogroups):
    full_path_nodes = []

    for _hap in haplogroups:
        _node = find(node, lambda x: x.id == _hap)
        if not _node:
            continue
        _path = [x.id for x in _node.iter_path_reverse()][::-1]
        full_path_nodes.extend([x for x in _path if x not in full_path_nodes])

    return full_path_nodes


raw_tsv = Path(sys.argv[1])
output_tsv = raw_tsv.name.replace(".raw.tsv", ".best.tsv")

#
# 1. Open the Raw TSV file and extract the best supported Haplogroup Nodes
#

df = pd.read_csv(raw_tsv, sep="\t", keep_default_na=False)

filtered = update_table(df)

best = extract_best(filtered)

#
# 2. Subset the raw TSV to contain the full path (update the PhyloTree relationship)
#

node = get_tree(df)

valid_set = get_full_path(node, set(best["Haplogroup"]))

# i) Prune the tree
for _node in PostOrderIter(node):
    if _node.id not in valid_set:
        _node.parent = None

# ii) Sort tree by their order in df
haplogroup_order = {name.strip(): i for i, name in enumerate(df["Haplogroup"])}
for _node in PostOrderIter(node):
    if _node.children:
        _node.children = sorted(
            _node.children, key=lambda c: haplogroup_order.get(c.id, -1)
        )

# iii) Filter df to contain only valid nodes (preserving df row order)
subset_df = df[df["Haplogroup"].str.strip().isin(valid_set)].copy()

# iv) Update the PhyloTree column with the re-rendered ASCII tree
phylotree_map = {
    row.node.id: f"{row.pre.rstrip()} {row.node.id}"
    for row in RenderTree(node, style=AsciiStyle)
}
subset_df["PhyloTree"] = subset_df["Haplogroup"].str.strip().map(phylotree_map)

subset_df.to_csv(output_tsv, sep="\t", index=False)

#
# 4. Output a summary for the final report
#
