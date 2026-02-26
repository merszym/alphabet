#! /usr/bin/env python3

from anytree import AnyNode, RenderTree, PostOrderIter, PreOrderIter
from anytree.render import AsciiStyle
from anytree.search import findall, find
from xml.dom.minidom import parse
import sys
import re
import copy
import itertools
from collections import Counter

# Penalty formula constants
MIN_POSITION_SUPPORT_PERCENT = 5   # minimum allele frequency to count a position as supported
INNER_NODE_BONUS_THRESHOLD = 20    # branch positions needed before inner-node penalty drops to 0
GAP_PENALTY_MULTIPLIER = 3         # weight applied to sum_of_gaps in the penalty
DELTA_PENALTY_DIVISOR = 3          # divisor for sequence/branch support delta penalty
BRANCH_SUPPORT_PENALTY_DIVISOR = 5 # divisor for branch position support penalty
NODE_SUPPORT_PENALTY_DIVISOR = 3   # divisor for node position support penalty

def check_position_coverage(poly, data, all_parent_positions=[]):
    def return_uncovered():
        return 0,0,0, None

    #if a position upstream was mutated, update the branch support here (+1 for support)
    is_remutation = False

    # ignore insertions (maybe do that later...)
    if '.' in poly:
        return return_uncovered()
    #haplogroup-info
    pos = re.search('[0-9]+', poly).group()
    base = re.search('[A-Zd]', poly).group()
    # pileup-info
    try:
        pile = data[pos] # e.g. CCTC
        cov = len(pile) # coverage
        if cov == 0:
            return return_uncovered()
        target = pile.count(base) # on target
        perc = round((target/cov)*100, 2)

        if poly.endswith("!"): #a remutation
            # check, if the mutation is already represented in the branch leading here
            # e.g. 16311T in L3 --> 16311T! in several U subgroups
            if pos+base in all_parent_positions:
                return return_uncovered()
            elif any(x.startswith(pos) for x in all_parent_positions):
                # e.g. 15301A in L3'4'5'6 is unsupported, but 15301G! in N is supported
                is_remutation = True

    except KeyError:
        return return_uncovered()
    ##
    ## compare the poly with the mpileup_data to return coverage
    ##
    ##
    return perc, target, cov, is_remutation


def get_position_weights(xml_tree):
    positions = []
    for xml_haplogroup in xml_tree.getElementsByTagName('haplogroup'):
        # Now parse all the positions (poly-tags) and update the dictionary
        for xml_child in xml_haplogroup.childNodes: #child here means XML childs --> for parsing the POLYs
            if xml_child.nodeType == xml_child.ELEMENT_NODE and xml_child.tagName == 'details':
            # Get data for each haplogroup-defining position 
            # Get the 'poly' elements directly under the 'details' element
                for xml_poly in xml_child.getElementsByTagName('poly'):
                    # extract position from XML
                    # poly = e.g. 1234T
                    poly = xml_poly.firstChild.data
                    positions.append(poly)
    counter = Counter(positions)
    weights = {pos:1/counter[pos] for pos in counter.keys()}
    return weights
    

#
# This is aweful spaghetti-code and I am sorry for this!!
#

# define input
xml_path = sys.argv[1]
pileup_path = sys.argv[2]
prefix = sys.argv[3]
show_best = int(sys.argv[4])

# open XML file
with open(xml_path) as xml_file:
    xml_tree = parse(xml_file)

weights = get_position_weights(xml_tree)

# open pileup file
pileup_data = {}
with open(pileup_path) as pileup_file:
    for _line in pileup_file:
        _cols = _line.split('\t')
        _sequence = _cols[4]
        _quality = _cols[5]
        # this line ignores the masked bases at the end
        _good_bases = ''.join([b for b, q in zip(_sequence, _quality) if q != '!'])
        # import position -> bases
        pileup_data[_cols[1]] = _good_bases.upper()

# create the anynode tree
# create dictionaries to be filled during parsing
# |
# x NODE: Position and reads
# |
# x NODE
# |
# V BRANCH (sum of nodes): Positions and reads
raw_data = {
    'branch_reads_support':0,
    'branch_reads_covered':0,
    'node_reads_support':0,
    'node_reads_covered':0,
    'branch_positions_covered': 0,
    'branch_positions_support': 0,
    'unique_branch_positions_covered': 0,
    'unique_branch_positions_support': 0,
    'node_positions_covered': 0,
    'node_positions_support': 0,
    'unique_node_positions_covered': 0,
    'unique_node_positions_support': 0,
    'branch_positions': [],
    'node_positions': [], #only positions from this haplogroup node
    'node_positions_rendered': [], #node_positions, but including the read coverage statistics
    'gaps_required':0, # how many intermediate nodes were skipped to come here 
    'sum_of_gaps':0,
    'penalty':-1, #for branches at the leaves, update later
}

node = AnyNode(id='mtMRCA', parent=None, data=raw_data.copy())
name_node_dict = {'mtMRCA':node}
max_support = 0

# walk through the XML and fill the anynode-tree
for xml_haplogroup in xml_tree.getElementsByTagName('haplogroup'):
    # thats the haplogroup label
    name = xml_haplogroup.getAttribute('name')
    
    # get the parent haplogroup name from the XML parent
    # to get the AnyTree node parent from the dict
    xml_parent_node = xml_haplogroup.parentNode
    if xml_parent_node.tagName == 'haplogroup':
        parent = xml_parent_node.getAttribute('name')
    else:
        continue
    
    # Now get the actual parent-node (AnyTree) from the dict 
    parent_node = name_node_dict[parent]

    # Add the data-dict for the current haplogroup node
    data = copy.deepcopy(raw_data)
    # check if the parent has position-support:
    parent_support = parent_node.data['node_positions_support'] > 0
    
    data.update({
        'branch_positions_covered': parent_node.data['branch_positions_covered'], # later: add node covered on top 
        'branch_positions_support': parent_node.data['branch_positions_support'], # later: add node support on top
        'unique_branch_positions_covered': parent_node.data['unique_branch_positions_covered'], # later: add node covered on top 
        'unique_branch_positions_support': parent_node.data['unique_branch_positions_support'], # later: add node support on top
        'branch_positions': parent_node.data['branch_positions'].copy(), # later: add node_positions on top
        'branch_reads_covered': parent_node.data['branch_reads_covered'], # later: add node reads covered on top
        'branch_reads_support': parent_node.data['branch_reads_support'], # later: add node reads support on top
        'gaps_required':0 if parent_support else parent_node.data['gaps_required'] + 1,
        'sum_of_gaps': parent_node.data['sum_of_gaps']
    })

    # Now parse all the positions (poly-tags) and update the dictionary
    for xml_child in xml_haplogroup.childNodes: #child here means XML childs --> for parsing the POLYs
        if xml_child.nodeType == xml_child.ELEMENT_NODE and xml_child.tagName == 'details':
        # Get data for each haplogroup-defining position 
        # Get the 'poly' elements directly under the 'details' element
            for xml_poly in xml_child.getElementsByTagName('poly'):
                # extract position from XML
                # poly = e.g. 1234T
                poly = xml_poly.firstChild.data
                data['node_positions'].append(poly)
                data['branch_positions'].append(poly)
                poly_weight = weights[poly]

                # now parse the positions and calculate coverage
                perc, target, cov, mutation = check_position_coverage(poly, pileup_data, parent_node.data['branch_positions'])

                parsed_poly = f"{'**' if mutation else ''}{poly} ({perc:.2f}% {target}/{cov})"
                data['node_positions_rendered'].append(parsed_poly)

                data['node_reads_support'] += target
                data['node_reads_covered'] += cov
                data['branch_reads_support'] += target
                data['branch_reads_covered'] += cov

                if cov > 0:
                    data['node_positions_covered'] += 1
                    data['branch_positions_covered'] += 1
                    if poly_weight == 1:
                        data['unique_branch_positions_covered'] += 1
                        data['unique_node_positions_covered'] += 1
                if perc > MIN_POSITION_SUPPORT_PERCENT:
                    data['node_positions_support'] += 1
                    data['branch_positions_support'] += 1
                    data['branch_positions_support'] += mutation
                    if poly_weight == 1:
                        data['unique_branch_positions_support'] += 1
                        data['unique_node_positions_support'] += 1
    
    if data['node_positions_support'] == 0:
        data['sum_of_gaps'] += 1

    if data['branch_positions_support'] > max_support:
        max_support = data['branch_positions_support']


    #add node to the tree
    tmp = AnyNode(id=name, parent=parent_node, data=data)
    name_node_dict.update({name:tmp})   

# now summarize stats on each node
min_penalty = []

for hap in PostOrderIter(node):
    # The penalty is used to get the best supported node. The smaller, the better

    branch_sequence_support = 0
    sequence_support = 0
    delta_penalty = 0
    branch_support_penalty = 5
    node_support_penalty = 0

    if hap.data['node_reads_covered'] > 0:
        sequence_support = hap.data['node_reads_support']/hap.data['node_reads_covered'] * 100
        
    if hap.data['branch_reads_covered'] > 0:
        branch_sequence_support = hap.data['branch_reads_support']/hap.data['branch_reads_covered'] * 100

        delta = max(sequence_support, branch_sequence_support) - min(sequence_support, branch_sequence_support)
        delta_penalty = delta // DELTA_PENALTY_DIVISOR

    if hap.data['unique_branch_positions_covered'] > 0:
        branch_support = hap.data['unique_branch_positions_support'] / hap.data['unique_branch_positions_covered'] * 100
        branch_support_penalty = (100 - branch_support) // BRANCH_SUPPORT_PENALTY_DIVISOR
    
    if hap.data['node_positions_covered'] > 0:
        node_support = hap.data['node_positions_support'] / hap.data['node_positions_covered'] * 100
        node_support_penalty = (100 - node_support) // NODE_SUPPORT_PENALTY_DIVISOR

    inner_node_penalty = max(0, INNER_NODE_BONUS_THRESHOLD - hap.data['branch_positions_support'])

    hap.data['penalty'] = (
            hap.data['sum_of_gaps'] * GAP_PENALTY_MULTIPLIER +
            int(delta_penalty) +
            int(branch_support_penalty) +
            int(node_support_penalty) +
            inner_node_penalty -
            int(hap.data['unique_branch_positions_support'] > 0)
        )
    
    min_penalty.append(hap.data['penalty'])

def print_header(file):
    print('\t'.join(
        [
            'Order',
            'Parent',
            'PhyloTree',
            'Penalty',
            'RequiredGaps',
            'SumOfGaps',
            'BranchPositionSupport#',
            'BranchPositionSupport%',
            'BranchPositionSupportUnique#',
            'BranchPositionSupportUnique%',
            'BranchSequenceSupport#',
            'BranchSequenceSupport%',
            'PositionSupport#',
            'PositionSupport%',
            'PositionSupportUnique#',
            'PositionSupportUnique%',
            'SequenceSupport#',
            'SequenceSupport%',
            'SupportDelta',
            'Positions'
        ]
    ), file=file)

def print_line(row, file, n):
    parent = '-'
    branch_position_support = 0
    branch_sequence_support = 0
    unique_branch_position_support = 0
    position_support = 0
    unique_position_support = 0
    sequence_support = 0
    branch_support_delta = 0

    if row.node.data['node_reads_covered'] > 0:
        sequence_support = row.node.data['node_reads_support']/row.node.data['node_reads_covered'] * 100
    
    if row.node.data['branch_positions_covered'] > 0:
        branch_position_support = row.node.data['branch_positions_support']/row.node.data['branch_positions_covered'] * 100
    
    if row.node.data['unique_branch_positions_covered'] > 0:
        unique_branch_position_support = row.node.data['unique_branch_positions_support']/row.node.data['unique_branch_positions_covered'] * 100

    if row.node.data['branch_reads_covered'] > 0:
        branch_sequence_support = row.node.data['branch_reads_support']/row.node.data['branch_reads_covered'] * 100
        try:
            branch_support_delta = max(sequence_support, branch_sequence_support) - min(sequence_support, branch_sequence_support)
        except Exception:
            branch_support_delta = 0

    if row.node.data['node_positions_covered'] > 0:
        position_support = row.node.data['node_positions_support']/row.node.data['node_positions_covered'] * 100
    
    if row.node.data['unique_node_positions_covered'] > 0:
        unique_position_support = row.node.data['unique_node_positions_support']/row.node.data['unique_node_positions_covered'] * 100

    if row.node.parent:
        parent = row.node.parent.id
    
    print('\t'.join(
            [
                str(n),
                parent,
                f"{row.pre.rstrip()} {row.node.id}",
                f"{row.node.data['penalty']:.1f}",
                f"{row.node.data['gaps_required']}",
                f"{row.node.data['sum_of_gaps']}",
                f"{row.node.data['branch_positions_support']}/{row.node.data['branch_positions_covered']}",
                f"{branch_position_support:.2f}%",
                f"{row.node.data['unique_branch_positions_support']}/{row.node.data['unique_branch_positions_covered']}",
                f"{unique_branch_position_support:.2f}%",
                f"{row.node.data['branch_reads_support']}/{row.node.data['branch_reads_covered']}",
                f"{branch_sequence_support:.2f}%",
                f"{row.node.data['node_positions_support']}/{row.node.data['node_positions_covered']}",
                f"{position_support:.2f}%",
                f"{row.node.data['unique_node_positions_support']}/{row.node.data['unique_node_positions_covered']}",
                f"{unique_position_support:.2f}%",
                f"{row.node.data['node_reads_support']}/{row.node.data['node_reads_covered']}",
                f"{sequence_support:.2f}%",
                f"{branch_support_delta:.2f}%",
                f"{'; '.join(row.node.data['node_positions_rendered'])}"
            ]
        ), file=file
    )

# and print summary files
rendered = list(RenderTree(node, style=AsciiStyle))

# First, print stats for every haplogroup!
with open(f"{prefix}.raw.tsv", 'w') as tree1:
    print_header(tree1)
    for n,row in enumerate(rendered, 1):
        print_line(row, tree1, n)

# print the best tree (output all nodes with the two lowest penalties)
sorted_penalties = sorted(list(set(min_penalty)))

# find the nodes that have the lowest penalty (lowest with the provided wiggle-room)
best_nodes = findall(node, filter_ = lambda node: any([node.data['penalty'] in sorted_penalties[:show_best]]))
keep = []
for _best_node in best_nodes:
    keep.extend([x.id for x in _best_node.path])

# now copy the full tree and filter it for rendering
best_tree = copy.deepcopy(node)
for hap in PostOrderIter(best_tree):
    # eliminate all nodes that are not in the best_nodes
    if hap.id in keep:
        continue
    else:
        hap.parent = None

tree2_render = list(RenderTree(best_tree, style=AsciiStyle))
with open(f"{prefix}.best.tsv", 'w') as tree2:    
    print_header(tree2)
    for n,row in enumerate(tree2_render, 1):
        print_line(row, tree2, n)


## Print the best node (min_penalty) stats to StdOut
best_nodes = findall(best_tree, filter_ = lambda node: node.data['penalty'] == min(min_penalty))
note = ""

# Option 1: Only one best
if len(best_nodes) == 1:
    best_node = best_nodes[0]
    note = "Lowest Penalty"

# Option 2: Multiple: Check if all are on the same branch ---
# A set of nodes is on one branch if one node is ancestor of all others
else: 
    def is_ancestor(a, b):
        return a in b.path

    # if nothing was found:
    if len(best_nodes)==0:
        note = "No Haplogroup detected"
        best_node = best_tree

    else:

        # Sort by depth (shortest path = most upstream)
        matches_sorted = sorted(best_nodes, key=lambda n: len(n.path))

        #most-upstream
        best_node = matches_sorted[0]

        #check if the best_node is ancestral to all other best nodes
        if all(is_ancestor(best_node, n) for n in matches_sorted):
            note = "Highest Haplogroup with Lowest Penalty"
            pass
    
        else:
            # Otherwise: find Lowest Common Ancestor (LCA)
            note = "Shared Ancestor of Haplogroups with Lowest Penalty"
            paths = [node.path for node in matches_sorted]

            # Walk level by level until paths diverge
            for nodes_at_level in itertools.zip_longest(*paths):
                #check on each level, if the paths are all the same
                if len(set(nodes_at_level))==1:
                    best_node = nodes_at_level[0]
                else:
                    break

# Now print the stats to sdtout
if best_node.data['branch_positions_covered'] > 0:
    _branch_support = best_node.data['branch_positions_support']/best_node.data['branch_positions_covered'] * 100
else:
    _branch_support = 0

print(
    "Phylotree", 
    "BranchSupport", 
    "Penalty", 
    "SumOfGaps", 
    "SequenceSupport",
    "Note", 
    sep='\t', 
    file=sys.stdout,
    end='\n'
)

print(
    best_node.id,
    f"{best_node.data['branch_positions_support']}/{best_node.data['branch_positions_covered']} ({_branch_support:.2f}%)",
    f"{best_node.data['penalty']}",
    f"{best_node.data['sum_of_gaps']}",
    f"{best_node.data['branch_reads_support']}/{best_node.data['branch_reads_covered']}",
    note,
    sep="\t",
    file=sys.stdout,
    end="\n"
)
