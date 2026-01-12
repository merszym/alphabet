#! /usr/bin/env python3

from anytree import AnyNode, RenderTree, PostOrderIter, LevelOrderGroupIter
from anytree.render import AsciiStyle
from xml.dom.minidom import parse
import sys
import re
import copy

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
raw_data = {
    'branch_reads_support':0,
    'branch_reads_covered':0,
    'node_reads_support':0,
    'node_reads_covered':0,
    'branch_positions_covered': 0,
    'branch_positions_support': 0,
    'node_positions_covered': 0,
    'node_positions_support': 0,
    'branch_positions': [],
    'node_positions': [], #only positions from this haplogroup node
    'node_positions_rendered': [], #node_positions, but including the read coverage statistics
    'penalty':-1 #for branches at the leaves, update later
}

node = AnyNode(id='mtMRCA', parent=None, data=raw_data.copy())
name_node_dict = {'mtMRCA':node}

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
    data.update({
        'branch_positions_covered': parent_node.data['branch_positions_covered'], # later: add node covered on top 
        'branch_positions_support': parent_node.data['branch_positions_support'], # later: add node support on top
        'branch_positions': parent_node.data['branch_positions'].copy(), # later: add node_positions on top
        'branch_reads_covered': parent_node.data['branch_reads_covered'], # later: add node reads covered on top
        'branch_reads_support': parent_node.data['branch_reads_support'], # later: add node reads support on top
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
                if perc > 10:
                    data['node_positions_support'] += 1
                    data['branch_positions_support'] += 1
                    data['branch_positions_support'] += mutation

    #add node to the tree
    tmp = AnyNode(id=name, parent=parent_node, data=data)
    name_node_dict.update({name:tmp})   

# now summarize stats on each node
for hap in PostOrderIter(node):
    #walk from the leaves up and set the penalty 
    if hap.data['node_positions_support'] != 0:
        hap.data['penalty'] = 0
        continue
    if hap.id == 'mtMRCA': # root
        hap.data['penalty'] = 0
        continue

    # introduce the penalty for nodes in a branch so that
    # covered, but dont have the target allel
    #   Tree                       Branch   Position  Penalty
    #  +-- L1c3b'c                 0/1      0/1       3 (Distance to covered tip
    #      |-- L1c3b               0/3      0/2       2 (Distance to covered tip)
    #      |   |-- L1c3b1          0/6      0/3       1 (Distance to covered tip)
    #      |   |   |-- L1c3b1a     1/8      1/2       0 (Tip)
    #      |   |   +-- L1c3b1b     0/8      0/2       -1 (Tip)
    #      |   +-- L1c3b2          0/10     0/7       -1 (Tip)
    #      +-- L1c3c               0/10     0/9       -1 (Tip)
    #
    # For each node, get the minimum Penalty from the children
    # keep -1 if -1 else add plus 1
    try:
        # check if a children exists that doesnt have -1 as penalty
        child_penalty = min(x.data['penalty'] for x in hap.children if x.data['penalty'] != -1)
    except:
        # keep the penalty at -1
        # set to -2, so that adding 1 makes it -1 again :D
        child_penalty = -2
        
    hap.data['penalty'] = child_penalty + 1

def print_header(file):
    print('\t'.join(
        [
            'Order',
            'Parent',
            'PhyloTree',
            'Penalty',
            'BranchSupport',
            'BranchSupportPercent',
            #'BranchSupportReads',
            #'NodeReadSupport',
            'PositionSupport',
            'SequenceSupport'
        ]
    ), file=file)

def print_line(row, file, n):
    if row.node.data['branch_positions_covered'] > 0:
        branch_support = row.node.data['branch_positions_support']/row.node.data['branch_positions_covered'] * 100
    else:
        branch_support = 0
    if row.node.parent:
        parent = row.node.parent.id
    else:
        parent = '-'
    print('\t'.join(
            [
                str(n),
                parent,
                f"{row.pre.rstrip()} {row.node.id}",
                f"{row.node.data['penalty']}",
                f"{row.node.data['branch_positions_support']}/{row.node.data['branch_positions_covered']}",
                f"{branch_support:.2f}%",
                # f"{row.node.data['branch_reads_support']}/{row.node.data['branch_reads_covered']}",
                # f"{row.node.data['node_reads_support']}/{row.node.data['node_reads_covered']}",
                f"{row.node.data['node_positions_support']}/{row.node.data['node_positions_covered']}",
                f"{'; '.join(row.node.data['node_positions_rendered'])}"
            ]
        ), file=file
    )

# and print summary files
rendered = list(RenderTree(node, style=AsciiStyle))

# First, print stats for every haplogroup!
with open(f"{prefix}.tree_all_groups.tsv", 'w') as tree1:
    print_header(tree1)
    for n,row in enumerate(rendered, 1):
        print_line(row, tree1, n)

#Then, prune the tree, so that:
# - nodes with less than 70% branch support are removed
# - nodes with 0 coverage are removed
# - nodes with penalty 4 or more are removed (actually, this is now mostly driven by the branch support...)
# 
# iterate through the tree and set all parents to None if these conditions apply

for hap in PostOrderIter(node):
    #walk from the leaves up and count the number of successful (leave) haplogroups on each node:
    #thats what the PostOrderIter does
    if hap.id == 'mtMRCA':
        continue # root
    
    if hap.data['branch_positions_covered'] > 0:
        branch_support = hap.data['branch_positions_support']/hap.data['branch_positions_covered'] * 100
    else:
        continue
        #
        # if branch_positions_covered == 0, then we know that we are at the very top, but maybe not yet at the root.
        # i.e. the very early positions (e.g. on L1'2'3'4'5'6) lack support. 
        # We dont want the full tree to be cropped, so skip them! We basically start the tree at a lower node in this case :)
        # 


    if branch_support < 70 or hap.data['penalty'] >= 4 or hap.data['penalty'] == -1:
        hap.parent = None
    
    # if a node was filtered because of too little support, also remove unsupported intermediate nodes
    if hap.data['node_positions_support'] == 0 and len(hap.children)==0:
        hap.parent = None

prune_rendered = list(RenderTree(node, style=AsciiStyle))
with open(f"{prefix}.tree_70perc_support.tsv", 'w') as tree2:    
    print_header(tree2)
    for n,row in enumerate(prune_rendered, 1):
        print_line(row, tree2, n)

#
#
# Do the same thing again, but now for 80% support and 3 gaps only
#
#


#Then, prune the tree, so that:
# - nodes with less than 70% branch support are removed
# - nodes with 0 coverage are removed
# - nodes with penalty 4 or more are removed (actually, this is now mostly driven by the branch support...)
# 
# iterate through the tree and set all parents to None if these conditions apply

for hap in PostOrderIter(node):
    #walk from the leaves up and count the number of successful (leave) haplogroups on each node:
    #thats what the PostOrderIter does
    if hap.id == 'mtMRCA':
        continue # root
    
    if hap.data['branch_positions_covered'] > 0:
        branch_support = hap.data['branch_positions_support']/hap.data['branch_positions_covered'] * 100
    else:
        continue

    if branch_support < 70 or hap.data['penalty'] >= 3:
        hap.parent = None
    
        # if a node was filtered because of too little support, also remove unsupported intermediate nodes
    if hap.data['node_positions_support'] == 0 and len(hap.children)==0:
        hap.parent = None

tree3_render = list(RenderTree(node, style=AsciiStyle))
with open(f"{prefix}.tree_70perc_support_3gaps.tsv", 'w') as tree3:    
    print_header(tree3)
    for n,row in enumerate(tree3_render, 1):
        print_line(row, tree3, n)