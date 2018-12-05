from ete3 import Tree

class Node(object):
  def __init__(self, left, right):
    self.left = left
    self.right = right

  def __repr__(self):
    return '(%s, %s)' % (self.left, self.right)
class MyTree(Tree):

    def __init__(self, *args, **kwargs):
        super(MyTree, self).__init__(*args, **kwargs)

    def __eq__(self, other):
        rf = self.robinson_foulds(other, unrooted_trees=True)
        return not bool(rf[0])

trees = [MyTree(newick) for newick in all_newicks]
# Given a tree and a label, yields every possible augmentation of the tree by
# adding a new node with the label as a child "above" some existing Node or Leaf.
def add_leaf(tree, label):
  yield Node(label, tree)
  if isinstance(tree, Node):
    for left in add_leaf(tree.left, label):
      yield Node(left, tree.right)
    for right in add_leaf(tree.right, label):
      yield Node(tree.left, right)

# Given a list of labels, yield each rooted, unordered full binary tree with
# the specified labels.
def enum_unordered(labels):
  if len(labels) == 1:
    yield labels[0]
  else:
    for tree in enum_unordered(labels[1:]):
      for new_tree in add_leaf(tree, labels[0]):
        yield new_tree


from itertools import permutations

def ingroup_generator(species, n):
    for perm in permutations(species, n):
        yield tuple([tuple(perm), tuple(s for s in species if s not in perm)])

def format_newick(s, outgroup=''):
    return '(' + ', '.join('({})'.format(', '.join(p)) for p in s) + ',({}));'.format(outgroup)

def sort_newick(t):
    if isinstance(t, str):
        return sorted(t)
    else:
        if all(isinstance(c, str) for c in t):
            return sorted(t)
        if all(isinstance(l, list) for l in t):
            return [sort_newick(l) for l in sorted(t, key=lambda k: sorted(k))]
        else:
            return [sort_newick(l) for l in t]


def canonical_newick(n):
    n = n.replace(',', '')
    p = nestedExpr().parseString(n).asList()
    s = sort_newick(p)
    return str(s)

species = ["A","B","C","D","E"]
outgroup = "E"
ingroup = [s for s in species if s != outgroup]

itertools_newicks= []
for n in range(1, len(ingroup)):
    for p in ingroup_generator(ingroup, n):
        itertools_newicks.append(format_newick(p, outgroup))

enum_newicks= []
for t in enum_unordered(ingroup):
    enum_newicks.append('({},({}));'.format(t, outgroup))

from collections import defaultdict

all_newicks = itertools_newicks + enum_newicks

d = defaultdict(list)
for newick in all_newicks:
    d[canonical_newick(newick)].append(newick)

for canonical, newicks in d.items():
    print(canonical)
    for newick in newicks:
        print('\t', newick)
    print("")


from pyparsing import nestedExpr


