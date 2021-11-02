from itertools import combinations 

import numpy as np
import gurobipy as gp
from gurobipy import GRB, quicksum


def solveMILP(n_nodes, n_edges, edge_pairs, w, size):

    # Create a new model
    m = gp.Model("graph")

    # Create variables
    x = m.addVars(n_nodes, vtype=GRB.BINARY, name="x")
    y = m.addVars(n_edges, vtype=GRB.BINARY, name="y")

    for n, pair in enumerate(edge_pairs):
        i, j = pair
        y[n].setAttr(GRB.Attr.VarName, f"y{i}{j}")

    # Set objective
    m.setObjective(quicksum(w[i]*y[i] for i in range(n_edges)), GRB.MINIMIZE)
    m.update()

    # Add edge constraints
    for pair in edge_pairs:
        i, j = pair
        y_ij = m.getVarByName(f"y{i}{j}")
        m.addConstr(x[i] + x[j] - 2 * y_ij >= 0, f"c{i}{j}>=0")
        m.addConstr(x[i] + x[j] - 2 * y_ij <= 1, f"c{i}{j}<=1")

    # Add minimum subgraph size constraint
    m.addConstr(quicksum(x[i] for i in range(n_nodes)) == size, "subgraph_size")

    # Optimize model
    m.optimize()
    print('Obj: %g' % m.objVal)



n_nodes = 50
n_edges = int(0.5 * (n_nodes * (n_nodes - 1)))
edge_pairs = list(combinations(range(n_nodes), 2))  # 01, 02, 03, 04, 12, 13, 14, 23, 24, 34
w = 100 * np.random.rand(n_edges)
size = 10

# write els-alipd like output:
outfile = "/home/robaina/Documents/TRAITS/repset_test/esl_output.txt"

with open(outfile, "w") as outfile:
    first_line = "# seqname1 seqname2 %id nid denomid %match nmatch denommatch\n"
    outfile.write(first_line)

    for n, (i, j) in enumerate(edge_pairs):
        edge = f"n{i} n{j}"
        line = f"n{i} n{j}  {w[n]}     39    100 100.00    100    100\n"
        outfile.write(line)




solveMILP(n_nodes, n_edges, edge_pairs, w, size)