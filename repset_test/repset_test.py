import subprocess
import pickle
from itertools import combinations

import numpy as np
import gurobipy as gp
from gurobipy import GRB, quicksum


def saveToPickleFile(python_object, path_to_file='object.pkl'):
    """
    Save python object to pickle file
    """
    out_file = open(path_to_file,'wb')
    pickle.dump(python_object, out_file)
    out_file.close()


def terminalExecute(command_str: str,
                    suppress_shell_output=False,
                    work_dir: str = None,
                    return_output=False) -> subprocess.STDOUT:
    """
    Execute given command in terminal through Python
    """
    if suppress_shell_output:
        stdout = subprocess.DEVNULL
    else:
        stdout = None
    output = subprocess.run(
        command_str, shell=True,
        cwd=work_dir, capture_output=return_output,
        stdout=stdout
        )
    return output
    

def computeRepsetOptimum(repset_out: str, esl_out: str) -> float:
    """
    Compute Repset optimal value based on PI data from esl-alipid output
    """
    with open(repset_out, "r") as repsetout:

        node_order = repsetout.read().split('\n')
        node_order.remove('')
        order_edges = [
            tuple(sorted(c, key=lambda x: int(x[1:].strip())))
            for c in combinations(node_order, 2)
            ]

    with open(esl_out, "r") as eslout:
        pi_data = {}
        lines = eslout.readlines()
        for line in lines[1:]:
            edge = tuple(line.split()[:2])
            pi = line.split()[2]
            pi_data[edge] = pi

    return sum([float(pi_data[edge]) for edge in order_edges])


def solveMILP(n_nodes, n_edges, edge_pairs, w, size):

    # Create a new model
    m = gp.Model("graph")

    # Create variables
    x = m.addVars(n_nodes, vtype=GRB.BINARY, name="x")
    y = m.addVars(n_edges, vtype=GRB.BINARY, name="y")

    for n, pair in enumerate(edge_pairs):
        i, j = pair
        y[n].setAttr(GRB.Attr.VarName, f"y{i}{j}")

    m.update()

    # Add edge constraints
    for pair in edge_pairs:
        i, j = pair
        y_ij = m.getVarByName(f"y{i}{j}")
        m.addConstr(x[i] + x[j] - 2 * y_ij >= 0, f"c{i}{j}>=0")
        m.addConstr(x[i] + x[j] - 2 * y_ij <= 1, f"c{i}{j}<=1")

    # Add minimum subgraph size constraint
    m.addConstr(quicksum(x[i] for i in range(n_nodes)) == size, "subgraph_size")
    
    # Set objective
    m.setObjective(quicksum(w[i]*y[i] for i in range(n_edges)), GRB.MINIMIZE)

    # Optimize model
    m.optimize()
    min_val = m.objval

    # Set objective
    m.setObjective(quicksum(w[i]*y[i] for i in range(n_edges)), GRB.MAXIMIZE)

    # Optimize model
    m.optimize()
    max_val = m.objval

    return (min_val, max_val)


# RUN TEST
# *********************************
n_samples = 50
repset_dist_to_milp = []  # position of repset's solution between milps min and max (relative to min)

outdir = "/home/robaina/Documents/TRAITS/repset_test/repset_results"
pifile = "/home/robaina/Documents/TRAITS/repset_test/esl_output.txt"

for n in range(n_samples):
    
    print(f"Running sample {n} out of {n_samples}")

    n_nodes = 50
    n_edges = int(0.5 * (n_nodes * (n_nodes - 1)))
    edge_pairs = list(combinations(range(n_nodes), 2))  # 01, 02, 03, 04, 12, 13, 14, 23, 24, 34
    w = 100 * np.random.rand(n_edges)
    size = 10

    # write els-alipd like output:
    with open(pifile, "w") as outfile:
        first_line = "# seqname1 seqname2 %id nid denomid %match nmatch denommatch\n"
        outfile.write(first_line)

        for n, (i, j) in enumerate(edge_pairs):
            edge = f"n{i} n{j}"
            line = f"n{i} n{j}  {w[n]}     39    100 100.00    100    100\n"
            outfile.write(line)
    
    # Get MILP optimal value
    milp_val = solveMILP(n_nodes, n_edges, edge_pairs, w, size)
    
    # Get repset optimal value
    cmd_str = f"python3 /home/robaina/Software/repset/repset.py --outdir {outdir} --pi {pifile} --size {size} --mixture 1"
    terminalExecute(cmd_str, suppress_shell_output=False)

    repset_val = computeRepsetOptimum(
        repset_out="/home/robaina/Documents/TRAITS/repset_test/repset_results/repset.txt",
        esl_out=pifile
        )
     
    milp_repset_ratio = (milp_val[1] - repset_val) / (milp_val[1] - milp_val[0])
    repset_dist_to_milp.append(milp_repset_ratio)


saveToPickleFile(repset_dist_to_milp, 'repset_dist_to_milp.pickle')
print('Finished!')