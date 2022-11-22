import faultdiagnosistoolbox as fdt
import sympy as sym
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import networkx as nx
from collections import Counter

modelDef={}
modelDef['type'] = 'Symbolic'
#modelDef['x']= ['p1','p2','p3', 'p4','q1','q2','q3','q4','dp1','dp2','dp3','dp4' ,'qin1' , 'qin2']
modelDef['x']= ['p1','p2','q1','dp1','qin1' , 'p3' ,'q2' , 'dp2','p4','q3','dp3','qin2', 'q4', 'dp4']
modelDef['f'] = ['f1','f2','f3','f4','f5','f6']
modelDef['z'] = ['u1','u2','y1' , 'y2' , 'y3' , 'y4' , 'y5' , 'y6']
modelDef['parameters'] = ['CT1' , 'CT2' , 'CT3' , 'CT4' , 'RP1' , 'RP2' , 'RP3' , 'RP4']

sym.var(modelDef['x']);
sym.var(modelDef['f']);
sym.var(modelDef['z']);
sym.var(modelDef['parameters']);
modelDef['rels'] = [
    -dp1+1/CT1*(qin1-q1-f1),
    -q1+(p1-p2)/(RP1+f2) ,
    fdt.DiffConstraint('dp1','p1'),
    -qin1+u1,
    -p1+y1,
    -q1+y2,
    -dp2+1/CT2*(q1-q2-f3),
    -q2+(p2-p3)/(RP2+f4),
    fdt.DiffConstraint('dp2','p2'),
    -p2+y3,
    -q2+y4,
    -dp3+1/CT3*(qin2+q2-q3),
    -q3+(p3-p4)/(RP3+f5),
    fdt.DiffConstraint('dp3','p3'),
    -qin2+u2,
    -q3+y5,
    -dp4+1/(CT4)*(q3-q4-f6),
    -q4+(p4/RP4),
    fdt.DiffConstraint('dp4','p4'),
    -p4+y6
    ]

four_tank = fdt.DiagnosisModel(modelDef, name='Four Tank System ')

# model informations
four_tank.Lint()


# Plot model
four_tank.PlotModel()
_ = plt.title('Figure 1 : FOUR TANK SYSTEM')
plt.show()

# STRUCTURAL ANALYSIS
equations = four_tank.e
unknown_variables = four_tank.x
known_variables = four_tank.z

# Incidence matrix
incidence_matrix_u_v = four_tank.X
incidence_matrix_k_v = four_tank.Z

# numpy matrix to pandas datzframe matrix
df_u_v = pd.DataFrame(incidence_matrix_u_v, index=equations, columns=unknown_variables)
df_k_v = pd.DataFrame(incidence_matrix_k_v, index=equations, columns=known_variables)

# Tripartite graph
B = nx.Graph()
B.add_nodes_from(df_k_v.columns, bipartite=0)
B.add_nodes_from(df_u_v.index, bipartite=1)
B.add_nodes_from(df_u_v.columns, bipartite=2)

s = df_u_v.stack()
B.add_edges_from(s[s >= 1].index)
ss = df_k_v.stack()
B.add_edges_from(ss[ss >= 1].index)

nodes = B.nodes()

# for each of the parts create a set
nodes_0 = set([n for n in nodes if B.nodes[n]['bipartite'] == 0])
nodes_1 = set([n for n in nodes if B.nodes[n]['bipartite'] == 1])
nodes_2 = set([n for n in nodes if B.nodes[n]['bipartite'] == 2])

# set the location of the nodes for each set
pos = dict()
pos.update((n, (1, i)) for i, n in enumerate(nodes_0))
pos.update((n, (2, i)) for i, n in enumerate(nodes_1))
pos.update((n, (3, i)) for i, n in enumerate(nodes_2))

# ploting result
nx.draw(B, pos=pos,
        node_color='lightgreen',
        edge_color='lightblue',
        with_labels=True)
_ = plt.title(' Figure 2 : Tripartite graph of the FOUR TANK SYSTEM ')
plt.show()


def is_valid_cycle(cycle):
    """ boolean function for testing a giving cycle
    and return true if is valid ant false otherwise
    @see : J.Tierman. An efficient search algorithm to find the elementary
    circuits of a graph. Comm. ACM vol 13, pages 722—726, 1970.

     Parameters
    ----------
    cycle : list of simple path

    Returns
    -------
    True if cycle is valid and None otherwise.

    """
    calculated_uknown_variables = []
    valid = True

    # convert a cycle to digraph/path_graph
    G = nx.path_graph(cycle, create_using=nx.DiGraph())

    # loop for each elem in the cycle or for each node in the path_graph G
    for elem in cycle:
        # succ_current_elem = list(G.successors(elem))[0]
        # if end of path then valid = true
        """
        if (succ_current_elem in known_variables):
            valid = True
            break
        else:
           """
        if (elem in unknown_variables):
            # each node in DiGraph has only one predecessor and one successor
            elem_equation = list(G.predecessors(elem))[0]
            # return all the unkown variable of the element equation from the
            # tripartite graph B
            # succ list may have known variable within
            succ = list(B[elem_equation])
            # intersection of succ list and unknown variable list
            u_v_of_elem_equation = list(set(succ) & set(unknown_variables))
            # if elem are the only unknown variable in the equation then
            # elem are added to calculated_uknown_variables
            if (len(u_v_of_elem_equation) == 1):
                calculated_uknown_variables.append(u_v_of_elem_equation[0])
            else:
                # intersection between the set of unknwon variables of elem equation and the set of
                # the calculated variable, if lenght var == 1 then the value of this variable can
                # be deducted from the other variables, break the loop otherwise
                var = list(set(u_v_of_elem_equation) & set(calculated_uknown_variables))
                if (len(var) == 1):
                    calculated_uknown_variables.append(var[0])
                else:
                    valid = None
                    u_c_var = elem
    if (valid):
        return valid, None
    else:
        return valid, u_c_var



def generateCycles():
    """ Return a list of valid cycles """
    valid_paths = []
    cal_paths = []
    unc_vars = []
    for v in nodes_0:
        for c in nodes_0:
            paths = nx.all_simple_paths(B, v, c)
            for path in paths:
                reversed_path = [path[-(i + 1)] for i in range(len(path))]
                valid, u_c_var = is_valid_cycle(path)
                if (valid) and (path not in valid_paths) and (reversed_path not in valid_paths):
                    valid_paths.append(path)
                elif (not valid):
                    unc_vars.append(u_c_var)
                    cal_path = path[:path.index(u_c_var) - 3]
                    #print(u_c_var)
                    if (cal_path not in cal_paths):
                        cal_paths.append(cal_path)

    return valid_paths, cal_paths


cycles, cp = generateCycles()
print("all cycles : ")
for count , c in enumerate(cycles):
    print(f"C{count + 1}: ", c)
print(len(cycles))

#degree of redundancy
def red_degree(paths):
    var = []
    for calPath in paths:
        var.append(calPath[-1])
    return Counter(var)
print('degree of redundancy')
rd = red_degree(cp)
#red_degree = Counter(uv)
print(rd)

print('cycles valid : ')
print(cycles)
print('Cycles sensivities :')

# sorting equation from each cycle
eq_sets = []
for c in cycles:
    eqs_labels = list(set(c) & set(equations))
    eqs = []
    for e in eqs_labels:
        # add equation number
        eqs.append(equations.index(e))
    eq_sets.append(eqs)

for k, mi in enumerate(eq_sets):
    print(f"equations in cycle {k + 1}: ", end='')
    print([four_tank.e[ei] for ei in np.sort(mi)], end=', sensitive to faults ')
    print([four_tank.f[fi] for fi in np.argwhere(four_tank.FSM([mi])[0]).reshape(-1)])

# residuals test selection
test_selection_arr = four_tank.TestSelection(eq_sets, 'aminc', four_tank.IsolabilityAnalysisArrs(eq_sets))
print(len(test_selection_arr))
print('result after test selection : ')
print(test_selection_arr)


def getResiduals(cycles):
    """Return the list of residuals for a set of cycles
     Parameters
     ----------
     cycles : list, of valid cycles
     """
    residuals = []
    for cycle in cycles:
        residuals.append(list(set(cycle) & set(known_variables)))
    return residuals


def SignatureMatrix(residuals):
    """Return the signature matrix for a set of residuals
    Parameters
    ----------
    residuals : list, of the residuals generated from cycles .
    """
    fsm_matrix = []
    i = 0

    for r in residuals:
        curr_row = []
        for v in known_variables:
            if r[0] == v:
                curr_row.append(1)
            elif r[1] == v:
                curr_row.append(1)
            else:
                curr_row.append(0)
            i = i + 1

        fsm_matrix.append(curr_row)

    return fsm_matrix


fsm_matrix = SignatureMatrix(getResiduals(cycles))

# ploting  signature matrix
#plt.subplot(5, 1, 3)
plt.spy(fsm_matrix, markersize=10, marker='o')
plt.xticks(np.arange(0, four_tank.nz()), known_variables)
plt.yticks(np.arange(0, len(fsm_matrix)),
           ["r " + str(k + 1) for k in np.arange(0, len(fsm_matrix))])
plt.gca().xaxis.tick_bottom()
_ = plt.title(' Figure 3 : Signature matrix ')
plt.show()

def FSM(eqs_sets, plot=False):
    """Return the fault signature matrix for a set of equation sets.

    Input
    -----
      eqs_sets : list of sets of equations, e.g., MSO sets
      Plot     : If True, plot the fault signature matrix (dafault: False)
    """
    r = np.zeros((len(eqs_sets), four_tank.F.shape[1]), dtype=np.int64)
    for idx, eqs in enumerate(eqs_sets):
        r[idx, :] = np.any(four_tank.F[eqs, :], axis=0)

    if plot:
        plt.spy(r, markersize=10, marker='o')
        plt.xticks(np.arange(0, four_tank.nf()), four_tank.f)
        plt.yticks(np.arange(0, len(eqs_sets)),
                   ["r " + str(k + 1) for k in np.arange(0, len(eqs_sets))])
        plt.gca().xaxis.tick_bottom()

    return r


# ploting fault signature matrix
#plt.subplot(5, 1, 4)
FSM(eq_sets, plot=True)
_ = plt.title(' Figure 3 : Fault Signature matrix ')
plt.show()
# print(known_variables)
# print(SignatureMatrix(getResiduals(generateCycles())))


# IsolabilityAnalysis
#plt.subplot(5, 1, 5)
four_tank.IsolabilityAnalysisFSM(four_tank.FSM(eq_sets), permute=False, plot=True)
_ = plt.title(' Figure 4 : Isolability matrix ')
plt.show()



