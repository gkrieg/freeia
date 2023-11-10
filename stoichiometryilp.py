
#from ilp import *
import cplex
from xml.dom import minidom
import time

'''
This file contains the ILP for finding a flux-maintaining set of hyperedge that create a sink from a source of minimum size.

'''

#first write the ILP, then solve it

def write_maximize_flux_objective(ilpfile):
    ilpfile.write('Maximize\n')
    ilpfile.write('1.0 psilon\n')

def write_stoichiometry_objective(ilpfile,H,source_edges,binary_vars=True,sourcedict=None):
    #Objective function is to minimize the number of edge variables used
    ilpfile.write('Minimize\n')
    for e in H.hyperedge_id_iterator():
        if binary_vars == True:
            ilpfile.write('+ {} '.format(b(e)))
        else:
            ilpfile.write('+ {} '.format(f(e)))
    if sourcedict != None:
        sourcecoef = 1.0/(len(sourcedict) + 1)
        for so in sourcedict:
            ilpfile.write('+ {} {} '.format(sourcecoef,sourcedict[so]))

    ilpfile.write('\n')

def write_stoichiometry_ilp(ilpfile,H,stoichiometries,targets,source_edges,positiveflux=False,sources=[],binary_vars=True,negative_regulators=None,reversible=None,sourcevars=False,sourcenum=None,minimizesources=False,nonparametric=False):
    if sourcevars == True and nonparametric == False:
        sourcedict = write_minimal_precursor_objective(ilpfile,sources,justgetsourcedict=True)
        if minimizesources == True:
            write_stoichiometry_objective(ilpfile,H,source_edges,binary_vars=binary_vars,sourcedict=sourcedict)
        else:
            write_stoichiometry_objective(ilpfile,H,source_edges,binary_vars=binary_vars)
        ilpfile.write('Subject to\n')
        write_minimal_precursor_source_constraints(ilpfile,H,sourcedict,sources)
        if sourcenum != None:
            write_source_num_constraint(ilpfile,H,sourcedict,sources,sourcenum)
    elif nonparametric == True:
        sourcedict = write_minimal_precursor_objective(ilpfile,sources,justgetsourcedict=True)
        write_maximize_flux_objective(ilpfile)
        ilpfile.write('Subject to\n')
        write_target_constraints(ilpfile,H,targets,stoichiometries,'psilon')
    else:
        write_stoichiometry_objective(ilpfile,H,source_edges,binary_vars=binary_vars)
        ilpfile.write('Subject to\n')
    write_vertex_flux_constraints(ilpfile,H,stoichiometries,targets,positiveflux=positiveflux,sources=sources)
    write_edge_flux_constraints(ilpfile,H,targets)
    if negative_regulators != None:
        write_negative_regulator_constraints(ilpfile,H,negative_regulators)
    #epsilon = find_epsilon(H,stoichiometries,targets,source_edges,positiveflux=positiveflux,sources=sources)
    #epsilon = .000005
    epsilon = .005
    #epsilon = .03
    #epsilon =  .0015
    #epsilon = .5
    if nonparametric != True:
        write_target_constraints(ilpfile,H,targets,stoichiometries,epsilon)
    if reversible != None:
        write_reversible_constraints(ilpfile,H,reversible)
    if binary_vars == True:
        write_binary_edge_constraints(ilpfile,H)
        if sourcenum != None:
            write_binary_source_constraints(ilpfile,sources,sourcedict)
    ilpfile.write('End')
    if sourcevars == True:
        return sourcedict

def write_minimal_precursor_ilp(ilpfile,H,stoichiometries,targets,source_edges,positiveflux=False,sources=[],binary_vars=True,negative_regulators=None,reversible=False,edgenum=None):
    sourcedict = write_minimal_precursor_objective(ilpfile,sources)
    ilpfile.write('Subject to\n')
    write_minimal_precursor_source_constraints(ilpfile,H,sourcedict,sources)
    write_vertex_flux_constraints(ilpfile,H,stoichiometries,targets,positiveflux=positiveflux,sources=sources)
    write_edge_flux_constraints(ilpfile,H,targets)
    if edgenum != None:
        write_edge_num_constraint(ilpfile,H,edgenum)
    if negative_regulators != None:
        write_negative_regulator_constraints(ilpfile,H,negative_regulators)
    #epsilon = find_epsilon(H,stoichiometries,targets,source_edges,positiveflux=positiveflux,sources=sources)
    #epsilon = .000005
    epsilon = .005
    write_target_constraints(ilpfile,H,targets,stoichiometries,epsilon)
    if reversible == True:
        write_reversible_constraints(ilpfile,H)
    if binary_vars == True:
        write_binary_edge_constraints(ilpfile,H)
        write_binary_source_constraints(ilpfile,sources,sourcedict)
    ilpfile.write('End')
    print('wrote mp file')

def write_minimal_precursor_source_constraints(ilpfile,H,sourcedict,sources):
    for so in sources:
        for e in H.get_forward_star(so):
            ilpfile.write('{} - {} >= 0\n'.format(sourcedict[so],f(e)))
    
def write_source_num_constraint(ilpfile,H,sourcedict,sources,sourcenum):
    for so in sources:
        ilpfile.write(' + {} '.format(sourcedict[so]))
    ilpfile.write(' <= {}\n'.format(sourcenum))

def write_edge_num_constraint(ilpfile,H,edgenum):
    for e in H.hyperedge_id_iterator():
        ilpfile.write(' + {} '.format(b(e)))
    ilpfile.write(' <= {}\n'.format(edgenum))

def write_minimal_precursor_objective(ilpfile,sources,justgetsourcedict=False):
    if justgetsourcedict == False:
        ilpfile.write('Minimize\n')
    sourcedict = {}
    sourcenum = 0
    for so in sources:
        sourcedict[so] = s(sourcenum)
        if justgetsourcedict == False:
            ilpfile.write(' + {} '.format(sourcedict[so]))
        sourcenum += 1
    if justgetsourcedict == False:
        ilpfile.write('\n')
    return sourcedict

def write_binary_source_constraints(ilpfile,sources,sourcedict):
    #ilpfile.write('Binary\n')
    for so in sources:
        ilpfile.write('{}\n'.format(sourcedict[so]))

def write_negative_regulator_constraints(ilpfile,H,negative_regulators):
    print('starting negative regulators')
    node_set = H.get_node_set()
    for e in H.hyperedge_id_iterator():
        if e in negative_regulators:
            for v in negative_regulators[e]:
                if v in node_set:
                    for f in H.get_backward_star(v):
                        ilpfile.write('{} + {} <= 1\n'.format(b(e),b(f)))


def write_vertex_flux_constraints(ilpfile,H,stoichiometries,targets,positiveflux=False,source_edges=None,sources=[],aedge=False):
    #print(targets)
    for v in H.node_iterator():
        if v not in targets and v not in sources:
            if source_edges == None:
                backedges = H.get_backward_star(v)
                forwardedges = H.get_forward_star(v)
                if len(backedges) > 0 or len(forwardedges) > 0:
                    #add all the head vertices
                    #ilpfile.write('{}: '.format(v))
                    for e in backedges:
                        if aedge == False:
                            ilpfile.write('+ {} {} '.format(stoichiometries[e][v],f(e)))
                        else:
                            ilpfile.write('+ {} {} '.format(stoichiometries[e][v],a(e)))
                    #subtract all the tail vertices
                    for e in forwardedges:
                        if aedge == False:
                            ilpfile.write('- {} {} '.format(stoichiometries[e][v],f(e)))
                        else:
                            ilpfile.write('- {} {} '.format(stoichiometries[e][v],a(e)))

                    if positiveflux == False:
                        ilpfile.write('= 0\n')
                    else:
                        ilpfile.write('>= 0\n')


def write_edge_flux_constraints(ilpfile,H,targets,no_upper_bound=False):
    targetedges = set()
    for t in targets:
        newtargetedges = set(H.get_backward_star(t))
        targetedges.update(newtargetedges)
    for e in H.hyperedge_id_iterator():
        if e not in targetedges:
            ilpfile.write('{} >= 0\n'.format(f(e)))
            if no_upper_bound == False:
                ilpfile.write('{} <= 1\n'.format(f(e)))

def write_binary_edge_constraints(ilpfile,H):
    for e in H.hyperedge_id_iterator():
        ilpfile.write('{} - {} >= 0\n'.format(b(e),f(e)))
    ilpfile.write('Binary\n')
    for e in H.hyperedge_id_iterator():
        ilpfile.write('{}\n'.format(b(e)))


def write_target_constraints(ilpfile,H,targets,stoichiometries,epsilon,aedge=False):
    #print(targets)
    #print('target constraints')
    for target in targets:
        for i,e in enumerate(H.get_backward_star(target)):
            #if i == 0:
                #ilpfile.write('{} {} '.format(stoichiometries[e][target],f(e)))
            #else:
            if aedge == False:
                ilpfile.write(' + {} {} '.format(stoichiometries[e][target],f(e)))
            else:
                ilpfile.write(' + {} {} '.format(stoichiometries[e][target],a(e)))
    if isinstance(epsilon,str):
        ilpfile.write(' - {}  >= 0\n'.format(epsilon))
    else:
        ilpfile.write(' >= {}\n'.format(epsilon))
    #ilpfile.write(' > {}\n'.format(0))

def find_epsilon(H,stoichiometries,targets,source_edges,positiveflux=False,sources=[]):
    #here we solve essentially the same ILP but have a different objective function
    epsilon_f_name = 'epsilon.ilp'
    epsilon_f = open(epsilon_f_name,'w')
    write_epsilon_ilp(epsilon_f,H,stoichiometries,targets,source_edges,positiveflux=positiveflux,sources=sources)
    epsilon_f.write('End')
    epsilon_f.close()
    maxsources = solve_epsilon_ilp(epsilon_f_name)
    epsilon_f2_name = 'operation.lp'
    epsilon_f2 = open(epsilon_f2_name,'w')
    epsilon = epsilon_operation(epsilon_f2,H,stoichiometries,targets,source_edges,epsilon_f2_name,positiveflux=positiveflux,sources=sources)
    print('epsilon',epsilon)
    return epsilon

def epsilon_operation(ilpfile,H,stoichiometries,targets,source_edges,ilpfilename,positiveflux=False,sources=[]):
    #write_epsilon_objective(ilpfile,H,source_edges)
    ilpfile.write('Maximize f_e1\n')
    ilpfile.write('Subject to\n')
    write_edge_flux_constraints(ilpfile,H,targets,no_upper_bound=True)
    write_vertex_flux_constraints(ilpfile,H,stoichiometries,targets,positiveflux=positiveflux,sources=sources)
    print('wrote vertex flux constraints')
    write_epsilon_target_constraints(ilpfile,H,targets,stoichiometries)
    ilpfile.close()
    ilp = cplex.Cplex()
    ilp.read(ilpfilename)
    ilp.objective.set_linear('f_e1',0)
    maxmaxsource = 0
    for e in H.hyperedge_id_iterator():
        if e not in sources and e not in targets:
            ilp.objective.set_linear(f(e),1)
            ilp.solve()
            maxsource = ilp.solution.get_objective_value()
            print('maxsource',maxsource)
            if maxsource > maxmaxsource:
                maxmaxsource = maxsource
            ilp.objective.set_linear(f(e),0)
    print('maxmaxsource',maxmaxsource)
    return 1/maxmaxsource

def write_epsilon_ilp(ilpfile,H,stoichiometries,targets,source_edges,positiveflux=False,sources=[]):
    write_epsilon_objective(ilpfile,H,source_edges)
    ilpfile.write('Subject to\n')
    write_edge_flux_constraints(ilpfile,H,targets,no_upper_bound=True)
    write_vertex_flux_constraints(ilpfile,H,stoichiometries,targets,positiveflux=positiveflux,sources=sources)
    print('wrote vertex flux constraints')
    write_epsilon_target_constraints(ilpfile,H,targets,stoichiometries)

def write_epsilon_target_constraints(ilpfile,H,targets,stoichiometries):
    for target in targets:
        print('target',target)
        print(H.get_backward_star(target))
        for e in H.get_backward_star(target):
            ilpfile.write('+ {} {} '.format(stoichiometries[e][target],f(e)))
    ilpfile.write(' = 1\n')

def write_epsilon_objective(ilpfile,H,source_edges):
    #Objective function is to maximize the flux of the source vertices, by maximizing the flux of the source_edges
    ilpfile.write('Maximize\n')
    for e in source_edges:
        ilpfile.write('+ {} '.format(f(e)))
    ilpfile.write('\n')

def write_reversible_constraints(ilpfile,H,reversible_reactions):
    processed_edges = set()
    for e in reversible_reactions:
        if e not in processed_edges:
            processed_edges.add(e)
            g = reversible_reactions[e]
            processed_edges.add(g)
            ilpfile.write('{} + {} <= 1\n'.format(b(e),b(g)))

def find_source_edges(H,sources,add_edges=False,stoichiometries=None,negative_regulators=None):
    #This function should either add in the source edges or just compile them
    source_edges = set()
    if add_edges == False:
        for source in sources:
            for e in H.get_forward_star(source):
                source_edges.add(e)
        return list(source_edges)
    else:
        for source in sources:
            eid = H.add_hyperedge(set(),set([source]))
            if stoichiometries != None:
                stoichiometries[eid] = {}
                stoichiometries[eid][source] = 1.0
            if negative_regulators != None:
                negative_regulators[eid] = 'None'
            source_edges.add(eid)
        return list(source_edges),H,stoichiometries,negative_regulators
        
def solve_epsilon_ilp(f):
    ilp = cplex.Cplex()
    ilp.read(f)
    ilp.solve()
    maxsource = ilp.solution.get_objective_value()
    print('maxsource',maxsource)
    #TODO: This is probably not exactly right
    epsilon = 1.0 / maxsource
    return epsilon

def solve_minimal_precursor_ilp(f,H,binary_vars=True):
    ilp = cplex.Cplex()
    ilp.read(f)
    ilp.solve()
    numsolsfound = 1
    outprefix = 'precursorilp'
    if binary_vars == True:
        number_of_sources = ilp.solution.pool.get_objective_value(0)
        ilp.solution.pool.write('%s-%d.sol' % (outprefix,numsolsfound),0)
    else:
        number_of_sources = ilp.solution.get_objective_value()
        ilp.solution.write('%s-%d.sol' % (outprefix,numsolsfound))
    variables = getILPSolution(H,[],outprefix,numsolsfound,number_of_sources,False)
    #ones = [var for var in variables if variables[var] == 1 and var[0]=='b']
    sourceones = [var for var in variables if variables[var] > .2 and var[0]=='s']
    ones = [var for var in variables if variables[var] > 0 and var[0]=='f']
    print(ones,'ones','lenones',len(ones))
    print(sourceones,'sourceones','lensourceones',len(sourceones))
    print('\n')
    return ones, sourceones

def initialize_nonparametric_ilp(ilpfile, cplexlog):
    ilp = cplex.Cplex()
    ilp.set_results_stream(cplexlog)
    ilp.set_warning_stream(cplexlog)
    ilp.set_error_stream(cplexlog)
    ilp.set_log_stream(cplexlog)
    ilp.read(ilpfile)
    return ilp

def check_nonparametric_reachability(ilp):
    ilp.solve()
    if ilp.solution.get_objective_value() <= 0.00001 or ilp.solution.pool.get_num() <= 0 or ilp.solution.pool.get_objective_value(0) <= 0.00001 or ilp.solution.get_objective_value() <= 0.00001:
        #This means the solution was infeasible
        print('solution infeasible')
        print(ilp.solution.get_status_string())
        print(ilp.solution.pool.get_objective_value(0))
        print('raising exception for nonreachability')
        raise Exception
    else:
        print(ilp.solution.pool.get_objective_value(0))

def find_nonparametric_factory(nonparailpfname, ilpfname, H, stoichiometries, targets, source_edges, sources, solnum=0, positiveflux=True, reversible=None, extra_stoichiometries=None):

    #write the nonparametric file
    nonparailpfile = open(nonparailpfname,'w')
    write_stoichiometry_ilp(nonparailpfile,H,stoichiometries,targets,source_edges,sources=sources,positiveflux=positiveflux,binary_vars=True,nonparametric=True, reversible=reversible)
    nonparailpfile.close()
    cplexlog = open(nonparailpfname + 'cplex.log','w')
    nonparailp = initialize_nonparametric_ilp(nonparailpfname, cplexlog)
    try:
        check_nonparametric_reachability(nonparailp)
    except:
        raise Exception


    #reachability check passed, so run Odinn to find a minedge factory
    odinnstart = time.time()
    ilpfile = open(ilpfname,'w')
    write_stoichiometry_ilp(ilpfile,H,stoichiometries,targets,source_edges,sources=sources,positiveflux=positiveflux,binary_vars=True, reversible=reversible)
    print('written ilp')
    ilpfile.close()
    try:
        ilp_soln, active_edges = solve_stoichiometry_ilp_2(ilpfname,H,sources,binary_vars=True,solnum=solnum, cplexlog=cplexlog)
        #ilp_soln, active_edges = solve_stoichiometry_ilp_2(ilpfname,H,sources,binary_vars=True,solnum=solnum)
        odinnend = time.time()
        print(f'time for odinn: {odinnend - odinnstart}')
        len_active_edges = len(active_edges)
        print('trying to decrement')
        if len_active_edges == 1:
            decremented_exists = False
        elif len_active_edges == 0: # this is the case where no factory really exists under the epsilon
            decremented_exists = True
        else:
            decremented_exists = solve_decremented_ilp_nonparametric(nonparailp, H, len_active_edges)
    except:
        #means that no factory was found with Odinn
        print(traceback.format_exc())
        print('no factory found with Odinn')
        len_active_edges = None
        # set the flag indicating that a smaller factory exists so that the full nonparametric algorithm is run
        decremented_exists = True
    if decremented_exists == False:
        print('no decremented factory exists, returning Odinn result')
        solved_time = time.time()
        ones = solve_nonparametric_with_known_k(nonparailp,H,sources,len_active_edges,solnum=solnum)
        #import pickle as pkl
        #node_dict = pkl.load(open('node_dict.pkl','rb')) 
        #for e in ones:
            #e = unb(e)
            #print(f'edge: {e}:\ntail: {[(node_dict[v],stoichiometries[e][v],extra_stoichiometries[e][v]) for v in H.get_hyperedge_tail(e)]}\nhead: {[(node_dict[w], stoichiometries[e][w], extra_stoichiometries[e][w]) for w in H.get_hyperedge_head(e)]}\n')
        return solved_time
    else:
        nonparailp = initialize_nonparametric_ilp(nonparailpfname, cplexlog)
        ilpres, ones = solve_full_nonparametric_ilp(nonparailp,H,sources,solnum=solnum)
        #import pickle as pkl
        #node_dict = pkl.load(open('node_dict.pkl','rb')) 
        #for e in ones:
            #e = unb(e)
            #print(f'edge: {e}:\ntail: {[(node_dict[v],stoichiometries[e][v]) for v in H.get_hyperedge_tail(e)]}\nhead: {[(node_dict[w], stoichiometries[e][w]) for w in H.get_hyperedge_head(e)]}\n')
        solved_time = time.time()
        return solved_time


def solve_decremented_ilp_nonparametric(ilp, H, starting_k):
    alledges = [b(e) for e in H.get_hyperedge_id_set()]
    k = starting_k - 1
    val = [1.0] * len(alledges)
    eq = cplex.SparsePair(ind=alledges,val=val)
    ilp.linear_constraints.add(lin_expr = [eq], senses = ['L'], rhs = [k], names = ['k{0}'.format(k)])
    ilp.solve()
    maxepsilon = ilp.solution.get_objective_value()
    print(k,maxepsilon)
    ilp.linear_constraints.delete(f'k{k}')
    if maxepsilon < .00001:
        print('decrementing did not find a better solution')
        return False
    else:
        return True

def solve_nonparametric_with_known_k(ilp,H,sources,k,solnum=0):
    alledges = [b(e) for e in H.get_hyperedge_id_set()]
    falledges = [f(e) for e in H.get_hyperedge_id_set()]
    val = [1.0] * len(alledges)
    eq = cplex.SparsePair(ind=alledges,val=val)
    ilp.linear_constraints.add(lin_expr = [eq], senses = ['L'], rhs = [k], names = ['k{0}'.format(k)])
    ilp.solve()
    maxepsilon = ilp.solution.get_objective_value()
    print(k,maxepsilon)
    variables = ilp.solution.get_values(alledges)
    fvariables = ilp.solution.get_values(falledges)
    ones = [alledges[i] for i in range(len(variables)) if variables[i] > 0]
    fones = [falledges[i] for i in range(len(fvariables)) if fvariables[i] > 0]
    fvals = [fvariables[i] for i in range(len(fvariables)) if fvariables[i] > 0]
    sourceset = set()
    for one in ones:
        for v in H.get_hyperedge_tail(unb(one)):
            if v in sources:
                sourceset.add(v)
    ksolutionstring = f'{fones} ones\n{fvals}\n{ones} ones lenones {len(ones)}\n{sourceset} sourceones lensourceones {len(sourceset)}\n\n'
    print('nonparametric solution')
    print(ksolutionstring)
    return ones

def solve_full_nonparametric_ilp(ilp,H,sources,solnum=0):
    maxepsilon = 0
    alledges = [b(e) for e in H.get_hyperedge_id_set()]
    falledges = [f(e) for e in H.get_hyperedge_id_set()]
    k = 0
    while maxepsilon < .00001:
        if k != 0:
            k *= 2
        else:
            k = 1
        val = [1.0] * len(alledges)
        eq = cplex.SparsePair(ind=alledges,val=val)
        ilp.linear_constraints.add(lin_expr = [eq], senses = ['L'], rhs = [k], names = ['k{0}'.format(k)])
        ilp.solve()
        print(ilp.solution.get_status_string())
        maxepsilon = ilp.solution.get_objective_value()
        print(k,maxepsilon)
        if maxepsilon < .00001:
            ilp.linear_constraints.delete(f'k{k}')
        else:
            variables = ilp.solution.get_values(alledges)
            fvariables = ilp.solution.get_values(falledges)
            ones = [alledges[i] for i in range(len(variables)) if variables[i] > 0]
            fones = [falledges[i] for i in range(len(fvariables)) if fvariables[i] > 0]
            fvals = [fvariables[i] for i in range(len(fvariables)) if fvariables[i] > 0]
            sourceset = set()
            for one in ones:
                for v in H.get_hyperedge_tail(unb(one)):
                    if v in sources:
                        sourceset.add(v)
            ksolutionstring = f'{fones} ones\n{fvals}\n{ones} ones lenones {len(ones)}\n{sourceset} sourceones lensourceones {len(sourceset)}\n\n'
    print('hit maximum, searching inbetween')
    return search_range_nonparametric_ilp(k // 2, k, ilp, H, sources, (ilp, ones), ksolutionstring)

def search_range_nonparametric_ilp(lowerk, higherk, ilp, H, sources, higherksolution, higherksolutionstring):
    #epsilon is >0 at lower k, but <=0 at higherk
    k = lowerk + (higherk - lowerk) // 2
    alledges = [b(e) for e in H.get_hyperedge_id_set()]
    val = [1.0] * len(alledges)
    eq = cplex.SparsePair(ind=alledges,val=val)
    ilp.linear_constraints.add(lin_expr = [eq], senses = ['L'], rhs = [k], names = ['k{0}'.format(k)])
    ilp.solve()
    print(ilp.solution.get_status_string())
    maxepsilon = ilp.solution.get_objective_value()
    print(k,maxepsilon)
    if maxepsilon < .00001: #no flux at this k
        if higherk - k == 1:
            # This means higherk is the lowest value of k for which there is a factory, so we should output its solution
            print(higherksolutionstring)
            return higherksolution
        else:
            ilp.linear_constraints.delete(f'k{k}')
            return search_range_nonparametric_ilp(k, higherk, ilp, H, sources, higherksolution, higherksolutionstring)
    else: #flux through at this k
        falledges = [f(e) for e in H.get_hyperedge_id_set()]
        variables = ilp.solution.get_values(alledges)
        fvariables = ilp.solution.get_values(falledges)
        ones = [alledges[i] for i in range(len(variables)) if variables[i] > 0]
        fones = [falledges[i] for i in range(len(fvariables)) if fvariables[i] > 0]
        fvals = [fvariables[i] for i in range(len(fvariables)) if fvariables[i] > 0]
        sourceset = set()
        for one in ones:
            for v in H.get_hyperedge_tail(unb(one)):
                if v in sources:
                    sourceset.add(v)
        ksolutionstring = f'{fones} ones\n{fvals}\n{ones} ones lenones {len(ones)}\n{sourceset} sourceones lensourceones {len(sourceset)}\n\n'
        ksolution = (ilp, ones)
        if k - lowerk == 1:
            print(ksolutionstring)
            return ksolution
        else:
            ilp.linear_constraints.delete(f'k{k}')
            return search_range_nonparametric_ilp(lowerk, k, ilp, H, sources, ksolution, ksolutionstring)

def solve_stoichiometry_ilp_2(fi,H,sources,binary_vars=True,solnum=0, cplexlog=None):
    ilp = cplex.Cplex()
    if cplexlog != None:
        ilp.set_results_stream(cplexlog)
        ilp.set_warning_stream(cplexlog)
        ilp.set_error_stream(cplexlog)
        ilp.set_log_stream(cplexlog)
    startread = time.time()
    ilp.read(fi)
    endread = time.time()
    print(f'time to read: {endread - startread}')
    ilp.solve()
    endsolve = time.time()
    print(f'time to solve: {endsolve - endread}')
    numsolsfound = 1

    alledges = [b(e) for e in H.get_hyperedge_id_set()]
    falledges = [f(e) for e in H.get_hyperedge_id_set()]
    variables = ilp.solution.get_values(alledges)
    fvariables = ilp.solution.get_values(falledges)
    ones = [alledges[i] for i in range(len(variables)) if variables[i] > 0.2]
    fones = [falledges[i] for i in range(len(fvariables)) if fvariables[i] > 0]
    fvals = [fvariables[i] for i in range(len(fvariables)) if fvariables[i] > 0]

    #outprefix = 'overallilp{}'.format(solnum)
    if binary_vars == True:
        number_of_edges = ilp.solution.pool.get_objective_value(0)
        #ilp.solution.pool.write('%s-%d.sol' % (outprefix,numsolsfound),0)
    else:
        number_of_edges = ilp.solution.get_objective_value()
        #ilp.solution.write('%s-%d.sol' % (outprefix,numsolsfound))
    timesolution = time.time()
    print(f'time to write solution: {timesolution - endsolve}')
    #variables = getILPSolution(H,[],outprefix,numsolsfound,number_of_edges,False)
    #ones = [var for var in variables if variables[var] > 0.2 and var[0]=='b']
    #fones = [var for var in variables if variables[var] > 0 and var[0]=='f']
    #fvals = [variables[v] for v in fones]
    sourceset = set()
    for one in ones:
        for v in H.get_hyperedge_tail(unb(one)):
            if v in sources:
                sourceset.add(v)
    print('Odinn results')
    print(fones,'ones')
    print(fvals)
    print(ones,'ones','lenones',len(ones))
    print(sourceset,'sourceones','lensourceones',len(sourceset))
    print('\n')
    return ilp, ones

def enumerate_minimal_factories(ilpname,H,stoichiometries,targets,source_edges,negative_regulators,positiveflux=False,sources=[],negative=False):
    ilpfile = open(ilpname,'w')
    if negative == False:
        sourcedict = write_stoichiometry_ilp(ilpfile,H,stoichiometries,targets,source_edges,sources=sources,positiveflux=positiveflux,binary_vars=True,sourcevars=True)
    else:
        sourcedict = write_stoichiometry_ilp(ilpfile,H,stoichiometries,targets,source_edges,sources=sources,positiveflux=positiveflux,binary_vars=True,sourcevars=True,negative_regulators=negative_regulators)
    ilpfile.close()
    ilp = cplex.Cplex()
    ilp.read(ilpname)
    validsolutionfound = False

def solve_stoichiometry_ilp_nonparametric(fil,H,sources,binary_vars=True,solnum=0,doubling=False,starting_k=None):
    ilp = cplex.Cplex()
    ilp.read(fil)
    maxepsilon = 0
    alledges = [b(e) for e in H.get_hyperedge_id_set()]
    falledges = [f(e) for e in H.get_hyperedge_id_set()]
    if starting_k != None:
        #decrement k and see if it is solvable
        k = starting_k
        maxepsilon = 1
        already_decremented = False
        while maxepsilon >= .00001:
            k -= 1
            val = [1.0] * len(alledges)
            eq = cplex.SparsePair(ind=alledges,val=val)
            ilp.linear_constraints.add(lin_expr = [eq], senses = ['L'], rhs = [k], names = ['k{0}'.format(k)])
            ilp.solve()
            print(ilp.solution.get_status_string())
            if binary_vars == True:
                #maxepsilon = ilp.solution.pool.get_objective_value(0)
                maxepsilon = ilp.solution.get_objective_value()
                print(k,maxepsilon)
            else:
                maxepsilon = ilp.solution.get_objective_value()
            if maxepsilon >= .00001:
                #store the current solution, since we might need to output it
                variables = ilp.solution.get_values(alledges)
                fvariables = ilp.solution.get_values(falledges)
                ones = [alledges[i] for i in range(len(variables)) if variables[i] > 0]
                fones = [falledges[i] for i in range(len(fvariables)) if fvariables[i] > 0]
                fvals = [fvariables[i] for i in range(len(fvariables)) if fvariables[i] > 0]
                sourceset = set()
                for one in ones:
                    for v in H.get_hyperedge_tail(unb(one)):
                        if v in sources:
                            sourceset.add(v)
                ilp.linear_constraints.delete(f'k{k}')
                already_decremented = True
            elif already_decremented == True:
                #print the previous solution, since it was the last viable one
                print(fones,'ones')
                print(fvals)
                print(ones,'ones','lenones',len(ones))
                print(sourceset,'sourceones','lensourceones',len(sourceset))
                print('\n')
                return
            else:
                print('decrementing did not find a better solution')
                return
    else:
        k = 0
        ilp.solve()
        if ilp.solution.pool.get_num() <= 0 or ilp.solution.pool.get_objective_value(0) <= 0.00001 or ilp.solution.get_objective_value() <= 0.00001:
            #This means the solution was infeasible
            print('solution infeasible')
            print(ilp.solution.get_status_string())
            print(ilp.solution.pool.get_objective_value(0))
            raise Exception
        else:
            print(ilp.solution.pool.get_objective_value(0))
        #while maxepsilon < .00001 and k < 1:
        if doubling == True:
            k = 1
            while maxepsilon < .00001:
                k *= 2
                val = [1.0] * len(alledges)
                eq = cplex.SparsePair(ind=alledges,val=val)
                ilp.linear_constraints.add(lin_expr = [eq], senses = ['L'], rhs = [k], names = ['k{0}'.format(k)])
                ilp.solve()
                print(ilp.solution.get_status_string())
                if binary_vars == True:
                    #maxepsilon = ilp.solution.pool.get_objective_value(0)
                    maxepsilon = ilp.solution.get_objective_value()
                    print(k,maxepsilon)
                else:
                    maxepsilon = ilp.solution.get_objective_value()
                ilp.linear_constraints.delete(f'k{k}')
            print('hit maximum, incrementing from previous')
            k /= 2
            maxepsilon = 0
            while maxepsilon < .00001:
                k += 1
                val = [1.0] * len(alledges)
                eq = cplex.SparsePair(ind=alledges,val=val)
                ilp.linear_constraints.add(lin_expr = [eq], senses = ['L'], rhs = [k], names = ['k{0}'.format(k)])
                ilp.solve()
                print(ilp.solution.get_status_string())
                if binary_vars == True:
                    #maxepsilon = ilp.solution.pool.get_objective_value(0)
                    maxepsilon = ilp.solution.get_objective_value()
                    print(k,maxepsilon)
                else:
                    maxepsilon = ilp.solution.get_objective_value()
                if maxepsilon < .00001:
                    ilp.linear_constraints.delete(f'k{k}')
            '''
            while maxepsilon >= .00001:
                k -= 1
                val = [1.0] * len(alledges)
                eq = cplex.SparsePair(ind=alledges,val=val)
                ilp.linear_constraints.add(lin_expr = [eq], senses = ['L'], rhs = [k], names = ['k{0}'.format(k)])
                ilp.solve()
                print(ilp.solution.get_status_string())
                if binary_vars == True:
                    #maxepsilon = ilp.solution.pool.get_objective_value(0)
                    maxepsilon = ilp.solution.get_objective_value()
                    print(k,maxepsilon)
                else:
                    maxepsilon = ilp.solution.get_objective_value()
                ilp.linear_constraints.delete(f'k{k}')
            '''
            #k += 1 #Need to get k back to where it was a good value
            #val = [1.0] * len(alledges)
            #eq = cplex.SparsePair(ind=alledges,val=val)
            #ilp.linear_constraints.add(lin_expr = [eq], senses = ['L'], rhs = [k], names = ['k{0}'.format(k)])
            #ilp.solve()
                
        else:
            while maxepsilon < .00001:
                k += 1
                val = [1.0] * len(alledges)
                eq = cplex.SparsePair(ind=alledges,val=val)
                ilp.linear_constraints.add(lin_expr = [eq], senses = ['L'], rhs = [k], names = ['k{0}'.format(k)])
                ilp.solve()
                print(ilp.solution.get_status_string())
                if binary_vars == True:
                    #maxepsilon = ilp.solution.pool.get_objective_value(0)
                    maxepsilon = ilp.solution.get_objective_value()
                    print(k,maxepsilon)
                else:
                    maxepsilon = ilp.solution.get_objective_value()
                if maxepsilon < .00001:
                    ilp.linear_constraints.delete(f'k{k}')
        outprefix = 'nonparailp{}'.format(solnum)
        if binary_vars == True:
            #variables = ilp.solution.pool.get_values(0,alledges)
            #fvariables = ilp.solution.pool.get_values(0,falledges)
            variables = ilp.solution.get_values(alledges)
            fvariables = ilp.solution.get_values(falledges)
            ones = [alledges[i] for i in range(len(variables)) if variables[i] > 0]
            fones = [falledges[i] for i in range(len(fvariables)) if fvariables[i] > 0]
            fvals = [fvariables[i] for i in range(len(fvariables)) if fvariables[i] > 0]
            sourceset = set()
            for one in ones:
                for v in H.get_hyperedge_tail(unb(one)):
                    if v in sources:
                        sourceset.add(v)
            print(fones,'ones')
            print(fvals)
            print(ones,'ones','lenones',len(ones))
            print(sourceset,'sourceones','lensourceones',len(sourceset))
            print('\n')
        else:
            #TODO This will not work right now. Leaving it in because I'm not sure why it was there in the first place
            ilp.solution.write('%s-%d.sol' % (outprefix,k))
        #variables = getILPSolution(H,[],outprefix,1,maxepsilon,False)
        #ones = [var for var in variables if variables[var] > 0.2 and var[0]=='b']
        #fones = [var for var in variables if variables[var] > 0 and var[0]=='f']
        #fvals = [variables[v] for v in fones]
        return ilp, ones

def solve_stoichiometry_ilp(f,H,sources,binary_vars=True,solnum=0):
    ilp = cplex.Cplex()
    ilp.read(f)
    ilp.solve()
    numsolsfound = 1
    outprefix = 'overallilp{}'.format(solnum)
    if binary_vars == True:
        number_of_edges = ilp.solution.pool.get_objective_value(0)
        ilp.solution.pool.write('%s-%d.sol' % (outprefix,numsolsfound),0)
    else:
        number_of_edges = ilp.solution.get_objective_value()
        ilp.solution.write('%s-%d.sol' % (outprefix,numsolsfound))
    variables = getILPSolution(H,[],outprefix,numsolsfound,number_of_edges,False)
    ones = [var for var in variables if variables[var] > 0.2 and var[0]=='b']
    fones = [var for var in variables if variables[var] > 0 and var[0]=='f']
    fvals = [variables[v] for v in fones]
    sourceset = set()
    for one in ones:
        for v in H.get_hyperedge_tail(unb(one)):
            if v in sources:
                sourceset.add(v)
    print('normal factory ilp results')
    print(fones,'ones')
    print(fvals)
    print(ones,'ones','lenones',len(ones))
    print(sourceset,'sourceones','lensourceones',len(sourceset))
    print('\n')
    return ilp, ones

def enumerate_minimal_factories(ilpname,H,stoichiometries,targets,source_edges,negative_regulators,positiveflux=False,sources=[],negative=False):
    ilpfile = open(ilpname,'w')
    if negative == False:
        sourcedict = write_stoichiometry_ilp(ilpfile,H,stoichiometries,targets,source_edges,sources=sources,positiveflux=positiveflux,binary_vars=True,sourcevars=True)
    else:
        sourcedict = write_stoichiometry_ilp(ilpfile,H,stoichiometries,targets,source_edges,sources=sources,positiveflux=positiveflux,binary_vars=True,sourcevars=True,negative_regulators=negative_regulators)
    ilpfile.close()
    ilp = cplex.Cplex()
    ilp.read(ilpname)
    validsolutionfound = False
    pastsolution = []
    i = 0
    knockoutedge = ''
    sourceedgesused = []
    sourcesused = []
    firstlen = 0
    while validsolutionfound == False:
        #find the next solution
        print('\nfinding new solution\n')
        ilp,pastsolution,sourcesused,ilpsourceedgesused = get_next_stoichiometry_solution(ilp,H,sources,ones=pastsolution)
        if i == 0:
            firstlen = len(pastsolution)
        else:
            if len(pastsolution) > firstlen:
                return pastsolution,sourcesused,i
        i += 1
        print('iterationnumis ',i)
        if ilp == 'infeasible':
            return '','',i
    return pastsolution,sourcesused,i

def second_order_negative_regulation(ilpname,H,stoichiometries,targets,source_edges,negative_regulators,positiveflux=False,sources=[],randnum=1):
    #initially needs to setup the ILP object
    ilpfile = open(ilpname,'w')
    sourcedict = write_stoichiometry_ilp(ilpfile,H,stoichiometries,targets,source_edges,sources=sources,positiveflux=positiveflux,binary_vars=True,sourcevars=True,negative_regulators=negative_regulators,minimizesources=True)
    ilpfile.close()
    ilp = cplex.Cplex()
    ilp.read(ilpname)
    #while loop that goes in two steps, 
    validsolutionfound = False
    pastsolution = []
    i = 0
    knockoutedge = ''
    sourceedgesused = []
    sourcesused = []
    knockoutsolutions = []
    while validsolutionfound == False:
        #find the next solution
        print('\nfinding new solution\n')
        #ilp,pastsolution,sourcesused,ilpsourceedgesused = get_next_stoichiometry_solution(ilp,H,sources,ones=pastsolution,knockoutedge=knockoutedge,sourceedgesused=sourceedgesused)
        #ilp,pastsolution,sourcesused,ilpsourceedgesused = get_next_stoichiometry_solution(ilp,H,sources,ones=pastsolution,knockoutedge=knockoutedge,sourcesused=sourcesused,sourcedict=sourcedict,randnum=randnum)
        ilp,pastsolution,sourcesused,ilpsourceedgesused = get_next_stoichiometry_solution(ilp,H,sources,ones=pastsolution,knockoutsolutions=knockoutsolutions,sourcedict=sourcedict,randnum=randnum)
        i += 1
        print('iterationnumis ',i)
        if ilp == 'infeasible':
            return '','',i
        else:
            #validsolutionfound,knockoutedge, sourceedgesused = verify_second_order_negative_regulation(H,pastsolution,stoichiometries,negative_regulators,targets,source_edges,sourcesused,positiveflux=positiveflux,sources=sources,randnum=randnum)
            validsolutionfound,knockoutsolutions = verify_second_order_negative_regulation(H,pastsolution,stoichiometries,negative_regulators,targets,source_edges,sourcesused,positiveflux=positiveflux,sources=sources,randnum=randnum)
    return pastsolution,sourcesused,i


def verify_second_order_negative_regulation(H,solution,stoichiometries,negative_regulators,targets,source_edges,solutionsources,positiveflux=True,sources=[],randnum=1):
    #First compile a list hyperedges that make negative regulators of the path hyperedges
    '''
    negreghyperedges = set()
    for be in solution:
        e = unb(be)
        if e in negative_regulators:
            negregs = negative_regulators[e]
            for v in negregs:
                negreghyperedges.update(H.get_backward_star(v))
        
    #for each hyperedge, try to maximize the flux through it while making one unit of target flux
    print(negreghyperedges)
    ilpname = '2Onegreg.ilp'
    for e in negreghyperedges:
        print('solving new negreg hyperedge feasibility ilp for hyperedge {}'.format(e))
        ilpfile = open(ilpname,'w')
        write_second_order_negreg_ilp(ilpfile,H,stoichiometries,targets,source_edges,e,solutionsources,positiveflux=positiveflux,sources=sources)
        ilpfile.close()
        feasible = solve_second_order_negreg_ilp(ilpname)
        print('done solving new negreg hyperedge feasibility ilp for hyperedge {}, got {}'.format(e,feasible))
        if feasible == True:
            return False,e
    #This means that none of the negative regulators were reachable
    return True,''
    '''
    negreghyperedges = set()
    feasiblesolutions = []
    for be in solution:
        e = unb(be)
        if e in negative_regulators:
            negregs = negative_regulators[e]
            for v in negregs:
                if v in H.get_node_set():
                    for negedge in H.get_backward_star(v):
                        if negedge not in negreghyperedges:
                            negreghyperedges.add(negedge)
                            ilpname = '2Onegreg.ilp'
                            print('solving new negreg hyperedge feasibility ilp for hyperedge {}'.format(negedge))
                            ilpfile = open(ilpname,'w')
                            write_second_order_negreg_ilp(ilpfile,H,stoichiometries,targets,source_edges,negedge,solutionsources,positiveflux=positiveflux,sources=sources)
                            ilpfile.close()
                            feasible, sourceedgesused = solve_second_order_negreg_ilp(ilpname,H,sources,randnum=randnum)
                            #print('done solving new negreg hyperedge feasibility ilp for hyperedge {}, got {}, which made negative regulator {} with hyperedges {}'.format(e,feasible,v,[(H.get_hyperedge_tail(r),H.get_hyperedge_head(r)) for r in sourceedgesused]))
                            if feasible == True:
                                #return False,e, sourceedgesused
                                feasiblesolutions.append((e,sourceedgesused))
    #This means that none of the negative regulators were reachable
    if len(feasiblesolutions) != 0:
        return False,feasiblesolutions
    else:
        return True,feasiblesolutions
        

def solve_second_order_negreg_ilp(f,H,sources,randnum=1):
    ilp = cplex.Cplex()
    ilp.read(f)
    ilp.solve()
    statusnum = ilp.solution.get_status()
    print('statusnum',statusnum)
    #statusnums are 1 = optimal 2 = unbounded    3 = infeasible  4 = inf or unbounded
    if statusnum == 3:
        #This means it is infeasible
        return False,set()
    else:
        ilp.solution.write('sonegreg{}-1.sol'.format(randnum))
        variables = getILPSolution(H,[],'sonegreg{}'.format(randnum),1,ilp.solution.get_objective_value(),False)
        ones = [var for var in variables if variables[var] > .001 and var[0]=='f']
        vals = [variables[v] for v in ones]
        sourceset = set()
        sourceedgesused = set()
        for one in ones:
            for v in H.get_hyperedge_tail(unf(one)):
                if v in sources:
                    sourceset.add(v)
                    sourceedgesused.add(unf(one))
        #print(ones,vals,[(H.get_hyperedge_tail(unb(o)),H.get_hyperedge_head(unb(o))) for o in ones])
        #print('sources used', sourceset)
        return True, sourceedgesused
    
def write_second_order_negreg_ilp(ilpfile,H,stoichiometries,targets,source_edges,targetedge,allowedsources,positiveflux=False,sources=[]):
    '''
    ilpfile.write('Maximize\n')
    ilpfile.write('{} \n'.format(f(targetedge)))
    ilpfile.write('Subject to\n')
    ilpfile.write('{} <= 100000\n'.format(f(targetedge)))
    write_edge_flux_constraints(ilpfile,H,targets,no_upper_bound=True)
    write_vertex_flux_constraints(ilpfile,H,stoichiometries,targets,positiveflux=positiveflux,sources=sources)
    write_epsilon_target_constraints(ilpfile,H,targets,stoichiometries)
    write_forbidden_source_constraints(ilpfile,H,allowedsources,source_edges,sources)
    ilpfile.write('\nEnd\n')
    '''
    allowed_edges = set()
    for e in source_edges:
        for v in H.get_hyperedge_tail(e):
            if v in sources and v not in allowedsources:
                #means v is a forbidden source
                break
        else:
            #means we went through and all tail vertices that were sources were allowed sources
            allowed_edges.add(e)

    ilpfile.write('Minimize\n')
    #for e in H.hyperedge_id_iterator():
    for e in allowed_edges:
        ilpfile.write(' + {} '.format(f(e)))
    ilpfile.write('\nSubject to\n')
    ilpfile.write('{} >= .005\n'.format(f(targetedge)))
    write_edge_flux_constraints(ilpfile,H,targets,no_upper_bound=False)
    write_vertex_flux_constraints(ilpfile,H,stoichiometries,targets,positiveflux=positiveflux,sources=sources)
    write_forbidden_source_constraints(ilpfile,H,allowedsources,source_edges,sources)
    epsilon = .005
    write_target_constraints(ilpfile,H,targets,stoichiometries,epsilon)
    ilpfile.write('\nEnd\n')


def write_forbidden_source_constraints(ilpfile,H,allowedsources,source_edges,sources):
    allowed_edges = set()
    for e in source_edges:
        for v in H.get_hyperedge_tail(e):
            if v in sources and v not in allowedsources:
                break
        else:
            #means we went through and all of the vertices were either not sources or were allowed
            allowed_edges.add(e)
    forbidden_edges = set(source_edges).difference(allowed_edges)
    for e in forbidden_edges:
        ilpfile.write('{} = 0\n'.format(f(e)))
    return allowed_edges

def get_next_stoichiometry_solution(ilp,H,sources,ones=[],knockoutedge='',sourceedgesused=[],sourcesused=[],sourcedict={},randnum=1,knockoutsolutions=[]):
    #add in new constraint and solve
    #ones already has a list of the binary variables that we need to exclude
    print('getting next stoichiometry solution')
    if ones != []:
        val = [1.0] * len(ones)
        eq = cplex.SparsePair(ind=ones,val=val)
        ilp.linear_constraints.add(lin_expr = [eq], senses = ['L'], rhs = [len(ones)-1])

    if knockoutedge != '' and sourceedgesused != []:
        print('knockoutedge',knockoutedge)
        print('sourceedgesused',sourceedgesused)
        #this constraint is just to make sure we aren't repeating a solution we know won't work. 
        #So we can't have the knockout edge if the source hyperedges are used
        sourceedgesused.add(knockoutedge)
        sourceedgesused = [b(e) for e in sourceedgesused]
        val = [1.0] * len(sourceedgesused)
        eq = cplex.SparsePair(ind=sourceedgesused,val=val)
        ilp.linear_constraints.add(lin_expr = [eq], senses = ['L'], rhs = [len(sourceedgesused)-1])

    if knockoutedge != '' and sourcesused != []:
        print('knockoutedge',knockoutedge)
        print('sourcesused',sourcesused)
        #this constraint is just to make sure we aren't repeating a solution we know won't work. 
        #So we can't have the knockout edge if the source hyperedges are used
        excludevars = [sourcedict[so] for so in sourcesused] + [b(knockoutedge)]
        val = [1.0] * len(excludevars)
        eq = cplex.SparsePair(ind=excludevars,val=val)
        ilp.linear_constraints.add(lin_expr = [eq], senses = ['L'], rhs = [len(excludevars)-1])

    if knockoutsolutions != []:
        for knockoutedge,sourceedgesused in knockoutsolutions:
            #find the sources that are used from the source edges used
            sourcesused = set()
            for edge in sourceedgesused:
                for v in H.get_hyperedge_tail(edge):
                    if v in sources:
                        sourcesused.add(v)
            excludevars = [sourcedict[so] for so in sourcesused] + [b(knockoutedge)]
            print('\n\nexcludevars: {}\n\n'.format(excludevars))
            val = [1.0] * len(excludevars)
            eq = cplex.SparsePair(ind=excludevars,val=val)
            ilp.linear_constraints.add(lin_expr = [eq], senses = ['L'], rhs = [len(excludevars)-1])

    #extract new solution
    ilp.solve()
    numsolsfound = 1
    outprefix = '2onegreg{}'.format(randnum)
    if ilp.solution.pool.get_num() > 0:
        number_of_edges = ilp.solution.pool.get_objective_value(0)
        ilp.solution.pool.write('%s-%d.sol' % (outprefix,numsolsfound),0)
        variables = getILPSolution(H,[],outprefix,numsolsfound,number_of_edges,False)
        ones = [var for var in variables if variables[var] > .2 and var[0]=='b']
        print(ones)
        sourceset = set()
        sourceedgesused = set()
        for one in ones:
            for v in H.get_hyperedge_tail(unb(one)):
                if v in sources:
                    sourceset.add(v)
                    sourceedgesused.add(unb(one))
        return ilp,ones,sourceset,sourceedgesused
    else:
        return 'infeasible','','',[]

def getILPSolution(H,nodeset,outprefix,num,objective,verbose):
    #print ('\nGetting ILP Solution for Solution # %d in Pool' % (num))

    # parse xml
    xml = minidom.parse('%s-%d.sol' % (outprefix,num))
    cplexsol = xml.getElementsByTagName('CPLEXSolution')[0]
   
    # get variables
    elements = cplexsol.getElementsByTagName('variables')[0].getElementsByTagName('variable')
    variables = {}
    for v in elements:
        variables[str(v.getAttribute('name'))] = float(v.getAttribute('value'))

    out = open('%s-%d.variables' % (outprefix,num),'w')
    out.write('# Objective = %s\n' % (objective))
    out.write('#name\tval\trounded_val\n')

    if verbose == True:
        print ('VARIABLES:')
    
    numrounded = 0
    for v in variables:
        rounded = int(variables[v]+0.5)
        if verbose==True and variables[v] != 0:
            print((v,variables[v],rounded))
        out.write('%s\t%f\t%d\n' % (v,variables[v],rounded))
    out.close()
    
    #print (' wrote variables to file %s' % ('%s-%d.variables' % (outprefix,num)))

    return variables
    
def f(e):
    return('f_{}'.format(e))

def unf(e):
    return e[2:]

def unb(e):
    return e[2:]

def b(e):
    return('b_{}'.format(e))

def s(v):
    return('s_{}'.format(v))

def a(name):
    return 'a_%s' % (name)
