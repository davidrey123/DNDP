
# coding: utf-8

# In[ ]:

import time
import copy
import numpy as np
import docplex.mp
from docplex.mp.model import Model
import pyomo.environ as pyo
from pyomo.environ import SolverFactory
from pyomo.opt import SolverStatus, TerminationCondition


# In[ ]:

#---create instance from file

def read_instance(net,NDP,B_prop,m,scal_time,scal_flow,timelimit):
    
    #---read network and instance data in extended TNTP format
    nodes,links,capa,fftt,alpha,beta,links2,cost = read_network_data(net,NDP)
    
    #---read trip data in TNTP format
    OD,orig,dest = read_trip_data(net)    
    
    N = list(nodes)
    A1 = list(links)
    A2 = list(links2)
    A = A1+A2
    O = list(orig)
    D = list(dest)    
    
    #---create node-destination demand matrix
    TD = 0
    d = {(i,s):0 for i in N for s in D}
    for r in O:
        for s in D:
            d[r,s] = OD[r,s]
            TD += d[r,s]    
    for s in D:
        d[s,s] = - sum(d[j,s] for j in O)    
    
    #---create link delay function parameters from TNTP data   
    #---link delay functional form: t[i,j] = T[i,j] + c[i,j]*(x**exp[i,j])
    T = {(i,j):fftt[i,j] for (i,j) in A} 
    c = {(i,j):fftt[i,j]*alpha[i,j]/(capa[i,j]**beta[i,j]) for (i,j) in A}    
    e = {(i,j):beta[i,j] for (i,j) in A}         
    
    #---create link cost matrix and budget
    g = {(i,j):cost[i,j] for (i,j) in A2} 
    TC = sum(g[i,j] for (i,j) in A2)
    B = B_prop*TC
    
    #---link delay function linear approximation using m uniform segments 
    #---maximum link flow is instance-specific: value is calibrated for Sioux Falls under base demand
    Mflow = 1e5*scal_flow 
    V = set([i for i in range(0,m+1)])          
    a = {(i,j,v):float() for (i,j) in A for v in V}    
    for (i,j) in A:
        cnt = 0
        step = Mflow/(len(V)-1)
        for v in V:
            a[i,j,v] = cnt*step
            cnt += 1  
    
    #---time and flow scaling        
    #---big-M coefficient is instance-specific: value is calibrated for Sioux Falls under base demand
    Mtt = 1e3*scal_time
    for (i,j) in A:
        T[i,j] = T[i,j]*scal_time
        c[i,j] = c[i,j]*scal_time/(scal_flow**e[i,j])
    for i in N:
        for s in D:
            d[i,s] = d[i,s]*scal_flow     
    
    print('Instance',)
    print('Total scaled demand',TD*scal_flow)
    print('Total cost',TC,'Budget',B)
    
    data = {'nodes':N,'links1':A1,'links2':A2,'links':A,'orig':O,'dest':D,'fftt':T,'coef':c,'exp':e,
            'approx':V,'alpha':a,'cost':g,'demand':d,'budget':B,'Mflow':Mflow,'Mtt':Mtt,'timelimit':timelimit}
    return data


# In[ ]:

#---read network and instance data
#---nodes and links are sets; links is a set of tuples
#---cap, fftt, alpha, beta and cost are dictionaries where the keys are links

def read_network_data(net,NDP):
    network_data = open(net+'_net_'+NDP+'.txt','r')
    lines_net = network_data.readlines()
    network_data.close()
    nb_nodes = int(lines_net[1].split("\t")[0].split(" ")[3])
    nb_links = int(lines_net[3].split("\t")[0].split(" ")[3])    
    nb_links2 = int(lines_net[4].split("\t")[0].split(" ")[4])
    cap = {};fftt = {};alpha = {};beta = {};cost = {};
    nodes = set();links = set();links2 = set();
    
    #---offset is network-specific
    if net == 'SiouxFalls':
        offset = 9
    for i in range(offset,offset+nb_links):
        a = int(lines_net[i].split("\t")[1])
        b = int(lines_net[i].split("\t")[2])
        cap[(a,b)] = float(lines_net[i].split("\t")[3])
        fftt[(a,b)] = float(lines_net[i].split("\t")[5])
        alpha[(a,b)] = float(lines_net[i].split("\t")[6])
        beta[(a,b)] = float(lines_net[i].split("\t")[7])
        nodes.add(a)
        nodes.add(b)
        links.add((a,b))
    for i in range(offset+nb_links,offset+nb_links+nb_links2):
        a = int(lines_net[i].split("\t")[1])
        b = int(lines_net[i].split("\t")[2])
        cap[(a,b)] = float(lines_net[i].split("\t")[3])
        fftt[(a,b)] = float(lines_net[i].split("\t")[5])
        alpha[(a,b)] = float(lines_net[i].split("\t")[6])
        beta[(a,b)] = float(lines_net[i].split("\t")[7])
        cost[(a,b)] = float(lines_net[i].split("\t")[11])
        nodes.add(a)
        nodes.add(b)
        links2.add((a,b))        
    return nodes,links,cap,fftt,alpha,beta,links2,cost


#---read trip data
#---OD_demand is a dictionary where the keys are OD pairs and the values are the demands

def read_trip_data(net):
    trip_data = open(net+'_trips.txt','r')
    trip_lines = trip_data.readlines()
    trip_data.close()
    nb_zones = int(trip_lines[0].split("\t")[0].split(" ")[3])
    total_flow = float(trip_lines[1].split("\t")[0].split(" ")[3])
    dest = 0
    line_nb = 0    
    OD_demand = {} 
    O = set()
    D = set()
    for line in trip_lines:
        if line_nb>=5:
            if line.split(" ")[0]=="Origin":                
                orig = int(line.split(" ")[1])
                orig_flow = 0
            elif len(line.split())>0:
                k=0
                bouh = int(len(line.split())/3)            
                for i in range(bouh):                
                    dest = int(line.split()[k])
                    demand = float(line.split()[k+2].split(";")[0])
                    orig_flow += demand
                    k += 3
                    OD_demand[(orig,dest)] = demand
                    O.add(orig)
                    D.add(dest)
        line_nb += 1
    return OD_demand,O,D


# In[ ]:

#---solve TAP as a convex NLP using Pyomo and IPOPT for a given y-vector

def TAP_cvx(data,yopt):
    time0 = time.time()
    N = data['nodes']
    A = data['links']
    A1 = data['links1']
    A2 = data['links2']
    D = data['dest']
    T = data['fftt']
    c = data['coef']
    e = data['exp']
    d = data['demand']
    Mflow = data['Mflow']
    Mtt = data['Mtt']

    TAP = pyo.ConcreteModel()    
    TAP.x = pyo.Var([(i,j) for (i,j) in A],domain=pyo.NonNegativeReals)
    TAP.xc = pyo.Var([(i,j) for (i,j) in A],[s for s in D],domain=pyo.NonNegativeReals)
    
    TAP.cons = pyo.ConstraintList()
    for i in N:
        for s in D:
            TAP.cons.add(sum(TAP.xc[i,j,s] for j in N if (i,j) in A) - sum(TAP.xc[j,i,s] for j in N if (j,i) in A) == d[i,s])

    TAP.flow = pyo.ConstraintList()
    for (i,j) in A:
        TAP.flow.add(sum(TAP.xc[i,j,s] for s in D) == TAP.x[i,j])
        
    TAP.dsgn = pyo.ConstraintList()
    for (i,j) in A2:
        TAP.dsgn.add(TAP.x[i,j] <= yopt[i,j]*Mflow)        
    
    TAP.obj = pyo.Objective(expr=(sum(TAP.x[i,j]*T[i,j] + c[i,j]/(e[i,j]+1)*((TAP.x[i,j])**(e[i,j]+1)) for (i,j) in A)))
    
    opt = SolverFactory("ipopt.exe")
    results = opt.solve(TAP)
    time_cvx = time.time() - time0

    if (results.solver.status == SolverStatus.ok) and (results.solver.termination_condition == TerminationCondition.optimal):
        x_cvx = {}        
        t_cvx = {}
        TSTT_cvx = 0
        for (i,j) in A:            
            x_cvx[i,j] = pyo.value(TAP.x[i,j])
            if (i,j) in A2 and yopt[i,j] < 1e-4:
                t_cvx[i,j] = Mtt
            else:
                t_cvx[i,j] = T[i,j] + c[i,j]*(x_cvx[i,j]**(e[i,j]))            
            TSTT_cvx += x_cvx[i,j]*T[i,j] + c[i,j]*(x_cvx[i,j]**(e[i,j]+1))
        print('TSTT_cvx',TSTT_cvx)
        return TSTT_cvx,time_cvx,x_cvx       

    elif (results.solver.termination_condition == TerminationCondition.infeasible):
        print('infeasible',results.solver.status,results.solver.termination_condition)  
        return -1,time_cvx,{}
    else:
        print('Solver Status',results.solver.status,results.solver.termination_condition)    
        return -1,time_cvx,{}


# In[ ]:

#---solve MKKT for instance data

def model_MKKT(data):
    MKKT = Model(name='MKKT',log_output=False)
    N = data['nodes']
    A = data['links']
    A1 = data['links1']
    A2 = data['links2']
    D = data['dest']
    V = data['approx']
    T = data['fftt']
    c = data['coef']
    e = data['exp']
    a = data['alpha']
    g = data['cost']
    B = data['budget']    
    d = data['demand']
    Mflow = data['Mflow']
    Mtt = data['Mtt']
    timelimit = data['timelimit']
    Vp = V.difference({0})
    
    #---primal follower variables
    xc = {(i,j,s): MKKT.continuous_var() for (i,j) in A for s in D}
    ll = {(i,j,v): MKKT.continuous_var() for (i,j) in A for v in V}
    lr = {(i,j,v): MKKT.continuous_var() for (i,j) in A for v in V}
    
    #---dual follower variables
    pi = {(i,s): MKKT.continuous_var() for i in N for s in D}
    beta = {(i,j): MKKT.continuous_var() for (i,j) in A}
    gamma = {(i,j): MKKT.continuous_var() for (i,j) in A}
    mu = {(i,j): MKKT.continuous_var() for (i,j) in A2}
    
    #---leader and linearization variables
    y = {(i,j): MKKT.binary_var() for (i,j) in A2}
    phi = {(i,j): MKKT.continuous_var() for (i,j) in A2}
    
    MKKT.add_constraint(sum(y[i,j]*g[i,j] for (i,j) in A2) <= B)    
    for i in N:
        for s in D:
            MKKT.add_constraint(sum(xc[i,j,s] for j in N if (i,j) in A) 
                               - sum(xc[j,i,s] for j in N if (j,i) in A) == d[i,s])
    for (i,j) in A:
        MKKT.add_constraint(sum(xc[i,j,s] for s in D) == sum(ll[i,j,v]*a[i,j,v-1] + lr[i,j,v]*a[i,j,v] for v in Vp))        
        MKKT.add_constraint(sum(ll[i,j,v] + lr[i,j,v] for v in V) == 1)         
    
    for (i,j) in A2:
        for s in D:
            MKKT.add_constraint(xc[i,j,s] <= y[i,j]*Mflow)
        
    for (i,j) in A1:
        for s in D:
            MKKT.add_constraint(pi[i,s] - pi[j,s] + beta[i,j] >= - T[i,j])
    for (i,j) in A2:
        for s in D:
            MKKT.add_constraint(pi[i,s] - pi[j,s] + beta[i,j] + mu[i,j] >= - T[i,j])
    for (i,j) in A:        
        for v in Vp:    
            MKKT.add_constraint(- beta[i,j]*a[i,j,v] + gamma[i,j] >= - (c[i,j]/(e[i,j]+1))*(a[i,j,v]**(e[i,j]+1)))
        MKKT.add_constraint(gamma[i,j] >= 0)
    
    for (i,j) in A2:
        MKKT.add_constraint(phi[i,j] <= mu[i,j]) 
        MKKT.add_constraint(phi[i,j] >= mu[i,j] - (1 - y[i,j])*Mtt)
        MKKT.add_constraint(phi[i,j] <= y[i,j]*Mtt)
        MKKT.add_constraint(phi[i,j] >= 0)

    # primal-dual constraint
    MKKT.add_constraint(sum(T[i,j]*sum(xc[i,j,s] for s in D) 
                           + (c[i,j]/(e[i,j]+1))*sum(ll[i,j,v]*(a[i,j,v-1]**(e[i,j]+1)) + lr[i,j,v]*(a[i,j,v]**(e[i,j]+1)) 
                                                     for v in Vp) for (i,j) in A) 
                       <= -(sum(pi[i,s]*d[i,s] for i in N for s in D) 
                           + sum(gamma[i,j] for (i,j) in A) + sum(phi[i,j]*Mflow for (i,j) in A2)))
        
    for (i,j) in A:        
        for s in D:
            MKKT.add_constraint(xc[i,j,s] >= 0)
        for v in V:
            MKKT.add_constraint(ll[i,j,v] >= 0)
            MKKT.add_constraint(lr[i,j,v] >= 0)            
    for (i,j) in A2:
        MKKT.add_constraint(mu[i,j] >= 0)                                
    
    MKKT.minimize(sum(T[i,j]*sum(xc[i,j,s] for s in D) for (i,j) in A)
                 + sum(c[i,j]*sum(ll[i,j,v]*(a[i,j,v-1]**(e[i,j]+1)) + lr[i,j,v]*(a[i,j,v]**(e[i,j]+1)) 
                                  for v in Vp) for (i,j) in A))
       
    MKKT.parameters.threads = 1
    MKKT.parameters.mip.display = 0
    MKKT.parameters.timelimit = timelimit
    sol = MKKT.solve()
    time_MKKT = MKKT.solve_details.time
    if MKKT.solve_details.status == 'integer infeasible':
        print('status\t%s' % MKKT.solve_details.status)
        return -1,-1,time_MKKT,{},{}
    
    print('\n---MKKT-------------------------------------')
    print('status\t%s' % MKKT.solve_details.status)
    print('time\t%.2f' % MKKT.solve_details.time)
    print('OPT\t%.3f' % MKKT.objective_value)
    UB_MKKT = MKKT.objective_value
    gap_MKKT = MKKT.solve_details.mip_relative_gap    
    yopt = {(i,j):MKKT.solution.get_value(y[i,j]) for (i,j) in A2}
    xcopt = {(i,j,s):MKKT.solution.get_value(xc[i,j,s]) for (i,j) in A for s in D}
    xopt = {(i,j):sum(xcopt[i,j,s] for s in D) for (i,j) in A}
    return UB_MKKT,gap_MKKT,time_MKKT,yopt,xopt


# In[ ]:

#---solve SO-relaxation with iterative Interdiction Cuts

def algo_SOIC(data):
    SOIC = Model(name='SOIC',log_output=False)
    N = data['nodes']
    A = data['links']
    A1 = data['links1']
    A2 = data['links2']
    D = data['dest']
    V = data['approx']
    T = data['fftt']
    c = data['coef']
    e = data['exp']
    a = data['alpha']
    g = data['cost']
    B = data['budget']    
    d = data['demand']
    Mflow = data['Mflow']
    timelimit = data['timelimit']
    Vp = V.difference({0})
 
    xc = {(i,j,s): SOIC.continuous_var() for (i,j) in A for s in D}
    y = {(i,j): SOIC.binary_var() for (i,j) in A2}
    ll = {(i,j,v): SOIC.continuous_var() for (i,j) in A for v in V}
    lr = {(i,j,v): SOIC.continuous_var() for (i,j) in A for v in V}
    
    SOIC.add_constraint(sum(y[i,j]*g[i,j] for (i,j) in A2) <= B)
    for i in N:
        for s in D:
            SOIC.add_constraint(sum(xc[i,j,s] for j in N if (i,j) in A) 
                               - sum(xc[j,i,s] for j in N if (j,i) in A) == d[i,s])
    for (i,j) in A2:        
        for s in D:            
            SOIC.add_constraint(xc[i,j,s] <= y[i,j]*Mflow)

    for (i,j) in A:
        SOIC.add_constraint(sum(xc[i,j,s] for s in D) == sum(ll[i,j,v]*a[i,j,v-1] + lr[i,j,v]*a[i,j,v] for v in Vp))      
        SOIC.add_constraint(sum(ll[i,j,v] + lr[i,j,v] for v in V) == 1) 
      
    for (i,j) in A:
        for s in D:
            SOIC.add_constraint(xc[i,j,s] >= 0)
        for v in V:
            SOIC.add_constraint(ll[i,j,v] >= 0)
            SOIC.add_constraint(lr[i,j,v] >= 0)
        
    SOIC.minimize(sum(T[i,j]*sum(xc[i,j,s] for s in D) 
                      + c[i,j]*sum(ll[i,j,v]*(a[i,j,v-1]**(e[i,j]+1)) 
                                   + lr[i,j,v]*(a[i,j,v]**(e[i,j]+1)) for v in Vp) for (i,j) in A))
       
    SOIC.parameters.threads = 1
    SOIC.parameters.mip.display = 0
    SOIC.parameters.timelimit = timelimit
    
    converged = False
    nit = 0
    print('\n---SOIC-------------------------------------')
    t0 = time.time()
    UB_SOIC = 1e9
    yUB = {}
    xUB = {}
    while converged == False:  
     
        #---solve SO-relaxation
        sol = SOIC.solve()        
        if SOIC.solve_details.status == 'integer infeasible':
            print('status\t%s' % SOIC.solve_details.status)
            break
 
        SO_TSTT = SOIC.objective_value
        yopt = {(i,j):SOIC.solution.get_value(y[i,j]) for (i,j) in A2}
        xcopt = {(i,j,s):SOIC.solution.get_value(xc[i,j,s]) for (i,j) in A for s in D}
               
        #---solve TAP with yopt
        UE_TSTT,xUE,UE_time = TAP_lin(data,yopt)      
        print('%d\t%.3f\t%.3f\t%.3f\t%.3f' % (nit,SO_TSTT,UE_TSTT,SOIC.solve_details.time,UE_time))
        
        #---update UB if new TAP solution better than incumbent
        if UE_TSTT < UB_SOIC:
            print('update UB')
            UB_SOIC = UE_TSTT
            yUB = yopt
            xUB = xUE
        
        #---check convergence
        if SO_TSTT >= UB_SOIC:
            converged = True
            break        
        if (time.time() - t0) >= timelimit:
            print('>> time limit')
            break
        
        #---if not converged, add interdiction cut
        SOIC.add_constraint(sum((1-yopt[i,j])*y[i,j] + yopt[i,j]*(1 - y[i,j]) for (i,j) in A2) >= 1)
        
        nit += 1
        
    time_SOIC = min(time.time() - t0, timelimit)
    if UB_SOIC > 0:
        if SO_TSTT >= UB_SOIC:
            gap_SOIC = 0.0
        else:
            gap_SOIC = (UB_SOIC - SO_TSTT)/UB_SOIC
    else:
        gap_SOIC = 1e9
    print('UB\t%.3f' % (UB_SOIC))
    print('GAP\t%.3f' % (gap_SOIC))
    print('TIME\t%.2f' % (time_SOIC))
    return UB_SOIC,gap_SOIC,time_SOIC,yUB,xUB


#---solve TAP with y-vector

def TAP_lin(data,yopt):
    TAP = Model(name='TAP',log_output=False)
    N = data['nodes']
    A = data['links']
    A1 = data['links1']
    A2 = data['links2']
    D = data['dest']
    V = data['approx']
    T = data['fftt']
    c = data['coef']
    e = data['exp']
    a = data['alpha']
    d = data['demand']    
    Mflow = data['Mflow']
    timelimit = data['timelimit']
    Vp = V.difference({0})
    
    xc = {(i,j,s): TAP.continuous_var() for (i,j) in A for s in D}    
    ll = {(i,j,v): TAP.continuous_var() for (i,j) in A for v in V}
    lr = {(i,j,v): TAP.continuous_var() for (i,j) in A for v in V}
    
    for i in N:
        for s in D:
            TAP.add_constraint(sum(xc[i,j,s] for j in N if (i,j) in A) 
                               - sum(xc[j,i,s] for j in N if (j,i) in A) == d[i,s])
    for (i,j) in A:
        TAP.add_constraint(sum(xc[i,j,s] for s in D) 
                           - sum(ll[i,j,v]*a[i,j,v-1] + lr[i,j,v]*a[i,j,v] for v in Vp) == 0) 
        TAP.add_constraint(sum(ll[i,j,v] + lr[i,j,v] for v in V) == 1)
        
    for (i,j) in A2:        
        for s in D:            
            TAP.add_constraint(xc[i,j,s] <= yopt[i,j]*Mflow)
                
    TAP.minimize(sum(T[i,j]*sum(xc[i,j,s] for s in D) 
                        + (c[i,j]/(e[i,j]+1))*sum(ll[i,j,v]*(a[i,j,v-1]**(e[i,j]+1)) + lr[i,j,v]*(a[i,j,v]**(e[i,j]+1)) 
                                                  for v in Vp) for (i,j) in A))
        
    for (i,j) in A:        
        for s in D:
            TAP.add_constraint(xc[i,j,s] >= 0)
        for v in V:
            TAP.add_constraint(ll[i,j,v] >= 0)
            TAP.add_constraint(lr[i,j,v] >= 0)            

    TAP.parameters.threads = 1
    TAP.parameters.simplex.display = 0
    TAP.parameters.timelimit = timelimit
    sol = TAP.solve()
    
    UE_time = TAP.solve_details.time
    xcopt = {(i,j,s):TAP.solution.get_value(xc[i,j,s]) for (i,j) in A for s in D}
    llopt = {(i,j,v):TAP.solution.get_value(ll[i,j,v]) for (i,j) in A for v in V}
    lropt = {(i,j,v):TAP.solution.get_value(lr[i,j,v]) for (i,j) in A for v in V}
    xUE = {(i,j):sum(xcopt[i,j,s] for s in D) for (i,j) in A}
    UE_TSTT = 0
    topt = {(i,j):float()}
    for (i,j) in A:
        topt[i,j] = (T[i,j] + c[i,j]*sum(llopt[i,j,v]*(a[i,j,v-1]**e[i,j]) + lropt[i,j,v]*(a[i,j,v]**e[i,j]) for v in Vp))        
        UE_TSTT += topt[i,j]*sum(xcopt[i,j,s] for s in D)  
    return UE_TSTT,xUE,UE_time


# In[ ]:

#---solve SO-relaxation with Branch and Bound

def algo_SOBB(data):
    converged = False
    nit = 0
    node = 0
    gap = 1.0*float('inf')
    UB = 1e9
    best_LB = 0
    y_SOBB = {}
    x_SOBB = {}
    tree = {0:{'parent':0,'children':[],'solved':False,'active':True,'LB':None,'UB':None,
               'fixed0':[],'fixed1':[],'yopt':{},'xUB':{},'xtopt':{}}}
    candidates = [0]
    t0 = time.time()
    timelimit = data['timelimit']
    print('\n---SOBB---------------------------------')
    while converged == False:  
    
        #---select subproblem SP to explore
        SP,can = node_selection(candidates,tree)     
        
        #---if SP is not solved, solve SP
        if SP['solved'] == False:
            SP = model_SOBB(data,SP)
            SP = TAP_lin_BB(data,SP)  
            SP['solved'] = True
            print('%d\t%.3f\t%.3f\t%.3f\t%.3f' % (nit,SP['LB'],SP['UB'],SP['LB_time'],SP['UB_time']))
        else:
            print('%d\t%.3f\t%.3f\t---\t---' % (nit,SP['LB'],SP['UB']))
            
        #---if new UB is better than incumbent, update UB and solution
        if SP['UB'] < UB:
            print('update UB')
            UB = SP['UB']
            y_SOBB = SP['yopt']
            x_SOBB = SP['xUB']
            
            #---prune tree with new UB
            for i in tree:
                if tree[i]['active'] == True and tree[i]['LB'] >= UB:
                    tree[i]['active'] = False       
                      
        #---if SP['LB'] < UB, branch on SP else prune tree
        if SP['LB'] < UB:
            tree = branch(SP,can,tree,data)
        
        #---update candidate list
        SP['active'] = False
        candidates = [i for i in tree if tree[i]['active']==True]
        if len(candidates) == 0:
            converged = True
            best_LB = UB
            gap = 0
        else:
            best_LB = min([tree[i]['LB'] for i in candidates])
            if (UB - best_LB)/UB < gap:
                gap = (UB - best_LB)/UB
                    
        if (time.time() - t0) >= timelimit:
            print('>> time limit')
            break

        nit += 1

    time_SOBB = min(time.time() - t0,timelimit)
    UB_SOBB = UB
    gap_SOBB = gap
    print('time\t%.2f' % (time_SOBB))
    print('OPT\t%.3f' % (UB_SOBB))
    print('GAP\t%.3f' % (gap_SOBB))
    return UB_SOBB,gap_SOBB,time_SOBB,y_SOBB,x_SOBB


#---node selection rule: best bound first search

def node_selection(candidates,tree):
    min_LB = min([tree[tree[i]['parent']]['LB'] for i in candidates])
    candidates = [i for i in candidates if tree[tree[i]['parent']]['LB']==min_LB]
    return tree[candidates[0]],candidates[0]


#---branching rule: branch on highest x[i,j]*t[i,j] value

def branch(SP,can,tree,data):
    fixed = SP['fixed0'] + SP['fixed1']  
    A2 = data['links2']      
    if len(fixed) == len(A2):
        return tree
    else:
        cnt = len(tree)-1      
        free = [(i,j) for (i,j) in A2 if (i,j) not in fixed]
        maxxt = 0
        for (i,j) in free:
            if SP['xtopt'][i,j] >= maxxt - 1e-3:
                maxxt = SP['xtopt'][i,j]
                ybr = (i,j)
                         
        fixed00 = copy.deepcopy(SP['fixed0'])
        fixed00.append(ybr)
        fixed01 = copy.deepcopy(SP['fixed1'])
        fixed10 = copy.deepcopy(SP['fixed0'])
        fixed11 = copy.deepcopy(SP['fixed1'])
        fixed11.append(ybr)
        
        if SP['yopt'][ybr] == 1:
            SP['children'].append(cnt+1)
            tree[cnt+1] = {'parent':can,'children':[],'solved':False,'active':True,'LB':SP['LB'],'UB':None,
                           'fixed0':fixed00,'fixed1':fixed01,'yopt':{},'xUB':{},'xtopt':{}} 
            SP['children'].append(cnt+2)
            tree[cnt+2] = {'parent':can,'children':[],'solved':True,'active':True,'LB':SP['LB'],'UB':SP['UB'],
                           'fixed0':fixed10,'fixed1':fixed11,'yopt':SP['yopt'],'xUB':SP['xUB'],'xtopt':SP['xtopt']}   
        else:
            SP['children'].append(cnt+1)
            tree[cnt+1] = {'parent':can,'children':[],'solved':True,'active':True,'LB':SP['LB'],'UB':SP['UB'],
                           'fixed0':fixed00,'fixed1':fixed01,'yopt':SP['yopt'],'xUB':SP['xUB'],'xtopt':SP['xtopt']} 
            
            #---if adding ybr does not violate budget, add child
            if data['cost'][ybr] + sum(data['cost'][i,j] for (i,j) in fixed) <= data['budget']:            
                SP['children'].append(cnt+2)
                tree[cnt+2] = {'parent':can,'children':[],'solved':False,'active':True,'LB':SP['LB'],'UB':None,
                               'fixed0':fixed10,'fixed1':fixed11,'yopt':{},'xUB':{},'xtopt':{}}
        return tree
    
    
#---solve SO-relaxation at subproblem SP    
 
def model_SOBB(data,SP):
    SOBB = Model(name='SOBB',log_output=False)
    N = data['nodes']
    A = data['links']
    A1 = data['links1']
    A2 = data['links2']
    D = data['dest']
    V = data['approx']
    T = data['fftt']
    c = data['coef']
    e = data['exp']
    a = data['alpha']
    g = data['cost']
    B = data['budget']    
    d = data['demand']
    Mflow = data['Mflow']
    timelimit = data['timelimit']
    Vp = V.difference({0})
 
    xc = {(i,j,s): SOBB.continuous_var() for (i,j) in A for s in D}
    y = {(i,j): SOBB.binary_var() for (i,j) in A2}
    ll = {(i,j,v): SOBB.continuous_var() for (i,j) in A for v in V}
    lr = {(i,j,v): SOBB.continuous_var() for (i,j) in A for v in V}
    
    for (i,j) in A2:
        if (i,j) in SP['fixed0']:
            SOBB.add_constraint(y[i,j] == 0)
        if (i,j) in SP['fixed1']:
            SOBB.add_constraint(y[i,j] == 1)
    
    SOBB.add_constraint(sum(y[i,j]*g[i,j] for (i,j) in A2) <= B)
    for i in N:
        for s in D:
            SOBB.add_constraint(sum(xc[i,j,s] for j in N if (i,j) in A) 
                                   - sum(xc[j,i,s] for j in N if (j,i) in A) == d[i,s])
    for (i,j) in A2:        
        for s in D:            
            SOBB.add_constraint(xc[i,j,s] <= y[i,j]*Mflow)

    for (i,j) in A:
        SOBB.add_constraint(sum(xc[i,j,s] for s in D) == sum(ll[i,j,v]*a[i,j,v-1] + lr[i,j,v]*a[i,j,v] for v in Vp))      
        SOBB.add_constraint(sum(ll[i,j,v] + lr[i,j,v] for v in V) == 1) 
      
    for (i,j) in A:
        for s in D:
            SOBB.add_constraint(xc[i,j,s] >= 0)
        for v in V:
            SOBB.add_constraint(ll[i,j,v] >= 0)
            SOBB.add_constraint(lr[i,j,v] >= 0)
        
    SOBB.minimize(sum(T[i,j]*sum(xc[i,j,s] for s in D) 
                      + c[i,j]*sum(ll[i,j,v]*(a[i,j,v-1]**(e[i,j]+1)) 
                                   + lr[i,j,v]*(a[i,j,v]**(e[i,j]+1)) for v in Vp) for (i,j) in A))
       
    SOBB.parameters.threads = 1
    SOBB.parameters.mip.display = 0
    SOBB.parameters.timelimit = timelimit
  
    sol = SOBB.solve() 
    SP['LB_time'] = SOBB.solve_details.time
    
    if SOBB.solve_details.status == 'integer infeasible':
        print('status\t%s' % SOBB.solve_details.status)
        SP['LB'] = -1
        SP['yopt'] = {}
        SP['xtopt'] = {}
    else:
        xcopt = {(i,j,s):SOBB.solution.get_value(xc[i,j,s]) for (i,j) in A for s in D}        
        llopt = {(i,j,v):SOBB.solution.get_value(ll[i,j,v]) for (i,j) in A for v in V}
        lropt = {(i,j,v):SOBB.solution.get_value(lr[i,j,v]) for (i,j) in A for v in V}
        xopt = {(i,j):sum(xcopt[i,j,s] for s in D) for (i,j) in A}           
        xtopt = {}
        for (i,j) in A:
            xtopt[i,j] = sum(T[i,j]*sum(xcopt[i,j,s] for s in D) 
                             + c[i,j]*sum(llopt[i,j,v]*(a[i,j,v-1]**(e[i,j]+1)) 
                                          + lropt[i,j,v]*(a[i,j,v]**(e[i,j]+1)) for v in Vp) for (i,j) in A)
        
        SP['LB'] = SOBB.objective_value
        SP['yopt'] = {(i,j):SOBB.solution.get_value(y[i,j]) for (i,j) in A2}
        SP['xtopt'] = xtopt
    return SP


#---solve TAP at subproblem SP using y-vector from SO-relaxation

def TAP_lin_BB(data,SP):
    TAP = Model(name='TAP',log_output=False)
    N = data['nodes']
    A = data['links']
    A1 = data['links1']
    A2 = data['links2']
    D = data['dest']
    V = data['approx']
    T = data['fftt']
    c = data['coef']
    e = data['exp']
    a = data['alpha']
    d = data['demand']    
    Mflow = data['Mflow']
    timelimit = data['timelimit']
    Vp = V.difference({0})
    
    yopt = SP['yopt']
    
    xc = {(i,j,s): TAP.continuous_var() for (i,j) in A for s in D}    
    ll = {(i,j,v): TAP.continuous_var() for (i,j) in A for v in V}
    lr = {(i,j,v): TAP.continuous_var() for (i,j) in A for v in V}
    
    for i in N:
        for s in D:
            TAP.add_constraint(sum(xc[i,j,s] for j in N if (i,j) in A) 
                               - sum(xc[j,i,s] for j in N if (j,i) in A) == d[i,s])
    for (i,j) in A:
        TAP.add_constraint(sum(xc[i,j,s] for s in D) 
                           - sum(ll[i,j,v]*a[i,j,v-1] + lr[i,j,v]*a[i,j,v] for v in Vp) == 0) 
        TAP.add_constraint(sum(ll[i,j,v] + lr[i,j,v] for v in V) == 1)
        
    for (i,j) in A2:        
        for s in D:            
            TAP.add_constraint(xc[i,j,s] <= yopt[i,j]*Mflow)
        
                
    TAP.minimize(sum(T[i,j]*sum(xc[i,j,s] for s in D) 
                        + (c[i,j]/(e[i,j]+1))*sum(ll[i,j,v]*(a[i,j,v-1]**(e[i,j]+1)) + lr[i,j,v]*(a[i,j,v]**(e[i,j]+1)) 
                                                  for v in Vp) for (i,j) in A))
        
    for (i,j) in A:        
        for s in D:
            TAP.add_constraint(xc[i,j,s] >= 0)
        for v in V:
            TAP.add_constraint(ll[i,j,v] >= 0)
            TAP.add_constraint(lr[i,j,v] >= 0)            

    TAP.parameters.threads = 1
    TAP.parameters.simplex.display = 0
    TAP.parameters.timelimit = timelimit
    sol = TAP.solve()
  
    xcopt = {(i,j,s):TAP.solution.get_value(xc[i,j,s]) for (i,j) in A for s in D}
    llopt = {(i,j,v):TAP.solution.get_value(ll[i,j,v]) for (i,j) in A for v in V}
    lropt = {(i,j,v):TAP.solution.get_value(lr[i,j,v]) for (i,j) in A for v in V}
    xUE = {(i,j):sum(xcopt[i,j,s] for s in D) for (i,j) in A}
    UE_TSTT = 0
    topt = {(i,j):float()}
    for (i,j) in A:
        topt[i,j] = (T[i,j] + c[i,j]*sum(llopt[i,j,v]*(a[i,j,v-1]**e[i,j]) + lropt[i,j,v]*(a[i,j,v]**e[i,j]) for v in Vp))        
        UE_TSTT += topt[i,j]*sum(xcopt[i,j,s] for s in D) 
    SP['UB'] = UE_TSTT
    SP['xUB'] = xUE
    SP['UB_time'] = TAP.solve_details.time
    return SP


# In[ ]:

#---budget sensitivity experiment

output = open("Benchmark_budget.txt",'w')
net = 'SiouxFalls'
gtime0 = time.time()

#---iterate over 10- and 20-link instance sets
for i in range(1,3):
    
    #---iterate over 10 instances
    for j in range(1,11):        
        
        NDP = 'NDP_'+str(i*10)+'_'+str(j)
        print('instance',NDP)
        
        #---read instance data
        data = read_instance(net,NDP,0,100,1e-0,1e-3,600)
        
        #---get total cost
        TC = sum(data['cost'][i,j] for (i,j) in data['cost'])
        
        #---iterate over 3 levels of budget: 25%, 50% and 75% of TC
        for k in range(1,4):            
            data['budget'] = TC*k/4
            print('>>> NDP',i*10,j,100*k/4)            
            UB_MKKT,gap_MKKT,time_MKKT,y_MKKT,x_MKKT = model_MKKT(data)
            UB_SOIC,gap_SOIC,time_SOIC,y_SOIC,x_SOIC = algo_SOIC(data)
            UB_SOBB,gap_SOBB,time_SOBB,y_SOBB,x_SOBB = algo_SOBB(data)
        
            #---determine true TSTT using convex local solver 
            if UB_MKKT <= UB_SOIC and UB_MKKT <= UB_SOBB: 
                TSTT_cvx,time_cvx,x_cvx = TAP_cvx(data,y_MKKT)
            elif UB_SOIC <= UB_MKKT and UB_SOIC <= UB_SOBB: 
                TSTT_cvx,time_cvx,x_cvx = TAP_cvx(data,y_SOIC)
            elif UB_SOBB <= UB_MKKT and UB_SOBB <= UB_SOIC: 
                TSTT_cvx,time_cvx,x_cvx = TAP_cvx(data,y_SOBB)
            
            output.write('%d\t%d\t%.0f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n' 
                         % (i*10,j,100*k/4,
                            UB_MKKT,gap_MKKT,time_MKKT,
                            UB_SOIC,gap_SOIC,time_SOIC,
                            UB_SOBB,gap_SOBB,time_SOBB,
                            TSTT_cvx))
            output.flush()
print('total time', (time.time() - gtime0)/3600,'hours')
output.close()


# In[ ]:

#---demand sensitivity experiment

output = open("Benchmark_demand.txt",'w')
net = 'SiouxFalls'
gtime0 = time.time()

#---iterate over 10- and 20-link instance sets
for i in range(1,3):
    
    #---iterate over 10 instances
    for j in range(1,11):        
        
        NDP = 'NDP_'+str(i*10)+'_'+str(j)
        print('instance',NDP)
        
        #---read and create instance data with 25% of budget
        data = read_instance(net,NDP,0.25,100,1e-0,1e-3,600)
        
        #---store base demand
        base_demand = copy.deepcopy(data['demand'])
   
        #---iterate over 3 levels of demand: 50%, 100% and 150%
        for k in range(1,4):            
            if k==1:
                factor = 0.5
            if k==2:
                factor = 1.0
            if k==3:
                factor = 1.5                
            print('>>> NDP',i*10,j,factor)
            for rs in data['demand']:                
                data['demand'][rs] = factor*base_demand[rs]   
                
            UB_MKKT,gap_MKKT,time_MKKT,y_MKKT,x_MKKT = model_MKKT(data)
            UB_SOIC,gap_SOIC,time_SOIC,y_SOIC,x_SOIC = algo_SOIC(data)
            UB_SOBB,gap_SOBB,time_SOBB,y_SOBB,x_SOBB = algo_SOBB(data)
        
            #---determine true TSTT using convex local solver 
            if UB_MKKT <= UB_SOIC and UB_MKKT <= UB_SOBB: 
                TSTT_cvx,time_cvx,x_cvx = TAP_cvx(data,y_MKKT)
            elif UB_SOIC <= UB_MKKT and UB_SOIC <= UB_SOBB: 
                TSTT_cvx,time_cvx,x_cvx = TAP_cvx(data,y_SOIC)
            elif UB_SOBB <= UB_MKKT and UB_SOBB <= UB_SOIC: 
                TSTT_cvx,time_cvx,x_cvx = TAP_cvx(data,y_SOBB)
            
            output.write('%d\t%d\t%.0f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n' 
                         % (i*10,j,100*factor,
                            UB_MKKT,gap_MKKT,time_MKKT,
                            UB_SOIC,gap_SOIC,time_SOIC,
                            UB_SOBB,gap_SOBB,time_SOBB,
                            TSTT_cvx))
            output.flush()
print('total time', (time.time() - gtime0)/3600,'hours')
output.close()

