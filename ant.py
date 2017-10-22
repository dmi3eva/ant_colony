import numpy as np

INF = 1E9
ans = INF
s = 1 #starting vertex
f = 5 #finishing vertex


#return adjacency matrix (read from file with description of edjes)
#First line in file: V E
#Others lines: "from to cost"
def gen_adj_matrix():
    f = open('graph.txt', 'r')
    content = f.read().split('\n')
    numbers = []
    ind = 0
    for line in content:
        numbers.append(line.split(' '))
        ind += 1
    V, E = int(numbers[0][0]), int(numbers[0][1])
    adj = []
    for i in range(V + 1):
        adj.append([0] * (V + 1))
    for i in range(1, E + 1):
        
        adj[int(numbers[i][0])][int(numbers[i][1])] = int(numbers[i][2])
        adj[int(numbers[i][1])][int(numbers[i][0])] = int(numbers[i][2])
    return V, E, adj #V - number of vertexes,  E - number of edges
    
#returns distances from vertex number s to vertex number f    
def dijkstra(s, f, V, adj):
    d = [INF] * (V + 1)
    processed = [0] * (V + 1)
    parent = [-1] * (V + 1) 
    d[s] = 0
    while (processed[f] == 0):
        
        min_value = INF
        for i in range(1, V + 1):
            if ((d[i] < min_value) & (processed[i] == 0)):
                min_value = d[i]
                cur = i
        processed[cur] = 1
        
        for i in range(1, V + 1):
            if ((processed[i] == 0) & (adj[cur][i] > 0)):
                if (d[i] > d[cur] + adj[cur][i]):
                    d[i] = d[cur] + adj[cur][i]
                    parent[i] = cur
        
    return d, parent
     
def dijkstra_path(s, f, parent): 
    cur = f
    path = []
    while (cur != s):
        path.append(cur)
        cur = parent[cur]
    path.append(s)  
    path.reverse()
    return path
    
def degree (vertex, adj, V):
    deg = 0
    for i in range(1, V + 1):
        if (adj[vertex][i] != 0):
            deg += 1
    return deg

#initialize starting parameters of the algorithm
def init_parameters(V, E, adj, alpha_v, ro_v):
    alpha = alpha_v #greed
    beta = 1 - alpha #herdness
    ro = ro_v #rate of evaporation of pheromone
    p = []
    nu = []
    tau = []
    for i in range(V + 1):
        p.append([0.0] * (V + 1))
        nu.append([0.0] * (V + 1))
        tau.append([0.0] * (V + 1))
    for i in range(1, V + 1):
        for j in range(1, V + 1):
            if (adj[i][j] != 0):
                nu[i][j] = 1.0 / adj[i][j]
                tau[i][j] = 1.0 / degree(i, adj, V)
    p = make_prob(p, adj, alpha, beta, tau, nu)
    return alpha, beta, ro, tau, nu, p
def make_prob(p, adj, alpha, beta, tau, nu):
    for i in range(1, V + 1):
        for j in range(1, V + 1):
            if (adj[i][j] != 0):
                denominator = 0
                for l in range(1, V + 1):
                    if (adj[i][l] != 0):
                        denominator += np.power(tau[i][l], alpha) * np.power(nu[i][l], beta)
                p[i][j] = np.power(tau[i][j], alpha) * np.power(nu[i][j], beta) / denominator
    return p    
    
class Ant:
    l = 0
    pos = -1
    path = []
    cost = []

#generate ants with zero ants mileage  
def generate_ants(number, s):
    ants = []
    for k in range(number):
        tmp = Ant() 
        tmp.pos = s
        ants.append(tmp)
    return ants
   
#returns step of particular ant   
def step(ans, ant, p, V, adj, alpha, beta, ro, tau, nu):
    done = 0
    nei = adj_list(ant.pos, V, adj)
    
    summ = 0
    for v in nei:
        summ += p[ant.pos][v]
    
    r = np.random.sample() * summ
    
    summ = 0
    for v in nei:
        summ += p[ant.pos][v]
        if (r <= summ):
            ant.l += adj[ant.pos][v]
            ant.pos = v
            ant.path.append(v)
            ant.cost.append(ant.l)
            done = 1
            break
    
    if (ant.pos == s):
        ant.l = 0
        ant.path = []
    if (ant.pos == f):
        if (ans > ant.l):
            ans = ant.l
        p, tau = recalculate_pheromone(ant, p, V, adj, alpha, beta, ro, tau, nu)
        ant.pos = s
        ant.l = 0
        ant.path = []
    return ans, p, tau

def recalculate_pheromone(ant, p, V, adj, alpha, beta, ro, tau, nu):
    for i in range(0, len(ant.path) - 1):
        fromv = ant.path[i]
        tov = ant.path[i + 1] 
        if (adj[fromv][tov] != 0):
            tau[fromv][tov] = (1 - ro) * tau[fromv][tov] + 1.0 / (adj[fromv][tov])
    p = make_prob(p, adj, alpha, beta, tau, nu)
    return p, tau
#returns list of all adjacent vartexes
def adj_list(vertex, V, adj):
    lst = []
    for i in range(1, V + 1):
        if (adj[vertex][i] != 0):
            lst.append(i)
    return lst
   

V, E, adj = gen_adj_matrix()
d, parent = dijkstra(s, f, V, adj)
path = dijkstra_path(s, f, parent)
ants = generate_ants(4, s)
print ("Path:") 
print(path)



print ("Answer is " + str(d[f]))    

experiments = 100
correct = 0
res = 0
for ex in range(1, experiments):
    alpha, beta, ro, tau, nu, p = init_parameters(V, E, adj, 0.9, 0.5)
    ans = INF
    max_step = 100
    step_num = 0
    i = 0
    while (ans != d[f] & i < max_step):
        i += 1
        if (d[f] == ans):
            break
        for j in range(0, len(ants)):
            ans_n, p, tau = step(ans, ants[j], p, V, adj, alpha, beta, ro, tau, nu)
            if (ans_n < ans):
                #print(str(ans_n) + " " + str(i))
                ans = ans_n
                step_num = i
                if (d[f] == ans):
                    
                    break
    if (d[f] == ans):    
        res += step_num
        correct += 1
    #print(str(step_num) + " " + str(ans))
print(res * 1.0 / correct)