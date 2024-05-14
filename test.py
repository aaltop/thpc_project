import numpy as np
import matplotlib.pyplot as plt

def new_route(num_of_locations, rng):

    route = np.arange(num_of_locations)[1:]
    return np.concatenate(([0], rng.permutation(route)))

def calculate_total_distance(route, distances):

    n = route.size
    total_distance = 0.0
    for i in range(n):

        from_index = route[i]
        to_index = route[(i+1)%n]

        total_distance += distances[from_index, to_index]

    return total_distance

def breed(route1, route2, distances, rng):

    child = np.zeros_like(route1)-1
    not_added = np.arange(len(route1))
    child[0] = 0
    parents = (route1, route2)
    for i in range(len(route1)-1):
        i += 1

        # take away the latest addition
        not_added = not_added[not_added != child[i-1]]
        if len(not_added) < 2:
            child[i] = not_added[0]
            return child

        prev_loc = child[i-1]
        # if false, pick parent1, else parent2
        from_parent =  int(distances[prev_loc, route1[i]] > distances[prev_loc, route2[i]])
        new_candidate = parents[from_parent][i]
        if not new_candidate in child:
            child[i] = new_candidate
            continue

        new_candidate = parents[not from_parent][i]
        if not new_candidate in child:
            child[i] = new_candidate
            continue

        child[i] = rng.choice(not_added)

    # NOTE: return is above

def mutate(route, rng):

    i1, i2 = rng.choice(np.arange(len(route)-1)+1, size=2, replace=False)
    route[i1], route[i2] = route[i2], route[i1]

def calculate_fitness(candidate_distances):

    one_minus_max = 1 - candidate_distances/candidate_distances.max()
    return one_minus_max/one_minus_max.sum()

def cull(fitness, num_survivors, rng):
    '''
    Return indices of survivors based on their <fitness>.
    '''

    return rng.choice(len(fitness), size=num_survivors, p=fitness, replace=False)

def create_candidates(num_candidates, num_locs, distances, rng):

    routes = tuple(new_route(num_locs, rng) for _ in range(num_candidates+1))
    routes_distances = np.array([calculate_total_distance(route,distances) for route in routes])
    fitness = calculate_fitness(routes_distances)
    is_fit = fitness > 0
    routes = tuple(route for i, route in enumerate(routes) if is_fit[i])
    fitness = fitness[is_fit]
    return routes, fitness

def find_optimal_route(distances, num_routes, num_bred, mutation_chance, generations, rng):

    if num_bred > (num_routes*(num_routes-1)):
        raise ValueError("The number of bred children each generation should be less than or equal to the number of unique parent combinations")

    num_locs = distances.shape[0]
    best_routes = [None]*generations

    routes, fitness = create_candidates(num_routes, num_locs, distances, rng)

    for gen in range(generations):
        have_not_bred = np.ones([len(routes)]*2,dtype=bool)
        have_not_bred[np.arange(len(routes)),np.arange(len(routes))] = False

        # pick partners in advance
        partners = -np.ones((2, num_routes**2), dtype=int)
        valid_partner = np.zeros(num_routes**2, dtype=bool)
        combined_fitness = np.zeros(num_routes**2)
        for i in range(num_routes):
            for j in range(num_routes):
                idx = i*num_routes + j

                valid_partner[idx] = not (i==j)
                if not valid_partner[idx]:
                    continue

                partners[0,idx] = i
                partners[1,idx] = j
                combined_fitness[idx] = fitness[i]+fitness[j]

        valid_idx = np.nonzero(valid_partner)[0]
        combined_fitness = (combined_fitness/combined_fitness.sum())
        partners_idx = rng.choice(valid_idx, size=num_bred, p=combined_fitness[valid_idx], replace=False)

        # breed the candidates to get children
        children = [0]*num_bred
        for b in range(num_bred):

            i1, i2 = partners[:, partners_idx[b]]
            
            # breed and mutate
            route1 = routes[i1]
            route2 = routes[i2]
            children[b] = breed(route1, route2, distances, rng)
            if mutation_chance < rng.random():
                mutate(children[b], rng)

        # choose new generation from among the children based on fitness
        children_distances = np.array([calculate_total_distance(child, distances) for child in children])
        children_fitness = calculate_fitness(children_distances)
        best_routes[gen] = children[np.argmax(children_fitness)]
        routes_idx = cull(children_fitness, num_routes, rng)
        routes = tuple(children[i] for i in routes_idx)
        fitness = children_fitness[routes_idx]

    return best_routes

def random_search(to_consider, distances, gens, rng):

    return np.array([calculate_total_distance(new_route(to_consider, rng), distances) for _  in range(gens)])

def random_breed(to_consider, distances, gens, rng):

    best_distances = np.ones(gens)
    for i in range(gens):

        route1 = new_route(to_consider, rng)
        route2 = new_route(to_consider, rng)
        child = breed(route1, route2, distances, rng)
        best_distances[i] = calculate_total_distance(child, distances)

    return best_distances

def random_breed_better(to_consider, distances, gens, rng):

    best_distances = np.ones(gens)
    best_routes = [None]*gens
    for i in range(gens):

        # if i == 0:
        #     route1 = new_route(to_consider, rng)
        # else:
        #     route1 = best_routes[np.argmin(distances)]

        route1 = new_route(to_consider, rng)

        route1_dist = calculate_total_distance(route1, distances)
        route2 = new_route(to_consider, rng)
        route2_dist = calculate_total_distance(route2, distances)
        child = breed(route1, route2, distances, rng)
        child_dist = calculate_total_distance(child, distances)
        routes = (route1, route2, child)
        routes_dist = (route1_dist, route2_dist, child_dist)
        best_idx = min((0,1,2), key=lambda x: routes_dist[x])
        best_routes[i] = routes[best_idx]
        best_distances[i] = routes_dist[best_idx]

    return best_distances

def test01():

    dat = np.loadtxt("wg59_dist.txt")
    print(dat.shape)

    to_consider = 11
    distances = dat[:to_consider, :to_consider]
    print(dat.shape)
    
    rng = np.random.default_rng()
    # print(new_order(to_consider, rng))

    for i in range(10):
        print("--------------------")
        route1 = new_route(to_consider, rng)
        print("First", route1, calculate_total_distance(route1, distances))
        route2 = new_route(to_consider, rng)
        print("Second", route2, calculate_total_distance(route2, distances))
        child = breed(route1, route2, distances, rng)
        print("Child", child, calculate_total_distance(child, distances))
        mutate(child, rng)
        print("Mutated", child, calculate_total_distance(child, distances))

def main():

    dat = np.loadtxt("wg59_dist.txt")
    print(dat.shape)

    to_consider = 11
    distances = dat[:to_consider, :to_consider]
    rng = np.random.default_rng()

    gens = 1000
    num_bred = 15
    best_in_gen = find_optimal_route(distances, 10, num_bred, 1/num_bred, gens, rng)
    best_distances = [calculate_total_distance(route, distances) for route in best_in_gen]
    random_distances = random_search(to_consider, distances, gens, rng)
    random_breed_distances = random_breed(to_consider, distances, gens, rng)
    random_breed_better_distances = random_breed_better(to_consider, distances, gens, rng)
    plt.plot(np.sort(best_distances), label="Breeding")
    plt.plot(np.sort(random_distances), label="Random")
    plt.plot(np.sort(random_breed_distances), label="Random Breed")
    plt.plot(np.sort(random_breed_better_distances), label="Random Breed Better")
    plt.legend()
    plt.show()

if __name__ == "__main__":

    main()