import numpy as np
import matplotlib.pyplot as plt
import os

def load_distances(filename):

    return np.loadtxt(filename)[:, 0]

def main():

    breed = load_distances("generated_data/breed.txt")
    random = load_distances("generated_data/random_search.txt")
    random_breed = load_distances("generated_data/random_breed.txt")
    random_breed_better = load_distances("generated_data/random_breed_better.txt")


    # get parallel files
    # ------------------------------

    # parallel algorithm with migration
    i = -1
    parallel_breeds1 = []
    while True:
        i += 1
        filename = f"generated_data/parallel_breed1_{i}.txt"
        if os.path.exists(filename):

            parallel_breeds1.append(load_distances(filename))
            continue
        
        break

    # parallel algorithm without migration (just the serial, but
    # in parallel processes)
    i = -1
    parallel_breeds2 = []
    while True:
        i += 1
        filename = f"generated_data/parallel_breed2_{i}.txt"
        if os.path.exists(filename):

            parallel_breeds2.append(load_distances(filename))
            continue
        
        break

    # =================================

    # plot all the data from the serial case
    idx = np.arange(len(breed))
    plt.plot(idx, np.sort(random), label=f"Random, min={np.min(random):.1f}", linewidth=2)
    idx += 1
    plt.plot(idx, np.sort(random_breed), label=f"Random Breed, min={np.min(random_breed):.1f}", linewidth=2)
    idx += 1
    plt.plot(idx, np.sort(random_breed_better), label=f"Random Breed Better, min={np.min(random_breed_better):.1f}", linewidth=2)
    idx += 1
    plt.plot(idx,np.sort(breed), label=f"Serial Breed, min={np.min(breed):.1f}", linewidth=2)

    # plot the parallel cases' data such that just one is plotted first,
    # then each loop adding the data from another process and choosing
    # from each generation the best distance
    parallel_breed = parallel_breeds1[0]
    for j in range(i):
        if 0 < j:
            parallel_breed = np.where(parallel_breeds1[j] < parallel_breed, parallel_breeds1[j], parallel_breed)
        idx += 1
        plt.plot(
            idx,np.sort(parallel_breed), 
            label=f"Migration Parallel Breed with {j+1} processes, min={np.min(parallel_breed):.1f}", 
            linewidth=2, linestyle="--"
        )

    parallel_breed = parallel_breeds2[0]
    for j in range(i):
        if 0 < j:
            parallel_breed = np.where(parallel_breeds2[j] < parallel_breed, parallel_breeds2[j], parallel_breed)
        idx += 1
        plt.plot(
            idx,np.sort(parallel_breed), 
            label=f"Parallel Breed with {j+1} processes, min={np.min(parallel_breed):.1f}",
            linewidth=2
        )

    plt.legend()
    plt.show()



if __name__ == "__main__":

    main()