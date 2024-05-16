import numpy as np
import matplotlib.pyplot as plt
import os


def main():

    breed = np.loadtxt("generated_data/breed.txt")
    random = np.loadtxt("generated_data/random_search.txt")
    random_breed = np.loadtxt("generated_data/random_breed.txt")
    random_breed_better = np.loadtxt("generated_data/random_breed_better.txt")


    # get parallel files
    i = -1
    parallel_breeds = []
    while True:
        i += 1
        filename = f"generated_data/parallel_breed{i}.txt"
        if os.path.exists(filename):

            parallel_breeds.append(np.loadtxt(filename))
            continue
        
        break

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
    parallel_breed = parallel_breeds[0]
    for j in range(i):
        if 0 < j:
            parallel_breed = np.where(parallel_breeds[j] < parallel_breed, parallel_breeds[j], parallel_breed)
        idx += 1
        plt.plot(idx,np.sort(parallel_breed), label=f"Parallel Breed with {j+1} processes, min={np.min(parallel_breed):.1f}", linewidth=2)

    plt.legend()
    plt.show()



if __name__ == "__main__":

    main()