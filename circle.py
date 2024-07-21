import numpy as np

def main():

    circ = np.array(
        [(np.cos(2*np.pi*i), np.sin(2*np.pi*i)) for i in np.linspace(0,1,51)[:-1]]
    )

    print(circ)

    np.savetxt("circle_points.txt", circ)

    print(np.linalg.norm(circ[0,:]))

    circ_dist = -1*np.ones([len(circ)]*2)
    for i in range(len(circ)):

        circ_dist[i,:] = np.linalg.norm(circ-circ[i,:], axis=1)

    with open("circle_dist.txt", "w") as f:
        print(len(circ_dist), file=f)
        np.savetxt(f, circ_dist, fmt="%10.9f")




if __name__ == "__main__":

    main()