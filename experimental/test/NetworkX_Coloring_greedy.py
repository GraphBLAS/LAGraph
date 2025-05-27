import networkx as nx
import time
from scipy.io import mmread

files = [
    "/raid/matrices/karate/karate.mtx",
    "/raid/matrices/Harvard500/Harvard500.mtx",
    "/raid/matrices/USpowerGrid/USpowerGrid.mtx",
    "/raid/matrices/as-Skitter/as-Skitter.mtx",
    "/raid/matrices/com-Youtube/com-Youtube.mtx",
    "/raid/matrices/com-LiveJournal/com-LiveJournal.mtx",
    "/raid/matrices/com-Orkut/com-Orkut.grb",
    "/raid/matrices/com-Friendster/com-Friendster.grb"
]

for f in files:
    print("==========================reading file", f)

    A = mmread(f)
    G = nx.Graph(A)

    print("done reading. exececuting...")

    t = time.perf_counter()
    d = nx.coloring.greedy_color(G)
    t = time.perf_counter() - t

    num_colors = len(set(d.values()))

    print("==========================time      :", t, "sec")
    print("==========================num colors:", num_colors)