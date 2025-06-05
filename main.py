from snappy import *
import csv
import sys
from datetime import datetime
from alexander_poly import alexander_presentation, alexander_nullity, alexander_polynomial
from helpers import eisermann

HEADER = ["PD Code", "Alexander Polynomial", "Satisfies Eisermann"]

def generate_links(num_links, file, max_crossings=30, components=2):
    writer = csv.writer(file)

    while i < 100:
        L = random_link(crossings = max_crossings, num_components=components, initial_map_gives_link=True, simplify="global")
        
        # Throw it out if it doesn't have enough components
        if (len(L.link_components) < components):
            continue

        # Check efficient obstructions first        
        if (L.signature() != 0) or (L.linking_matrix()[0][0] != 0):
            continue
    
        # Check Alexander nullity
        M = alexander_presentation(L)
        if alexander_nullity(L) != components - 1:
            continue

        # Keep track of the polynomial and Eisermann
        a_poly = alexander_polynomial(M, components - 1)
        satisfies_eisermann = eisermann(L)
        
        writer.writerow([L.PD_code(), a_poly, satisfies_eisermann])
        i += 1


def main():

    file_name = "Links_" + str(datetime.now())
    num_links = 0
    max_crossings = 30
    num_components = 2

    if len(sys.argv) < 3:
        print("Too few command line arguments")
        return
    if len(sys.argv) >= 3:
        num_links = int(sys.argv[2])
    if len(sys.argv) >= 4:
        file_name = sys.argv[3]
    if len(sys.argv) >= 5:
        max_crossings = int(sys.argv[4])
    if len(sys.argv) >= 6:
        num_components = int(sys.argv[5])
    

    with open(file_name, 'w', newline='') as out_file:
        generate_links(num_links, out_file, max_crossings=max_crossings, components=num_components)