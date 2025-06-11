from snappy import *
import csv
import sys
from datetime import datetime
import test

HEADER = ["PD Code", "Num Components", "Num Crossings", "Signature 0", "Linking Matrix 0", "Alexander Polynomial", "Satisfies Eisermann"]

def generate_links(num_links, file, max_crossings=30, components=2):
    writer = csv.writer(file)

    for i in range(num_links):
        L = random_link(crossings = max_crossings, num_components=components, initial_map_gives_link=True, simplify="global")
        if len(L.link_components) < 2:
            continue
        test(L , writer)


def parse_commands():
    # Defaults
    file_name = "Links_" + str(datetime.now()) + ".csv"
    num_links = 0
    max_crossings = 30
    num_components = 2
    
    num_links_input = input("How many links would you like to generate?\n")
    if num_links_input != '':
        num_links = int(num_links_input)

    file_name_input = input("Where do you want the data to be stored?\n")
    if file_name_input != '':
        file_name = file_name_input
    
    max_crossings_input = input("What is the maximum number of crossing you want?\n")
    if max_crossings_input != '':
        max_crossings = int(max_crossings_input)

    num_components_input = input("How many components do you want in your links?\n")
    if num_components_input != '':
        num_components = int(num_components_input)

    return file_name, num_links, max_crossings, num_components

def main():
    file_name, num_links, max_crossings, num_components = parse_commands()

    with open(file_name, 'w', newline='') as out_file:
        generate_links(num_links, out_file, max_crossings=max_crossings, components=num_components)

if __name__ == '__main__':
    main()