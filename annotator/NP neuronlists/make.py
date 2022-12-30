from os import getcwd
from os.path import join, split
from glob import glob
import json

def main(path):

    neuron_lists = {}

    files = glob(join(path, "*Neurons.txt"))

    for file in files:

        (_, listname) = split(file)
        listname = listname[:-4]
        
        with open(file, 'r') as f:
            neurons = f.readlines()

        neurons = [x.strip() for x in neurons]
        neuron_lists[listname] = neurons

    with open(join(path, 'neuronlists.json'), 'w') as outfile:
        json.dump(neuron_lists, outfile)

if __name__ == "__main__":
    path = getcwd()
    main(path)