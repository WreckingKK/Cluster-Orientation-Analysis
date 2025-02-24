import time
from io import StringIO
from xml.etree import cElementTree
from pandas import read_csv
import networkx as nx
import os
import sys
from joblib import Parallel, delayed
from defect_to_xml_clip import data_extract_BtoC


class Xml(object):
    def __init__(self, filename):
        self.root = cElementTree.ElementTree(file=filename).getroot()[0]
        self.G = nx.Graph()
        self.nodes = {}

        tag = ["bond", "type"]
        for ta in tag:
            e = self.root.find(ta)
            self.nodes[ta] = read_csv(StringIO(e.text), sep='\\s+', header=None).values

        bonds = self.nodes.get('bond', None)
        if bonds is not None:
            for bond in bonds:
                ptc1, ptc2 = bond[1], bond[2]
                self.G.add_edge(ptc1, ptc2)

    def run(self):
        return self.nodes, self.G


def data_extract(file_name, k):
    print("four direction : ", file_name, k)

    xml = Xml(file_name)
    nodes, G = xml.run()

    """ 存信息 """
    type_i = {k: v[0] for k, v in zip(range(len(nodes["type"])), nodes["type"])}

    """ B to C """
    data = {k: v for k, v in type_i.items()}

    for i in G.degree:
        if i[1] != 1:
            if data[i[0]] == 'B':
                data[i[0]] = 'C'

    data_extract_BtoC(file_name, data)


def separate_file(filepath):
    files = []
    for file_name in os.listdir(filepath):
        if file_name.startswith('p.') and file_name.endswith('.xml'):
            input_path = os.path.join(filepath, file_name)
            files.append(input_path)
    return files


def run():
    file_l = sys.argv[1:]
    Back_Value = Parallel(n_jobs=1)([delayed(data_extract)(v, k) for k, v in enumerate(file_l)])
    print('Back_Value: ', Back_Value)


if __name__ == '__main__':
    start_time = time.time()
    run()
    end_time = time.time()
    duration = end_time - start_time
    print('duration_time', duration)
