"""
Cluster Division using NetworkX:

Use networkx to divide the graph into multiple clusters.

Iterate through each cluster and proceed to step 2.

Recursive Search for Parallel Monomers:

Starting from a single monomer, recursively search for parallel monomers and mark them as the same class.

Proceed to step 3.

Defect Condition Definition:

Define the conditions for a defect.

Determine if the smaller clusters qualify as defects.

Proceed to step 4.

Angle Calculation for Non-Defect Clusters:

Calculate the angles for non-defect clusters and store them in a list.

Output Image

"""

import time
import os
import sys
import math
import numpy as np
from io import StringIO
from xml.etree import cElementTree
from pandas import read_csv
import networkx as nx
from joblib import Parallel, delayed
from defect_to_xml_clip import data_extract_Clip


class Xml(object):
    def __init__(self, filename):
        self.root = cElementTree.ElementTree(file=filename).getroot()[0]
        self.G = nx.Graph()
        self.nodes = {}

        tag = ["bond", "type", "position", "box"]
        for ta in tag:
            e = self.root.find(ta)
            if ta == "box":
                self.nodes[ta] = np.array([float(e.attrib['lx']), float(e.attrib['ly']), float(e.attrib['lz'])])
            else:
                self.nodes[ta] = read_csv(StringIO(e.text), sep='\\s+', header=None).values

        bonds = self.nodes.get('bond', None)
        if bonds is not None:
            for bond in bonds:
                ptc1, ptc2 = bond[1], bond[2]
                self.G.add_edge(ptc1, ptc2)

    def run(self):
        return self.nodes, self.G

def pbc(p, box):
    return p - box * np.round(p / box)

def VEC(a, b, box):
    # unwrapping vector ab
    return pbc(a - b, box)

def cos_angle(v1, v2):
    return math.degrees(
        math.acos((v1[0] * v2[0] + v1[1] * v2[1]) / (np.sqrt(np.sum(v1 ** 2)) * np.sqrt(np.sum(v2 ** 2)))))

def info_link_bond(type_i, G):
    # initialize
    link_points = {k: [] for k in range(len(type_i))}

    # find the link particle
    for i in G.edges:
        link_points[i[0]].append(i[1])
        link_points[i[1]].append(i[0])
    return link_points

def push_search(li, link_points_original, type_i, position_i, box):
    def find_1_monomer(li, type_i, link_points): #li: each bonded cluster
        # find A_A
        A_A = []
        find = False
        for i in li:
            if type_i[i] == 'A':
                for j in link_points[i]:
                    if type_i[j] == 'A':
                        A_A.append(i)
                        A_A.append(j)
                        find = True
                        break
            if find:
                break
        temB = []
        temC = []
        for i in A_A:
            for j in link_points[i]:
                if type_i[j] == 'B':
                    temB.append(j)
                elif type_i[j] == 'C':
                    temC.append(j)

        return [A_A, temB, temC]

    def find_1(type_i, link_points, idx):
        #  find A_A
        A_A = []
        find = False
        for i in link_points[idx]:
            if type_i[i] == 'A':
                for j in link_points[i]:
                    if type_i[j] == 'A':
                        A_A.append(i)
                        A_A.append(j)
                        find = True
                        break
            if find:
                break

        temB = []
        temC = []
        for i in A_A:
            for j in link_points[i]:
                if type_i[j] == 'B':
                    temB.append(j)
                elif type_i[j] == 'C':
                    temC.append(j)

        return [A_A, temB, temC]

    def push_(monomers, link_points_original, tag_c, tag_aa):
        precision = 10
        parallel_branch = []
        parallel_monomer = []
        for i0 in monomers:
            for i1 in i0[2]:  #i0[2]:"C" particle list in each monomer
                if tag_c[i1] == 0:
                    tag_c[i1] = 1
                    for j1 in link_points_original[i1]:
                        # find the monomer linked with i0
                        if type_i[j1] == 'C':
                            tag_c[j1] = 1
                            monomer1 = find_1(type_i, link_points_original, j1)

                            parallel = False
                            a1 = i0[0][0]
                            a11 = i0[0][1]
                            a1_ = monomer1[0][0]
                            a11_ = monomer1[0][1]   #AA of two bonded monomers

                            tem = sorted([a1_, a11_])
                            tem1 = sorted([a1, a11])
                            if tag_aa[str(tem[0]) + "-" + str(tem[1])] == 0:
                                AA1_ = VEC(position_i[a1], position_i[a11], box)
                                AA2_ = VEC(position_i[a1_], position_i[a11_], box)
                                angle = cos_angle(AA1_, AA2_)
                                # print("angle : ", angle, 180 - angle)
                                if angle < precision or 180 - angle < precision:
                                    parallel = True
                                if parallel:
                                    tag_aa[str(tem[0]) + "-" + str(tem[1])] = tag_aa[str(tem1[0]) + "-" + str(tem1[1])]
                                    #parallel: {'A6-A9':A0,A3},{'A24-A27':A0,A3}.....{'all pl aa':A0,A3}
                                else:
                                    tag_aa[str(tem[0]) + "-" + str(tem[1])] = monomer1[0]
                                    #not parallel: {'A12-A15':A12,A15}
                                parallel_monomer.append(monomer1)
                            else:
                                if tag_aa[str(tem[0]) + "-" + str(tem[1])] == tag_aa[str(tem1[0]) + "-" + str(tem1[1])]:
                                    pass
                                else:
                                    AA1_ = VEC(position_i[a1], position_i[a11], box)
                                    AA2_ = VEC(position_i[a1_], position_i[a11_], box)
                                    angle = cos_angle(AA1_, AA2_)
                                    if angle < precision or 180 - angle < precision:
                                        tem2 = tag_aa[str(tem[0]) + "-" + str(tem[1])] + tag_aa[
                                            str(tem1[0]) + "-" + str(tem1[1])]
                                        parallel_branch.append(tem2)
        return parallel_monomer, parallel_branch

    def push_defect(monomers, link_points_original, type_i, tag_add, tag_c):
        # 开始 递归
        connect_monomer = []
        push_num = 0
        for i0 in monomers:
            for i1 in i0[2]: # traverse c particles in defects
                if tag_c[i1] == 0:
                    tag_c[i1] = 1
                    for j1 in link_points_original[i1]:
                        # get the next monomer linked with defect c particle
                        if type_i[j1] == 'C':
                            if tag_add[j1] == -1:
                                if tag_c[j1] == 0:
                                    tag_c[j1] = 1
                                    monomer1 = find_1(type_i, link_points_original, j1)
                                    monomer1 += [tag_add[j1]]
                                    connect_monomer.append(monomer1)
                                    push_num += 1
                            else:
                                i0 += [tag_add[j1]]
                                connect_monomer.append(i0)
        return connect_monomer, push_num


    # 1. find the monomer 1
    monomer = find_1_monomer(li, type_i, link_points_original)
    # 2. find the parallel AA cluster as key,value in tag_aa
    tag_c = {}
    for i in li:
        if type_i[i] == 'C':
            tag_c[i] = 0
    tag_aa = {}
    for i in li:
        if type_i[i] == 'A':
            for j in link_points_original[i]:
                if type_i[j] == 'A':
                    tem = sorted([i, j])
                    str_aa = str(tem[0]) + "-" + str(tem[1])
                    tag_aa[str_aa] = 0
                    break

    # start
    tem = sorted([monomer[0][0], monomer[0][1]]) #A-A
    tag_aa[str(tem[0]) + "-" + str(tem[1])] = monomer[0]
    parallel_monomer = [monomer]
    parallel_monomers = [monomer]
    parallel_branches = []
    while True:
        parallel_monomer, parallel_branch = push_(parallel_monomer, link_points_original, tag_c, tag_aa)
        parallel_branches += parallel_branch
        parallel_monomers += parallel_monomer
        if len(parallel_monomer) == 0:
            break

    # print("parallel_branches : ", parallel_branches)

    # link the nodes of branch AA bond
    inner_G = nx.Graph()
    for i in parallel_branches:
        inner_G.add_edge(i[0], i[1])
        inner_G.add_edge(i[1], i[2])
        inner_G.add_edge(i[2], i[3])

    # reverse the dic of tag_aa{key:value}, re_tag_aa{value:key1,key2,key3....}
    re_tag_aa = {}
    for k, v in tag_aa.items():
        tem = sorted(v)
        re_tag_aa[str(tem[0]) + "-" + str(tem[1])] = []

    for k, v in tag_aa.items():
        tem = sorted(v)
        re_tag_aa[str(tem[0]) + "-" + str(tem[1])].append(k)

    # construct new clusters included branch A-A bond
    clusters = []
    for i in nx.connected_components(inner_G):
        tem_l = []
        tem_d = {str(j): 0 for j in i}
        for j1 in tem_d:
            if tem_d[j1] == 0:
                for j2 in re_tag_aa:
                    str1 = j2.split('-')  # str ['11820', '11823']
                    if j1 == str1[0] or j1 == str1[1]:
                        tem_d[str1[0]] += 1
                        tem_d[str1[1]] += 1
                        tem_l += re_tag_aa[j2]
                        re_tag_aa[j2] = []

        clusters.append([len(tem_l), tem_l])

    # re_tag_aa: clusters without branch A-A bond
    # thus, merge re_tag_aa into clusters
    for k, v in re_tag_aa.items():
        if len(v) > 0:
            clusters.append([len(v), v])   # all clusters with different orientations [[len,['A-A','A-A'....]]]


    parallel_monomers_d = {}
    for i in parallel_monomers:
        tem = sorted(i[0])
        parallel_monomers_d[str(tem[0]) + "-" + str(tem[1])] = i

    re_clusters = []
    for i in clusters:
        tem = []
        for j in i[1]:
            tem.append(parallel_monomers_d[j])
        re_clusters.append([i[0], tem]) # all clusters with different orientations [[len,[[[A,A],[B],[C,C,C]],[[A,A],[B,B],[C,C]]]...]]


    re_clusters = sorted(re_clusters)
    max_cluster = re_clusters[-1]
    add_cluster = []
    defect_cluster = []

    if max_cluster[0] < 10:
        defect_cluster.append(max_cluster)
        for i in re_clusters[:-1]:
            defect_cluster.append(i)
    else:
        add_cluster.append(max_cluster)
        for i in re_clusters[:-1]:
            if i[0] / max_cluster[0] > 0.5:
                add_cluster.append(i)
            else:
                if i[0] >= 10:
                    add_cluster.append(i)
                else:
                    defect_cluster.append(i)

    de_sum = 0
    de_cluster = []
    for i in defect_cluster:
        de_sum += i[0]
        de_cluster += i[1]

    tag_add = {}
    for i in range(len(add_cluster)):
        id_ = add_cluster[i][1] # each cluster
        for j in id_:  # each monomer in a cluster [[14939, 14940],[],[14944, 14945, 14941, 14942]]
            for j1 in j: # each type in a monomer  [14944, 14945, 14941, 14942]
                for j2 in j1: # each particle in a monomer  14944
                    tag_add[j2] = i
    # tag_add:{each particle:dependent cluster id}

    tag_c = {}
    for i in de_cluster:
        ic = i[2]
        for j in ic:
            tag_c[j] = 0
        for j in i:
            for j1 in j:
                tag_add[j1] = -1
    # tag_c:{each c particle in defect part:0}
    # tag_add:{each particle in defect part:-1}
    defect_monomersss = []
    for i in defect_cluster:
        defect_monomers = i[1]
        defect_monomerss = []
        push_nums = 0
        while True:
            defect_monomers, push_num = push_defect(defect_monomers, link_points_original, type_i, tag_add, tag_c)
            defect_monomerss += defect_monomers
            push_nums += push_num
            if push_num == 0:
                break
        if defect_monomerss:
            defect_monomersss.append(defect_monomerss)
            # print(len(defect_monomerss))
            # print(defect_monomerss)
            # print('---------------------------------')

    defect_offshoot = {}
    for i in range(len(add_cluster)):
        defect_offshoot[i] = []

    # print(parallel_monomers_d)

    for i in defect_monomersss:
        tail_ = {}
        for j in i:
            for j1 in j[3:]:
                if j1 != -1:
                    tail_[j1] = 1

        for j in tail_.keys():
            cluster_d = {}
            for j1 in i:
                j10 = j1[0]
                tem = sorted(j10)
                cluster_d[str(tem[0]) + "-" + str(tem[1])] = 1

            for key in cluster_d:
                defect_offshoot[j].append(parallel_monomers_d[key])
    return add_cluster, [de_sum, de_cluster], defect_offshoot
def angle_four(add_cluster, position_i, box):
    # The angle of AA relative to the x-axis
    re_add_clusters = []
    x = np.array([1, 0])

    for i in add_cluster:
        angles = []
        for j in i[1]:
            a1 = j[0][0]
            a2 = j[0][1]

            aa1 = VEC(position_i[a1], position_i[a2], box)
            if aa1[0] > 0: #<180
                aa1 = VEC(position_i[a2], position_i[a1], box)
            angles.append(cos_angle(aa1, x))
        print("angle : ", max(angles), min(angles))
        re_add_clusters.append([i[0], max(angles), min(angles), i[1]])
    return re_add_clusters


def Show(type_i, clusters, defects, file_name):
    data = {key: type_i[key] for key in type_i}

    tag = 0
    for k, v in clusters.items():
        tag1 = 0              
        for i1 in v:
            add_cluster = i1[-1]
            for i in add_cluster:
                for j1 in i[0]:
                    data[j1] = type_i[j1] + str(tag) + "_" + str(tag1)
                for j1 in i[1]:
                    data[j1] = type_i[j1] + str(tag) + "_" + str(tag1)
                for j1 in i[2]:
                    data[j1] = type_i[j1] + str(tag) + "_" + str(tag1)
            tag1 += 1
        tag += 1
        
    for k, v in defects.items():
        for j in v[1]:
            for j1 in j[0]:
                data[j1] = type_i[j1] + 'F'
            for j1 in j[1]:
                data[j1] = type_i[j1] + 'F'
            for j1 in j[2]:
                data[j1] = type_i[j1] + 'F'

    data_extract_Clip(file_name, data)


def data_extract(file_name, k):
    print("four direction : ", file_name, k)

    xml = Xml(file_name)
    nodes, G = xml.run()
    all_mono = len(xml.nodes["type"]) / 6
    sim_time = float(xml.root.attrib['time_step']) * 0.0005


    type_i = {k: v[0] for k, v in zip(range(len(nodes["type"])), nodes["type"])}
    position_i = {k: v[:2] for k, v in zip(range(len(nodes["position"])), nodes["position"])}
    box = nodes["box"]
    box = box[:2]

    link_points_original = info_link_bond(type_i, G)


    clusters = {}
    defects = {}
    offshoots = {}
    size_all_cls = 0
    size_all_def = 0
    free_mono = 0
    num_cls = 0
    k = 0
    for i in nx.connected_components(G):
        if len(i) >= 120:
            add_cluster, defect_cluster, defect_offshoot = push_search(i, link_points_original, type_i, position_i, box)

            for j in add_cluster:
                size_all_cls += j[0]
            size_all_def += defect_cluster[0]
            num_cls += len(add_cluster)

            add_cluster = angle_four(add_cluster, position_i, box)
            clusters["cluster" + str(k)] = add_cluster
            defects["defect" + str(k)] = defect_cluster  # [de_sum, de_cluster]
            offshoots["offshoot" + str(k)] = defect_offshoot
            k += 1

        else:
            free_mono += len(i) / 6
            
    Show(type_i, clusters, defects, file_name)
    return xml.root.attrib['time_step'], clusters, defects, offshoots, size_all_cls/all_mono, size_all_def/all_mono, int(free_mono)/all_mono, sim_time, num_cls

# def jaccard_simi(id):

def separate_file(filepath):
    files = []
    for file_name in os.listdir(filepath):
        if file_name.startswith('BtoC_p.') and file_name.endswith('.xml'):
            input_path = os.path.join(filepath, file_name)
            files.append(input_path)
    return files


def run():
    file_l = sys.argv[1:]

    
    Back_Value = Parallel(n_jobs=1)([delayed(data_extract)(v, k) for k, v in enumerate(file_l)])
    # print('Back_Value: ', Back_Value)

    File1 = None
    File2 = None
    File3 = None
    File4 = None

    try:
        File1 = open("sql_alz.log", "w")
        # all clusters
        File2 = open("Clip_data.log", "w")
        File3 = open("cls_def.log", "w")
        File4 = open("def_corr.log", "w")

        File1.write('time\tAll_cls\tAll_def\tFree_mono\tNum_cls\n')
        File3.write('cls_size\tdef_size\n')
        File4.write('time\tcorr_id\tdef_num\n')

        for Back in Back_Value:
            time_step = Back[0]
            clusters = Back[1]
            defects = Back[2]
            offshoots = Back[3]
            size_all_cls = Back[4]
            size_all_def = Back[5]
            free_mono = Back[6]
            sim_time = Back[7]
            num_cls = Back[8]

            File1.write(f"{sim_time}\t{size_all_cls}\t{size_all_def}\t{free_mono}\t{num_cls}\n")

            File2.write("timestep: " + str(time_step) + "\n")

            #File3.write("time:" + str(sim_time) + "\n")
            #File3.write('cls_size\tdef_size\n')
            File4.write(f"{sim_time}\t")

            k1 = 0
            all_def_num = 0
            all_def_id = []
            while k1 < len(clusters):
                clustervo = clusters["cluster" + str(k1)]
                defectv = defects["defect" + str(k1)]
                offshootv = offshoots["offshoot" + str(k1)]

                k2 = 0
                File2.write("cluster" + str(k1) + "\n")

                for cvo in clustervo:
                    # File 2
                    File2.write("No " + str(k1) + " " + str(k2) + "\n")
                    File2.write(str(cvo[0]) + " " + str(cvo[1]) + " " + str(cvo[2])
                                + " " + str((cvo[1] + cvo[2]) / 2) + "\n")
                    File3.write(f"{cvo[0]}\t")

                    clustervi = cvo[3]
                    for j0 in range(len(clustervi)):

                        for j1 in range(len(clustervi[j0])):
                            for j2 in clustervi[j0][j1]:
                                File2.write(str(j2) + " ")

                            if j1 != len(clustervi[j0]) - 1:
                                File2.write(",")

                        if j0 != len(clustervi) - 1:
                            File2.write(" | ")

                    File2.write("\n")

                    # offshoot
                    offo = offshootv[k2]
                    File2.write("Offshoot " + str(k1) + " " + str(k2) + "\n")
                    File2.write(str(len(offo)) + "\n")
                    File3.write(f"{len(offo)}\n")
                    if len(offo) > 0:
                        for j0 in range(len(offo)):
                            offo[j0] = offo[j0][:3]

                            for j1 in range(len(offo[j0])):
                                for j2 in offo[j0][j1]:
                                    File2.write(str(j2) + " ")

                                if j1 != len(offo[j0]) - 1:
                                    File2.write(",")

                            if j0 != len(offo) - 1:
                                File2.write(" | ")
                        File2.write("\n")

                    k2 += 1

                dvnum = defectv[0]
                File2.write("defect" + str(k1) + "\n")
                File2.write(str(dvnum) + "\n")

                all_def_num += dvnum

                defectvi = defectv[1]
                for j0 in range(len(defectvi)):
                    if dvnum > 0 :

                        defectvi[j0] = defectvi[j0][:3]

                        id_v = defectvi[j0][0]
                        all_def_id.append(id_v)

                        for j1 in range(len(defectvi[j0])):
                            for j2 in defectvi[j0][j1]:
                                File2.write(str(j2) + " ")

                            if j1 != len(defectvi[j0]) - 1:
                                File2.write(",")

                        if j0 != len(defectvi) - 1:
                            File2.write(" | ")

                    File2.write("\n")

                k1 += 1

            File4.write(f"{all_def_id}\t{all_def_num}\n")
    finally:
        File1.close()
        File2.close()
        File3.close()
        File4.close()

if __name__ == '__main__':
    start_time = time.time()
    run()
    end_time = time.time()
    duration = end_time - start_time
    print('duration_time', duration)
