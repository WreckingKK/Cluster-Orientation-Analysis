import math
import os
from shutil import copyfile
import xml.etree.cElementTree as ET
from io import StringIO
import sys

sys.setrecursionlimit(100000)  

class XMLs:
    def __init__(self, file_name=None):
        self.all_beads_point = None
        self.tree = ET.parse(file_name)
        self.filename = file_name

    @staticmethod
    def __gain_config(tree):
        gala_node = tree.getroot()
        configuration = list(gala_node)[0]
        return configuration

    def __get_one_title_node(self, title):
        return self.__gain_config(self.tree).find(title)


    def __record_all_bead_point_and_index(self):
        self.all_beads_point = []
        t = self.__get_one_title_node('type')
        p = self.__get_one_title_node('position')
        t_IO = StringIO(t.text)
        p_IO = StringIO(p.text)
        t_IO.readline()
        p_IO.readline()

        n = int(p.attrib['num'])
        k = 0
        while k < n:
            t_l = t_IO.readline()
            p_l = p_IO.readline()
            t_l = t_l.replace('\n', '')
            trans_p_l = []
            p_l = p_l.replace('\n', '').split(' ')
            for p_i in p_l:
                if p_i:
                    trans_p_l.append(p_i)
                trans_p_l.append(k)
                trans_p_l.append(t_l)
                self.all_beads_point.append(trans_p_l)
            k += 1

    def first(self, data):
        t = self.__get_one_title_node('type')
        t_IO = StringIO(t.text)
        n = int(t.attrib['num'])

        k = 0
        new_t = StringIO()
        new_t.write(t_IO.readline())
        while k < n:
            new_t.write(data[k] + '\n')
            k += 1
        t.text = new_t.getvalue()
        self.tree.write(self.filename)

def cp_filename(f_i):
    f_i = f_i.split('\\')
    f_i[-1] = 'cp_' + f_i[-1]
    cp_f_i = '\\'.join(f_i)
    return cp_f_i
def new_filename(f_i, head_n):
    f_i = f_i.split('\\')
    f_i[-1] = head_n + f_i[-1]
    new_f_i = '\\'.join(f_i)
    return new_f_i

def data_extract_Clip(f_i, data):
    print("Clip to xml : ", f_i)
    cp_f_i = cp_filename(f_i)

    if copyfile(f_i, cp_f_i):
        xml = XMLs(cp_f_i)
        # xml.second(data)
        xml.first(data)
        old_file = open(cp_f_i)
        f = open(new_filename(f_i, 'Clip_'), 'w')
        f.write('<?xml version="1.0" encoding="UTF-8"?>\n' + old_file.read())
        f.close()
        old_file.close()
        os.remove(cp_f_i)

def data_extract_Simple(f_i, data):
    print("four to xml : ", f_i)
    cp_f_i = cp_filename(f_i)

    if copyfile(f_i, cp_f_i):
        xml = XMLs(cp_f_i)
        xml.first(data)
        old_file = open(cp_f_i)
        f = open(new_filename(f_i, 'S'), 'w')
        f.write('<?xml version="1.0" encoding="UTF-8"?>\n' + old_file.read())
        f.close()
        old_file.close()
        os.remove(cp_f_i)

def data_extract_Direciton(f_i, data):
    print("four to xml : ", f_i)
    cp_f_i = cp_filename(f_i)

    if copyfile(f_i, cp_f_i):
        xml = XMLs(cp_f_i)
        xml.first(data)
        old_file = open(cp_f_i)
        f = open(new_filename(f_i, 'X'), 'w')
        f.write('<?xml version="1.0" encoding="UTF-8"?>\n' + old_file.read())
        f.close()
        old_file.close()
        os.remove(cp_f_i)
        
def data_extract_BtoC(f_i, data):
    print("Clip to xml : ", f_i)
    cp_f_i = cp_filename(f_i)

    if copyfile(f_i, cp_f_i):
        xml = XMLs(cp_f_i)
        # xml.second(data)
        xml.first(data)
        old_file = open(cp_f_i)
        f = open(new_filename(f_i, 'BtoC_'), 'w')
        f.write('<?xml version="1.0" encoding="UTF-8"?>\n' + old_file.read())
        f.close()
        old_file.close()
        os.remove(cp_f_i)