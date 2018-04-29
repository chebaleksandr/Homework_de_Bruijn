from Bio import SeqIO
import Bio
from collections import defaultdict
from graphviz import Digraph

class Vertex: 
#
# Класс отвечающий за вершины графа, имеет атрибуты:
# seq - последовательность; 
# coverage - покрытие, показывает сколько раз встретился подобный узел
# in_edges - входящие ребра ( т.е. данная вершина является суффиксом некоторго ребра )
# out_edges - исходящие ребра ( т.е. данная вершина является префиксом некоторго ребра )
#
# Методы класса:
# increase_coverage - увеличивает покрытие на 1
#
   
    def __init__(self, seq):
        self.seq = seq
        self.coverage = 1
        self.in_edges = {}
        self.out_edges = {}
        
    def increase_coverage(self):
        self.coverage += 1


class Edge:
#
# Класс отвечающий за ребра графа, имеет атрибуты:
# seq - последовательность; 
# coverage - покрытие, показывает сколько раз встретился данное ребро
# in_edges - входящие ребра ( т.е. данная вершина является суффиксом некоторго ребра )
# out_edges - исходящие ребра ( т.е. данная вершина является префиксом некоторго ребра )
#
# Методы класса:
# increase_coverage - увеличивает покрытие на 1
#    
    def __init__(self,k1,k2):
        self.seq = k1 + k2[-1]
        self.start = k1
        self.end = k2
        self.n = 2
        self.coverage = 0
    
    def calc_coverage(self,c1,c2):
        self.coverage = (c1+c2)/2


class Graph:

    def __init__(self,k):
        self.vertices = {}
        self.edges_list = {}
        self.k = k
        
    def add_read(self,read):
        read_lng = len(read)
        if read_lng < self.k:
            return
            
        kmer = read[:k]
        if kmer in self.vertices:
            self.vertices[kmer].increase_coverage()
        else:
            self.vertices[kmer] = Vertex(kmer)
        
        for next_kmer_indx in range(1,read_lng-k+1,1):
            next_kmer = read[next_kmer_indx:(next_kmer_indx+k)]
            if next_kmer in self.vertices:
                self.vertices[next_kmer].increase_coverage()
            else:
                self.vertices[next_kmer] = Vertex(next_kmer)
            
            new_edge = Edge(kmer,next_kmer)
            
            new_edge_name = kmer + next_kmer[-1]
            #self.edges_list[name_edges] = [kmer, next_kmer]
            self.edges_list[new_edge_name] = new_edge
            
            self.vertices[next_kmer].in_edges[kmer]  = [new_edge]
            
            self.vertices[kmer].out_edges[next_kmer] = [new_edge]

            kmer = next_kmer
    
    def calc_init_edge_coverage(self):
        
        for current_vertex in self.vertices.keys():
            for next_vertex in self.vertices[current_vertex].out_edges.keys():
                self.vertices[current_vertex].out_edges[next_vertex][0].calc_coverage(self.vertices[current_vertex].coverage,self.vertices[next_vertex].coverage)
    
    def vizualize(self,file,mark):
        draw_g = Digraph()
        if mark == 'complex':
            for i in my_graph.vertices:
                draw_g.node(i,i)
            for j in my_graph.edges_list:
                e=my_graph.edges_list[j]
                draw_g.edge(e.start,e.end, label = str(j))
            draw_g.render(file)
        elif mark == 'simple':
            for i in my_graph.vertices:
                draw_g.node(i,str(my_graph.vertices[i].coverage))
            for j in my_graph.edges_list:
                e=my_graph.edges_list[j]
                draw_g.edge(e.start,e.end, label = 'Coverage: {}, length: {}'.format(my_graph.edges_list[j].coverage,len(my_graph.edges_list[j].seq)))
            draw_g.render(file)
            draw_g.render(file)
        else:
            print('Check mark!')
        

if __name__ == '__main__':
    
    dataset = '/home/aleksandr/Documents/my_fasta_exmpl.fasta'

    k = 6
    
    my_graph = Graph(k)
    
    with open(dataset, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            read = str(record.seq)
            my_graph.add_read(read)
#            my_graph.add_read( str(record.reverse_complement().seq) )

    my_graph.calc_init_edge_coverage()
    
    for v in my_graph.vertices:
        print('Vertex: {}, coverage: {}'.format(v,my_graph.vertices[v].coverage))
        for e in my_graph.vertices[v].out_edges:
            print('-> Out edge: {}'.format(e))
        for e in my_graph.vertices[v].in_edges:
            print('-> In edge: {}'.format(e))
            
    
    my_graph.vizualize('/home/aleksandr/Documents/round-table.gv','simple')
