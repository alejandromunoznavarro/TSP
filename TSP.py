#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Alejandro Muñoz Navarro
"""

import os
import sys
import math
from matplotlib import pyplot as plt
import numpy as np

import random 

def leer_fichero(file_path):
    
    # Lista para guardar las posiciones
    positions = []
    if not os.path.isfile(file_path):
        print("ERROR - no se ha encontrado el fichero")
        sys.exit()

    with open(file_path, "r") as file:
        for line in file:
            if line is not None:
                if line[0].isdigit():

                    line = line.strip()
                    line_elements = line.split()
                    positions.append((float(line_elements[1]), float(line_elements[2])))

    # Extrae el fichero
    file_name = os.path.basename(file_path)

    return positions, file_name

def get_distances(positions:list):
    distances = np.zeros((len(positions),len(positions)))
    npuntos = len(positions)
    for i in range(npuntos):
        for j in range(npuntos):
            if positions[j][0]==positions[i][0] and positions[j][1] == positions[i][1]:
                aux = np.inf
            else:
                aux = math.sqrt(math.pow((positions[j][0]-positions[i][0]),2)+math.pow((positions[j][1]-positions[i][1]),2))
            distances[i][j] = aux
    return distances

def pintar(positions,path):
    for position in positions:
        plt.plot(position[0],position[1], marker="o", color="red")
    x = []
    y = []
    for i in path:
        position = positions[i]
        x.append(position[0])
        y.append(position[1])
    
    plt.plot(x,y,linestyle='solid',color='green')
    plt.show()
    
class Grafo(object):
    def __init__(self, distances: list, num_hormigas, alpha, beta,Q):
        self.distances = distances
        self.npuntos = len(distances)
        self.visitados = []
        self.permitidos = [i for i in range(0,self.npuntos)]
        self.alpha = alpha
        self.beta = beta
        self.Q = Q
        
        # Inicializamos el grafo con un nodo aleatorio
        inicial = random.randint(0,self.npuntos-1)
        # Añadimos el nodo inicial a la lista de visitados
        self.visitados.append(inicial)
        # Ponemos el puntero del nodo actual en el nodo inicial
        self.actual = inicial
        # Eliminamos el nodo inicial de la lista de permitidos
        self.permitidos.remove(inicial)
        # Inicializamos la distancia total a 0
        self.distancia_total = 0

        # Calculamos el camino escogiendo aquel camino que sea el más corto
        while len(self.permitidos) != 0:
            minimo = np.inf
            siguiente = -1
            for i in self.permitidos:
                if self.distances[self.actual][i]<minimo:
                    minimo = self.distances[self.actual][i]
                    siguiente = i
            self.permitidos.remove(siguiente)
            self.visitados.append(siguiente)
            self.distancia_total += minimo
            self.actual = siguiente

        # Añadimos de nuevo el nodo inicial para cerrar el camino
        self.visitados.append(inicial)
        self.distancia_total += self.distances[self.actual][inicial]

        # Inicializamos las feromonas como (num_hormigas*Q)/(npuntos * distancia_camino_minimo)
        self.feromonas = [[(num_hormigas*self.Q)/(self.npuntos*self.distancia_total) for i in range(self.npuntos)] for j in range(self.npuntos)] 

class Colonia(object):
    def __init__(self,num_hormigas,rho):
        self.num_hormigas = num_hormigas
        self.rho = rho
        
    # Actualiza las feromonas
    # τxy = (1-ρ)τxy + ΣΔτxy
    def actualiza_feromonas(self, grafo:Grafo, hormigas: list):
        for i, fila in enumerate(grafo.feromonas):
            for j, columna in enumerate(fila):
                grafo.feromonas[i][j] *= 1-self.rho
                for hormiga in hormigas:
                    grafo.feromonas[i][j] += hormiga.feromonas[i][j]
    
    def resolver(self,grafo:Grafo):
        minimo = np.inf
        solucion = []
        sinMejora = 0
        
        # COndición de parada
        while sinMejora<20:
            hormigas = [Hormiga(self,grafo) for id_hormiga in range(self.num_hormigas)]
            for hormiga in hormigas:
                for nodo in range(grafo.npuntos - 1):
                    hormiga.selecciona_nodo()
                    
                hormiga.distancia_total += grafo.distances[hormiga.visitados[-1]][hormiga.visitados[0]]
                hormiga.visitados.append(hormiga.visitados[0])
                if hormiga.distancia_total < minimo:
                    minimo = hormiga.distancia_total
                    solucion = hormiga.visitados
                    sinMejora = 0
                else:
                    sinMejora += 1
                # Calcula las delta_feromonas (Δτxy) para el cálculo que realizaremos a continuación
                hormiga.delta_feromonas()
                
            # Una vez que la hormiga llega a su destino, se debe actualizar las feromonas de cada camino 
            # respecto al desgaste o evaporación ρ y si pasó la hormiga por el camino dejando feromonas.
            self.actualiza_feromonas(grafo,hormigas)
            
        return solucion, minimo
    
class Hormiga(object):
    def __init__(self, colonia: Colonia, grafo: Grafo):
        #self.colonia = colonia # colonia
        self.grafo = grafo
        self.distancia_total = 0.0
        self.visitados = []
        self.feromonas = []
        self.permitidos = [i for i in range(grafo.npuntos)]
        
        # Situamos a la hormiga en un nodo aleatorio del grafo
        inicio = random.randint(0, grafo.npuntos - 1)
        # Añadimos el nodo a la lista de visitados por la hormiga
        self.visitados.append(inicio)
        # Situamos el puntero del nodo actual en el nodo de inicio
        self.actual = inicio
        # Eliminamo el nodo inicial de la lista de permitidos
        self.permitidos.remove(inicio)
    
    # Selecciona un nodo
    # Para calcular el posible camino a tomar de una hormiga, se utiliza la Ecuación:
    # Pxy = (τ(e)*N(e))/(Στ(e)*N(e)), donde τ(e) son las feromonas depositadas por ese
    # camino (A a B), N(e) es la visibilidad del camino dado por N(e) = 1/peso_del_camino,
    # y el denominador de esta ecuación siendo la suma del cálculo anterior (τ(e)*N(e)) de
    # todos los posibles que tiene esa hormiga
    def selecciona_nodo(self):
        N = 1/self.grafo.distances
        numerador = []
        denominador = 0
        for i in self.permitidos:
            aux = (self.grafo.feromonas[self.actual][i]**self.grafo.alpha)*(N[self.actual][i]**self.grafo.beta)
            denominador += aux
            numerador.append(aux)
        P = numerador/denominador
        peso_minimo = -1
        nodo_minimo = -1
        for i,p in enumerate(P):
            if p>peso_minimo:
                peso_minimo = p
                nodo_minimo = self.permitidos[i]
        
        self.permitidos.remove(nodo_minimo)
        self.visitados.append(nodo_minimo)
        self.distancia_total += self.grafo.distances[self.actual][nodo_minimo]
        self.actual = nodo_minimo
    
    # Calcula las Δτxy como Q/Lk si se usó xy ó 0 en caso contrario
    def delta_feromonas(self):
        self.feromonas = [[0 for j in range(self.grafo.npuntos)] for i in range(self.grafo.npuntos)]
        for i in range(0,len(self.visitados)-1):
            nodo = self.visitados[i]
            siguiente = self.visitados[i+1]
            # Q es el parámetro de aprendizaje y Lk (distancia total) es el costo del camino por hormiga
            self.feromonas[nodo][siguiente] = self.grafo.Q/self.distancia_total
            
def buscar(d,m,r,Q,a,b):

    colonia = Colonia(m,r)
    grafo = Grafo(d,m,a,b,Q)
    camino,longitud = colonia.resolver(grafo)
    return camino,longitud

positions, fname = leer_fichero("./pbk411.tsp")
for position in positions:
    plt.plot(position[0],position[1], marker="o", color="red")
plt.show()
distances = get_distances(positions)

# Parámetros mapa 411

m = 10
r = 0.5
q = 0.6
a = 2
b = 20

# Parámetros mapa 29
"""
m = 20
r = 0.25
q = 0.6
a = 0.35
b = 1
"""
mejor = np.inf
mejorc = []
media = 0
size = 5
for i in range(size):
    d = distances.copy()

    camino,longitud = buscar(d,m,r,q,a,b)
    
    if longitud<mejor:
        mejorc = camino
        mejor = longitud
    media += longitud
        
longitud_media = media/size

pintar(positions,mejorc)
print('Longitud_media:',longitud_media)
print('Mejor_longitud:',mejor)


# m = 20, ρ = 0,25, α = 0,35, Q = 0,6, β = 1.





