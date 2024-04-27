# -*- coding: utf-8 -*-
"""
Modul, welches die ganze Logik der Punkte abkapselt
"""
import csv
import matplotlib.pyplot as plt
import random
import math

class Point:
    def __init__(self, coords, isCenter = False, clusterAffiliation = -1):
        # Dimensionalität
        self.dim = len(coords)
        
        # liste von Coordinaten
        # Speichern als Tupel für Performance
        self.coordinates = coords
                
        # Zentrum
        self.isCenter = isCenter
        
        # Clusterzugehörigkeit
        self.cluster = clusterAffiliation
        
        
    # To-String-Methode für schönes anschauen
    def toStringLong(self):
        stringRep = "Dim: " + str(self.dim)
        stringRep +=";\t Center: " + str(self.isCenter)
        stringRep +=";\t Cluster: "  + str(self.cluster)  
        stringRep +=";\t Coords: " + str(self.coordinates)
        """
        
        +";\tCoords: ["
        for i in range(0, self.dim - 1):
            stringRep += str(self.coordinates[i]) + ", "
        stringRep += str(self.coordinates[self.dim - 1]) + "]"
        """
        return stringRep
    
    # To-String-Metode für die Datei später
    # Format: dim,color,coord1,coord2,...,coordn
    def toString(self):
        stringRep = str(self.dim) + ","
        for i in range(0, self.dim - 1):
            stringRep += str(self.coordinates[i]) + ","
        stringRep += str(self.coordinates[self.dim - 1])
        return stringRep
    
    
    # Euclidische Distanz zwischen zwei Punkten
    def dist_to(self, other):
        
        # Sanaty check: Dimensionen Stimmen
        if(len(self.coordinates) != len(other.coordinates)):
            print(f"Dimensions not matching: {len(self.coordinates)} != {len(other.coordinates)}")
            # Fehlercode
            return -1
        
        # Initialize distance to zero
        distance = 0
        
        # Loop through each coordinate in the point2
        for coord1, coord2 in zip(self.coordinates, other.coordinates):
            # Add the square of the difference of coordinates to the distance
            distance += (coord2 - coord1) ** 2
        
        # Take the square root of the sum to get the Euclidean distance
        distance = math.sqrt(distance)
        
        return distance



def show_points(listOfPoints, showRadii = False, radii = None):
    construct_plot(listOfPoints, showRadii, radii).show()
    return

def save_plot(fileName, listOfPoints, showRadii = False, radii = None):
    construct_plot(listOfPoints, showRadii, radii).savefig(fileName)
    return


def construct_plot(listOfPoints, showRadii = False, radii = None):
    fig, ax = plt.subplots(figsize=(5,5))
    
    # Ersten zwei Coordinaten
    firstCoordinate = [point.coordinates[0] for point in listOfPoints]
    secondCoordinate = [point.coordinates[1] for point in listOfPoints]
    
    # Einfärben
    colors = [point.cluster for point in listOfPoints]

 
    ax.scatter(firstCoordinate, secondCoordinate, alpha=0.5, c=colors)
    
    # Zentren Markieren
    for i in range(0, len(listOfPoints)):
        if(listOfPoints[i].isCenter):
            ax.scatter(listOfPoints[i].coordinates[0], listOfPoints[i].coordinates[1], marker="+", s=250, c="black")
            # ggf. Radius hinzufügen
            if(showRadii):
                circle = plt.Circle((listOfPoints[i].coordinates[0], listOfPoints[i].coordinates[1]),
                                    color="black", fill=False)
                circle.set_radius(radii[listOfPoints[i].cluster])
                ax.add_artist(circle)
            
    
    return plt




def read_points(fileName, withMeta = False):
    csvFile = open(fileName, "r")
    rawDataList = list(csv.reader(csvFile, delimiter=","))
    
    # Falls Methadaten schon vorhanden (center und Cluster)
    if (withMeta):
        return parse_points_from_list_with_meta(rawDataList)
    
    pointsData = parse_points_from_list(rawDataList)
    
    return pointsData



def parse_points_from_list_with_meta(rawDataList):
    # Wie viele Punkte werden eingelesen
    numberOfPoints = len(rawDataList)
    
    # Leere Liste für die Punkte als Rückgabe
    listOfPoints = []
    
    # Starte bei 1 da erste Zeile die Überschrifften enthält
    for i in range(1, numberOfPoints):
        # Metadata
        dim = int(rawDataList[i][0])
        center = bool(int(rawDataList[i][1]))
        cluster = int(rawDataList[i][2])
        
        # Coordinates
        coords = []
        for k in range(3, 3+dim):
            coords.append(float(rawDataList[i][k]))
            
        # Liste in Tupel umwandeln
        coordTupel = tuple(coords)
        
        # Point via Constructor
        newPoint = Point(dim, center, cluster, coordTupel)
        
        # add Point to list
        listOfPoints.append(newPoint)
    
    
    return listOfPoints



def parse_points_from_list(rawDataList):
    # Wie viele Punkte werden eingelesen
    numberOfPoints = len(rawDataList)
    
    # Leere Liste für die Punkte als Rückgabe
    listOfPoints = []
    
    # Starte bei 1 da erste Zeile die Überschrifften enthält
    for i in range(1, numberOfPoints):
        # Metadata
        dim = int(rawDataList[i][0])
        
        # Coordinates
        coords = []
        for k in range(1, 1+dim):
            coords.append(float(rawDataList[i][k]))
            
        # Liste in Tupel umwandeln
        coordTupel = tuple(coords)
        
        # Point via Constructor
        newPoint = Point(coordTupel)
        
        # add Point to list
        listOfPoints.append(newPoint)
    
    
    return listOfPoints



def create_random_point(dim = 2, borderLow = 0, borderHigh = 1):
    
    # Zufällige Koordinaten bestimmen
    coords = []
    
    for i in range(0, dim):
        coords.append((random.random() * (borderHigh - borderLow)) + borderLow)
        
    return Point(coords)

