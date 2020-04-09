import config
from Drug import Drug
import json
import numpy as np
import sys


def getAllDrugsAsDict(path):
    allDrugs = {}
    drugs = getAllDrugs(path)
    for drug in drugs:
        allDrugs[drug.name] = drug

    return allDrugs

def getAllDrugs(path):
    drugs = {}
    try:
        file = open(path, "r")
        jsonDict = json.loads(file.read())
        for key, value in jsonDict.items():
            drug = Drug()
            drug.loadJSON(value)
            drugs[drug.name] = drug

        file.close()
    except:
        print("There was a problem reading the matrix from: " + str(path))
        print("Unexpected error:", sys.exc_info()[0])
    return drugs


def loadFromFiles(drugsListPath, targetsListPath):
    allDrugs = []
    allTargets = []
    try:
        with open(drugsListPath) as json_file:
            allDrugs = json.load(json_file)
        with open(targetsListPath) as json_file:
            allTargets = json.load(json_file)
    except:
        print("There was a problem reading the Drugs and Targets from: " + str(drugsListPath) + "\nand: " + str(
            targetsListPath))
        print("Unexpected error:", sys.exc_info()[0])
    return allDrugs, allTargets


def loadMatrix(path):
    matrix = []
    try:
        matrix = np.load(path)
    except IOError as e:
        print("There was a problem reading the matrix from: " + str(path))
        print("Error: " + str(e))
    except:
        print("There was a problem reading the matrix from: " + str(path))
        print("Unexpected error:", sys.exc_info()[0])
    return matrix


def saveToFile(path, listOfThings):
    try:
        listOfThingsJSON = json.dumps(listOfThings)
        file = open(path, "w")
        file.write(listOfThingsJSON)
        file.close()
    except:
        print("There was a problem writing to the file: " + str(path))
        print("Unexpected error:", sys.exc_info()[0])


def saveDrugsToFiles(path, drugs):
    try:
        file = open(path, "w")
        jsonString = json.dumps(drugs, default=lambda o: o.__dict__, sort_keys=False, indent=4)
        file.write(jsonString)
        file.close()
    except:
        print("There was a problem writing to the file: " + str(path))
        print("Unexpected error:", sys.exc_info()[0])


def loadGODictionary(path):
    dictGO = {}
    try:
        with open(path) as json_file:
            dictLoaded = json.load(json_file)

        for key, value in dictLoaded.items():
            dictGO[key] = value[:-1].split(",")

    except:
        print("There was a problem reading from the file: " + str(path))
        print("Unexpected error:", sys.exc_info()[0])

    return dictGO

def loadAndBuildDictionary(path, position, splitCharacter):
    newDict = {}
    try:
        file = open(path, "r")
        for line in file:
            lineSplit = line.split(splitCharacter)
            newDict[lineSplit[0]] = lineSplit[position]

        file.close()
    except:
        print("There was a problem writing to the file: " + str(path))
        print("Unexpected error:", sys.exc_info()[0])

    return newDict
