import config
from Drug import Drug
import json
import numpy as np
import sys


def get_all_drugs(path):
    drugs = {}
    try:
        file = open(path, "r")
        json_dict = json.loads(file.read())
        for key, value in json_dict.items():
            drug = Drug()
            drug.load_json(value)
            drugs[drug.name] = drug

        file.close()
    except:
        print("There was a problem reading the matrix from: " + str(path))
        print("Unexpected error:", sys.exc_info()[0])
    return drugs


def load_from_files(drugs_list_path, targets_list_path):
    all_drugs = []
    all_targets = []
    try:
        with open(drugs_list_path) as json_file:
            all_drugs = json.load(json_file)
        with open(targets_list_path) as json_file:
            all_targets = json.load(json_file)
    except:
        print("There was a problem reading the drugs and targets from: " + str(drugs_list_path) + "\nand: " + str(
            targets_list_path))
        print("Unexpected error:", sys.exc_info()[0])
    return all_drugs, all_targets


def load_matrix(path):
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


def save_to_file(path, list_of_things):
    try:
        list_of_things_json = json.dumps(list_of_things)
        file = open(path, "w")
        file.write(list_of_things_json)
        file.close()
    except:
        print("There was a problem writing to the file: " + str(path))
        print("Unexpected error:", sys.exc_info()[0])


def save_drugs_to_files(path, drugs):
    try:
        file = open(path, "w")
        json_string = json.dumps(drugs, default=lambda o: o.__dict__, sort_keys=False, indent=4)
        file.write(json_string)
        file.close()
    except:
        print("There was a problem writing to the file: " + str(path))
        print("Unexpected error:", sys.exc_info()[0])


def load_go_dictionary(path):
    dict_go = {}
    try:
        with open(path) as json_file:
            dict_loaded = json.load(json_file)

        for key, value in dict_loaded.items():
            dict_go[key] = value[:-1].split(",")

    except:
        print("There was a problem reading from the file: " + str(path))
        print("Unexpected error:", sys.exc_info()[0])

    return dict_go


def load_and_build_dictionary(path, position, split_character):
    new_dict = {}
    try:
        file = open(path, "r")
        for line in file:
            line_split = line.split(split_character)
            new_dict[line_split[0]] = line_split[position]

        file.close()
    except:
        print("There was a problem writing to the file: " + str(path))
        print("Unexpected error:", sys.exc_info()[0])

    return new_dict
