import json
from Protein import Protein


class Drug:
    def __init__(self, name="", status="", tar=None, enz=None, tra=None, car=None, atc=None):
        self.name = name
        self.status = status
        self.tar = tar if tar is not None else {}
        self.enz = enz if enz is not None else {}
        self.tra = tra if tra is not None else {}
        self.car = car if car is not None else {}
        self.atc = atc if atc is not None else []

    def __str__(self) -> str:
        return "Drug: " + str(self.name) + "\nStatus" + str(self.status)

    list_of_protein_types = ["tar", "enz", "tra", "car"]

    def has_activity_value(self):
        if len(self.tar) > 0 or len(self.enz) > 0 or len(self.tra) > 0 or len(self.car):
            return True
        return False

    def has_target_and_activity_value(self):
        return True if self.tar else False

    def load_json(self, obj):
        if isinstance(obj, dict):
            self.name = obj["name"]
            self.status = obj["status"]
            self.atc = obj["atc"]
            for type in self.list_of_protein_types:
                type_dict = {}
                for k, p in obj[type].items():
                    prt = Protein()
                    prt.load_json(p)
                    type_dict[k] = prt
                obj[type] = type_dict
            self.tar = obj["tar"]
            self.enz = obj["enz"]
            self.tra = obj["tra"]
            self.car = obj["car"]

        else:
            print("The given drug is not a dictionary and can not be parsed!")
