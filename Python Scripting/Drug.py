import json
from Protein import Protein


class Drug:
    def __init__(self, name="", status="", tar={}, enz={}, tra={}, car={}, atc=[]):
        self.name = name
        self.status = status
        self.tar = tar
        self.enz = enz
        self.tra = tra
        self.car = car
        self.atc = atc

    def __str__(self) -> str:
        return "Drug: " + str(self.name) + "\nStatus" + str(self.status)

    listOfProteins = ["tar", "enz", "tra", "car"]

    def hasActivityValue(self):
        if len(self.tar) > 0 or len(self.enz) > 0 or len(self.tra) > 0 or len(self.car):
            return True
        return False

    def hasTargetAndActivityValue(self):
        return True if self.tar else False

    def loadJSON(self, obj):
        if isinstance(obj, dict):
            self.name = obj["name"]
            self.status = obj["status"]
            self.atc = obj["atc"]
            for type in self.listOfProteins:
                tempList = {}
                for k, p in obj[type].items():
                    prt = Protein()
                    prt.loadJSON(p)
                    tempList[k] = prt
                obj[type] = tempList
            self.tar = obj["tar"]
            self.enz = obj["enz"]
            self.tra = obj["tra"]
            self.car = obj["car"]

        else:
            print("The given Drug is not a dictionary and can not be parsed!")
