class Protein:
    def __init__(self, name= "", mutation = "", species = "HUMAN", type = "", value = None, dbSource = "", experimentalValue = False):
        self.name = name
        self.mutation = mutation
        self.species = species
        self.type = type
        self.value = value
        self.dbSource = dbSource
        self.experimentalValue = experimentalValue

    def __str__(self) -> str:
        return "Protein: " + str(self.name) + " Species: " + str(self.species) + " Value: " + str(self.value)

    def loadJSON(self, obj):
        if isinstance(obj, dict):
            self.name = obj["name"]
            self.mutation = obj["mutation"]
            self.species = obj["species"]
            self.type = obj["type"]
            self.value = obj["value"]
            self.dbSource = obj["dbSource"]
            self.experimentalValue = obj["experimentalValue"]
        else:
            print("The given Protein is not a dictionary and can not be parsed!")