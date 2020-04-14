class Protein:
    def __init__(self, name="", mutation="", species="HUMAN", type="", value=0.0, db_source="",
                 experimental_value=False, less_greater_than_sign=0):
        self.name = name
        self.mutation = mutation
        self.species = species
        self.type = type
        self.value = value
        self.db_source = db_source
        self.experimental_value = experimental_value
        self.less_greater_than_sign = less_greater_than_sign

    def __str__(self) -> str:
        return "Protein: " + str(self.name) + " Species: " + str(self.species) + " Value: " + str(self.value)

    def load_json(self, obj):
        if isinstance(obj, dict):
            self.name = obj["name"]
            self.mutation = obj["mutation"]
            self.species = obj["species"]
            self.type = obj["type"]
            self.value = obj["value"]
            self.db_source = obj["db_source"]
            self.experimental_value = obj["experimental_value"]
            self.less_greater_than_sign = obj["less_greater_than_sign"]
        else:
            print("The given Protein is not a dictionary and can not be parsed!")
