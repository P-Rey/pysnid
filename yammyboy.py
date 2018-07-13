import yaml

def yaml_loader(filepath):
    with open(filepath, "r") as file_descriptor:
        data = yaml.load(file_descriptor)
    return data

def yaml_dump(filepath, data):
    with open(filepath, "w") as file_descriptor:
        yaml.dump(data, file_descriptor)

if __name__ == "__main__":
    filepath ="template_directory.yml"
    data = yaml_loader(filepath)

#    suffix = str(i)
    
    path = "supernova_"+str(1)
    items = data.get(path)
    junk, rhs= items.split("label:", 1)
    sn_type, junk = rhs.split("type:", 1)
    del junk
#    fuck, this = lhs.split(":",1)
#    lhs1,rhs1 = rhs.split("subtype:",1)
