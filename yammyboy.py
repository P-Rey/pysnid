import yaml

def yaml_loader(filepath):
    with open(filepath, "r") as file_descriptor:
        data = yaml.load(file_descriptor)
    return data

def yaml_dump(filepath, data):
    with open(filepath, "w") as file_descriptor:
        yaml.dump(data, file_descriptor)

if __name__ == "__main__":
    filepath = "example.yml"
    data = yaml_loader(filepath)
    print(data)
    path = "supernova_"+str(1)
    items = data.get(path)
    lhs,rhs = items.split("type:",1)
#    for item_name, item_value in items.iteritems():
#         print(item_name, item_value)