import pyg4ometry
import json
import itertools
import pyg4ometry.geant4 as _g4
import numpy as np
import pickle
import periodictable
from materialcreator import MaterialCreator


def pickleit(obj, path):
    outfile = open(path, 'wb')
    pickle.dump(obj, outfile)
    outfile.close()

def unpickleit(path):
    infile = open(path, 'rb')
    output = pickle.load(infile)
    infile.close()
    return output

class SingleBin:
    def __init__(self, bin_name, reg, world, material_creator):
        self.bin_name = bin_name
        self.reg = reg
        self.world = world
        self.material_creator = material_creator
    
    def set_volume(self, volume):
        self.volume = volume
    
    def set_materialZ(self, materialZ):
        # "ZAF": "Z1A1.00794F0.043;Z6A12.0107F0.793;Z7A14.0068F0.007;Z8A15.9994F0.121;Z14A28.0854F0.036;",
        self.materialZ = materialZ
        self.material = _g4.MaterialPredefined("G4_Fe")
    
    def set_materialZAF(self, material, ZAFstring):
        self.material = None
        if len(ZAFstring) == 0:
            return
        self.material = self.material_creator.create_material(ZAFstring, material[4])
        # molarmassdensity = material[4]
        # elements_string = ZAFstring.split(";")
        # elements = []
        # for each_element_string in elements_string:
        #     if len(each_element_string) == 0:
        #         continue
        #     each_element_string = each_element_string.replace("Z", "")
        #     each_z = float(each_element_string.split("A")[0])
        #     each_a = float(each_element_string.split("A")[1].split("F")[0])
        #     each_f = float(each_element_string.split("F")[1])
        #     elements.append((each_z, each_a, each_f))
        # molarmass = sum([each[1] * each[2] for each in elements])
        # massdensity = molarmassdensity * 1000 * 1000 # g/m3

        # self.material = _g4.MaterialCompound("material_" + self.bin_name, massdensity, len(elements), self.reg)
        # for i in range(len(elements)):
        #     each_z, each_a, each_f = elements[i]
        #     each_symble = "NONE"
        #     each_name = "NONE"
        #     for el in periodictable.elements:
        #         if el.number == int(each_z):
        #             each_symble = el.symbol
        #             each_name = el.name
        #             break
        #     element = _g4.ElementSimple(each_name, each_symble, int(each_z), each_a)
        #     self.material.add_element_massfraction(element, each_f*each_a/molarmass)




        # ne = _g4.ElementSimple("nitrogen","N",7,14.01)
        # oe = _g4.ElementSimple("oxygen","O",8,16.0)
        # wm.add_element_massfraction(ne,0.7)
        # wm.add_element_massfraction(oe,0.3)

        # print(ZAFstring)
        # # exit(1)
        # if isinstance(material, str):
        #     self.material = _g4.MaterialPredefined(material)
        # else:
        #     self.material = material
    
    def set_PhiRange(self, phimin, phimax):
        self.phimin = phimin
        self.phimax = phimax
        self.dphi = abs(phimax - phimin)
    
    def set_ZRange(self, zmin, zmax):
        self.zmin = zmin
        self.zmax = zmax
        self.dz = abs(zmax - zmin)
    
    def set_Rrange(self, rmin, rmax):
        self.rmin = rmin
        self.rmax = rmax
        self.dr = abs(rmax - rmin)
    
    def set_translation(self, x, y, z):
        self.translation = [x, y, z]

    def set_rotation(self, x, y, z):
        self.rotation = [x, y, z]
    
    def create_volume(self):
        if self.dz == 0 or self.dr == 0 or self.dphi == 0 or self.material is None:
            return
        if self.rmin < 0.0001:
            self.rmin = 0.0001
        self.solid = pyg4ometry.geant4.solid.Tubs(f"{self.bin_name}_Solid_Tubs", self.rmin, self.rmax, self.dz, self.phimin, self.dphi - 0.001, self.reg)
        self.logical = pyg4ometry.geant4.LogicalVolume(self.solid, self.material, f"{self.bin_name}_Logical", self.reg)
        self.physical = pyg4ometry.geant4.PhysicalVolume(self.rotation, self.translation,self.logical, f"{self.bin_name}_Physical", self.world, self.reg)



class DetectorBinCollection:
    def __init__(self, name, reg, world, jsonmapping, material_creator):
        self.name = name
        self.reg = reg
        self.world = world
        self.jsonmapping = jsonmapping
        self.material_creator = material_creator
        self.bins = []
        self.get_meta_physical()
        self.get_meta_binning()
        self.generate_allbins()
    
    def get_meta_physical(self):
        meta_map = self.jsonmapping["binUtility"]
        self.translation = None
        self.rotation = None
        if "transform" in meta_map:
            self.translation = meta_map["transform"]["translation"]
            self.rotation = meta_map["transform"]["rotation"]

        if self.translation is None:
            self.translation = [0, 0, 0]
        if self.rotation is None:
            self.rotation = [0, 0, 0]

    def get_meta_binning(self):
        meta_map = self.jsonmapping["binUtility"]
        self.binPhi = []
        self.binR = []
        self.binZ = []
        self.binorder = []
        self.bin_collection = []
        for each_binmetamap in meta_map["binningdata"]:
            each_binning_tem = np.linspace(each_binmetamap["min"], each_binmetamap["max"], each_binmetamap["bins"]+1)
            each_binning = [[each_binning_tem[i], each_binning_tem[i+1]] for i in range(len(each_binning_tem)-1)]
            if each_binmetamap["value"] == "binR" or each_binmetamap["value"] == "AxisR":
                self.binR = each_binning
                self.binorder.append(["binR", self.binR])
            elif each_binmetamap["value"] == "binPhi" or each_binmetamap["value"] == "AxisPhi":
                self.binPhi = each_binning
                self.binorder.append(["binPhi", self.binPhi])
            elif each_binmetamap["value"] == "binZ" or each_binmetamap["value"] == "AxisZ":
                self.binZ = each_binning
                self.binorder.append(["binZ", self.binZ])
            else:
                raise ValueError(f"Unknown binning type {each_binmetamap['value']}")


    def process_singlebin(self, name, binR, binPhi, binZ, material, each_material_discription_ZAF=None):
        singlebin = SingleBin(f"{self.name}_{name}", self.reg, self.world, self.material_creator)
        singlebin.set_Rrange(binR[0], binR[1])
        singlebin.set_ZRange(binZ[0], binZ[1])
        singlebin.set_PhiRange(binPhi[0], binPhi[1])
        singlebin.set_materialZAF(material, each_material_discription_ZAF) # TODO: Change this to actual material
        singlebin.set_translation(self.translation[0], self.translation[1], self.translation[2])
        singlebin.set_rotation(self.rotation[0], self.rotation[1], self.rotation[2])
        singlebin.create_volume()
        return singlebin


    def generate_allbins(self):
        material_map = self.jsonmapping["data"]
        if len(material_map[0]) != len(self.binorder[0][1]):
            raise ValueError("Number of bins in binningdata and material data do not match")
        if len(material_map) != len(self.binorder[1][1]):
            raise ValueError("Number of bins in binningdata and material data do not match")

        for i in range(len(material_map)):
            each_bins = []
            for j in range(len(material_map[0])):
                binname = f"bin_{i}_{j}"
                each_material_map = material_map[i][j]
                each_material_discription = each_material_map["material"]
                each_material_discription_ZAF = each_material_map["ZAF"]
                each_thickness = float(each_material_map["thickness"])
                binPhi = binZ = [-each_thickness/2, each_thickness/2]
                binR = [0, each_thickness]
                if self.binorder[0][0] == "binR":
                    binR = self.binorder[0][1][j]
                if self.binorder[0][0] == "binPhi":
                    binPhi = self.binorder[0][1][j]
                if self.binorder[0][0] == "binZ":
                    binZ = self.binorder[0][1][j]
                if self.binorder[1][0] == "binR":
                    binR = self.binorder[1][1][i]
                if self.binorder[1][0] == "binPhi":
                    binPhi = self.binorder[1][1][i]
                if self.binorder[1][0] == "binZ":
                    binZ = self.binorder[1][1][i]
                eachbin = self.process_singlebin(binname, binR, binPhi, binZ, each_material_discription, each_material_discription_ZAF)
                each_bins.append(eachbin)
            self.bin_collection.append(each_bins)

class Detector:
    def __init__(self, name, jsonmapping):
        self.name = name
        self.jsonmapping = jsonmapping
        self.reg  = pyg4ometry.geant4.Registry()
        self.world_solid   = pyg4ometry.geant4.solid.Box("worldsolid", 10000, 10000, 10000, self.reg)
        self.world_logical   = pyg4ometry.geant4.LogicalVolume(self.world_solid,"G4_Galactic","worldlogical", self.reg)
        self.reg.setWorld(self.world_logical.name)
        self.material_creator = MaterialCreator(self.reg)

        self.generate_layers()
    
    def check_layer(self, layer):
        json_bins = layer["value"]["material"]["data"] 
        for each_bin in list(itertools.chain(*json_bins)):
            if float(each_bin["thickness"]) > 0:
                return True
        return False

    def generate_layers(self):
        print("Generating detector")
        self.layers = {}
        i = 0
        count = 0
        for each_layer in self.jsonmapping:
            name = f"dic_{i}_"
            if "layer" in each_layer:
                name += f"layer{each_layer['layer']}_"
                if each_layer['layer'] != 2:
                    continue
            if "approach" in each_layer:
                name += f"approach{each_layer['approach']}_"
                print(f"Skipping approach layer {name}")
                continue
            if "boundary" in each_layer:
                name += f"boundary{each_layer['boundary']}_"


            name = name[:-1]
            if self.check_layer(each_layer):
                if name in self.layers:
                    print(f"Layer {name} already exists")
                count += 1
                if count == 2 or True:
                    print(f"Creating layer {name}")
                    each_bin_collection = DetectorBinCollection(name, self.reg, self.world_logical, each_layer["value"]["material"], self.material_creator)
                    self.layers[name] = each_bin_collection
                # each_bin_collection = DetectorBinCollection(name, self.reg, self.world_logical, each_layer["value"]["material"])
                # self.layers[name] = each_bin_collection
            else:
                print(f"Skipping layer {name}")
            i += 1

        print("Detector generation complete")

    def save(self, filename):
        print("Saving detector as gdml")
        w = pyg4ometry.gdml.Writer()
        w.addDetector(self.reg)
        filename = filename.replace(".gdml", "")
        w.write(f"{filename}.gdml")

    def visualize(self):
        print("Visualizing detector")
        v = pyg4ometry.visualisation.VtkViewer()
        v.addLogicalVolume(self.world_logical)
        v.addAxes(20)
        v.view()

def main():
    # Load the JSON file
    with open("material-map.json") as f:
        atlas_geometry = json.load(f)
    atlas_geometry = atlas_geometry["Surfaces"]["entries"]
    atlas_detector = Detector("ATLAS", atlas_geometry)

    atlas_detector.save("ATLAS.gdml")
    # atlas_detector.visualize()


    # print(atlas_geometry)

    # reg  = pyg4ometry.geant4.Registry()
    # world_solid   = pyg4ometry.geant4.solid.Box("worldsolid", 5000, 5000, 5000, reg)
    # world_logical   = pyg4ometry.geant4.LogicalVolume(world_solid,"G4_Galactic","wroldlogical",reg)
    # reg.setWorld(world_logical.name)

    # testbin = SingleBin("testbin", reg, world_logical)
    # testbin.set_Rrange(700, 1000)
    # testbin.set_ZRange(-100, 100)
    # testbin.set_PhiRange(0, 3.14159)
    # testbin.set_materialZ(26)
    # testbin.set_translation(0, 0, 0)
    # testbin.set_rotation(0, 0, 0)
    # testbin.create_volume()

    # w = pyg4ometry.gdml.Writer()
    # w.addDetector(reg)
    # w.write('file.gdml')

    # v = pyg4ometry.visualisation.VtkViewer()
    # v.addLogicalVolume(world_logical)
    # v.addAxes(20)
    # v.view()

if __name__ == "__main__":
    main()
