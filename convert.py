import pyg4ometry
import json
import itertools
import pyg4ometry.geant4 as _g4
import numpy as np
import pickle
import periodictable
from materialcreator import MaterialCreator
from solidcreator import SolidCreator

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
    def __init__(self, bin_name, reg, world, material_creator, solid_creator):
        self.bin_name = bin_name
        self.reg = reg
        self.world = world
        self.material_creator = material_creator
        self.solid_creator = solid_creator

    def set_volume(self, volume):
        self.volume = volume
    
    def set_materialZ(self, materialZ):
        self.materialZ = materialZ
        self.material = _g4.MaterialPredefined("G4_Fe")
    
    def set_materialZAF(self, material, ZAFstring):
        self.material = None
        if len(ZAFstring) == 0:
            return
        self.material = self.material_creator.create_material(ZAFstring, material[4] * material[2] * 1000) # convert from /mm3 to /cm3
    
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
        self.solid = self.solid_creator.create_tub(self.rmin+ 0.00001, self.rmax- 0.00001, self.dz- 0.00001, self.dphi- 0.00001)
        #print(self.bin_name, self.translation, self.zmax, self.zmin)
        self.translation[2] += (self.zmax + self.zmin)/2
        self.set_rotation(0, 0, self.phimin)
        # self.set_translation(0, 0, (self.zmin + self.zmax)/2)
        # self.set_translation(0, 0, 0)
        # self.solid = pyg4ometry.geant4.solid.Tubs(f"{self.bin_name}_Solid_Tubs", self.rmin, self.rmax, self.dz, self.phimin, self.dphi - 0.001, self.reg)
        self.logical = pyg4ometry.geant4.LogicalVolume(self.solid, self.material, f"{self.bin_name}_Logical", self.reg)
        self.physical = pyg4ometry.geant4.PhysicalVolume(self.rotation, self.translation, self.logical, f"{self.bin_name}_Physical", self.world, self.reg)



class DetectorBinCollection:
    def __init__(self, name, reg, world, materialmap, material_creator, solid_creator, r):
        self.name = name
        self.reg = reg
        self.world = world
        self.materialmap = materialmap
        self.material_creator = material_creator
        self.solid_creator = solid_creator
        self.r = r
        self.bins = []
        self.get_meta_physical()
        self.get_meta_binning()
        self.generate_allbins()
    
    def get_meta_physical(self):
        meta_map = self.materialmap["binUtility"]
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
        meta_map = self.materialmap["binUtility"]
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
        singlebin = SingleBin(f"{self.name}_{name}", self.reg, self.world, self.material_creator, self.solid_creator)
        singlebin.set_Rrange(binR[0], binR[1])
        singlebin.set_ZRange(binZ[0], binZ[1])
        singlebin.set_PhiRange(binPhi[0], binPhi[1])
        singlebin.set_materialZAF(material, each_material_discription_ZAF)
        singlebin.set_translation(self.translation[0], self.translation[1], self.translation[2])
        singlebin.set_rotation(self.rotation[0], self.rotation[1], self.rotation[2])
        singlebin.create_volume()
        return singlebin


    def generate_allbins(self):
        material_map = self.materialmap["data"]
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
                each_thickness = 1
                binPhi = binZ = [-each_thickness/2, each_thickness/2]
                binR = [0, each_thickness]
                isrepresentative = False
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
                if self.binorder[0][0] in ("binPhi", "binZ") and self.binorder[1][0] in ("binPhi", "binZ"):
                    binR = [self.r - each_thickness/2, self.r + each_thickness/2]
                if not isrepresentative:
                    eachbin = self.process_singlebin(binname, binR, binPhi, binZ, each_material_discription, each_material_discription_ZAF)
                    each_bins.append(eachbin)
            self.bin_collection.append(each_bins)

class Detector:
    def __init__(self, name, materialmap, geometryR):
        self.name = name
        self.materialmap = materialmap
        self.geometryR = geometryR
        self.reg  = pyg4ometry.geant4.Registry()
        self.material_creator = MaterialCreator(self.reg)
        self.solid_creator = SolidCreator(self.reg, tolerance=0)

        self.new_logicals = []
        self.new_physicals = []

        #self.generate_layers()
    
    def create_world(self):
        self.world_solid = pyg4ometry.geant4.solid.Box("worldsolid", 50000, 50000, 100000, self.reg)
        self.create_air()
        self.world_logical = pyg4ometry.geant4.LogicalVolume(self.world_solid, self.air, "worldlogical", self.reg)
        self.reg.setWorld(self.world_logical.name)

    def use_world(self, logical):
        self.world_solid = logical.solid
        self.reg.transferSolid(logical.solid)
        self.reg.transferMaterial(logical.material)
        self.world_logical = pyg4ometry.geant4.LogicalVolume(self.world_solid, logical.material, "worldlogical", self.reg)
        self.reg.setWorld(self.world_logical.name)

    def add_volume(self, physical):
        rotation = [physical.rotation.x.eval(), physical.rotation.y.eval(), physical.rotation.z.eval(), physical.rotation.unit]
        position = [physical.position.x.eval(), physical.position.y.eval(), physical.position.z.eval(), physical.position.unit]
        logical = physical.logicalVolume
        self.reg.transferSolid(logical.solid)
        self.reg.transferMaterial(logical.material)
        self.reg.transferLogicalVolume(logical)
        #self.reg.transferPhysicalVolume(physical)
        new_physical = pyg4ometry.geant4.PhysicalVolume(rotation, position, logical, physical.name, self.reg.getWorldVolume(), self.reg)
        # self.new_logicals.append(logical)
        # self.new_physicals.append(new_physical)

    def create_air(self):
        self.c12 = _g4.Isotope("C12", 6, 12, 12)
        self.c13 = _g4.Isotope("C13", 6, 13, 13.0034)
        self.carbon = _g4.ElementIsotopeMixture("worldcarbon", "C", 2)
        self.carbon.add_isotope(self.c12, 0.9893)
        self.carbon.add_isotope(self.c13, 0.0107)

        self.n14 = _g4.Isotope("N14", 7, 14, 14.0031)
        self.n15 = _g4.Isotope("N15", 7, 15, 15.0001)
        self.nitrogen = _g4.ElementIsotopeMixture("worldnitrogen", "N", 2)
        self.nitrogen.add_isotope(self.n14, 0.99632)
        self.nitrogen.add_isotope(self.n15, 0.00368)

        self.o16 = _g4.Isotope("O16", 8, 16, 15.9949)
        self.o17 = _g4.Isotope("O17", 8, 17, 16.9991)
        self.o18 = _g4.Isotope("O18", 8, 18, 17.9992)
        self.oxygen = _g4.ElementIsotopeMixture("worldoxygen", "O", 3)
        self.oxygen.add_isotope(self.o16, 0.99757)
        self.oxygen.add_isotope(self.o17, 0.00038)
        self.oxygen.add_isotope(self.o18, 0.00205)

        self.ar36 = _g4.Isotope("Ar36", 18, 36, 35.9675)
        self.ar38 = _g4.Isotope("Ar38", 18, 38, 37.9627)
        self.ar40 = _g4.Isotope("Ar40", 18, 40, 39.9624)
        self.argon = _g4.ElementIsotopeMixture("worldargon", "Ar", 3)
        self.argon.add_isotope(self.ar36, 0.003365)
        self.argon.add_isotope(self.ar38, 0.000632)
        self.argon.add_isotope(self.ar40, 0.996003)

        self.air = _g4.MaterialCompound("worldair", 0.00120479, 4, self.reg)
        self.air.add_element_massfraction(self.carbon, 0.000124000124000124)
        self.air.add_element_massfraction(self.nitrogen, 0.755267755267755)
        self.air.add_element_massfraction(self.oxygen, 0.231781231781232)
        self.air.add_element_massfraction(self.argon, 0.0128270128270128)

        self.air.set_temperature(293.15)

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
        for each_layer in self.materialmap:
            name = f"dic_{i}_"
            if "volume" in each_layer:
                name += f"volume{each_layer['volume']}_"
            # if each_layer["volume"] != 2:
            #     continue

            if "layer" in each_layer:
                name += f"layer{each_layer['layer']}_"

                # if each_layer['layer'] not in (2, 4):
                #     continue
            if "approach" in each_layer:
                name += f"approach{each_layer['approach']}_"
            else:
                continue

            if "boundary" in each_layer:
                name += f"boundary{each_layer['boundary']}_"
                continue

            # if int(each_layer['boundary']) != 16:
            #     continue


            name = name[:-1]
            if self.check_layer(each_layer):
                if name in self.layers:
                    print(f"Layer {name} already exists")
                count += 1
                if True:
                    print(f"Creating layer {name}")
                    r = self.geometryR.getR(int(each_layer['volume']), int(each_layer['layer']), int(each_layer['approach']))
                    each_bin_collection = DetectorBinCollection(name, self.reg, self.world_logical, each_layer["value"]["material"], self.material_creator, self.solid_creator, r)
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

class CylinderR:
    def __init__(self, geometrymapfile):
        with open(geometrymapfile) as f:
            geometrymap = json.load(f)
        self.Rdic = {}
        geometrymap = geometrymap["Surfaces"]["entries"]

        for each in geometrymap:
            if "layer" not in each or "approach" not in each:
                continue
            if each["value"]["bounds"]["type"] != "CylinderBounds":
                continue
            key = self.getindex(int(each["volume"]), int(each["layer"]), int(each["approach"]))
            if key in self.Rdic:
                print(f"Error: CylinderR key {key} already exists in Rdic")
                exit(1)
            self.Rdic[key] = float(each["value"]["bounds"]["values"][0])
    
    def getR(self, volumeid, layerid, approachid):
        key = self.getindex(volumeid, layerid, approachid)
        if key not in self.Rdic:
            return -1
        return self.Rdic[key]

    def getindex(self, volumeid, layerid, approachid):
        return volumeid * 1000000 + layerid * 1000 + approachid

def main():
    print("Reading ITk gdml")
    itk = pyg4ometry.gdml.Reader("ITk.gdml")
    reg = itk.getRegistry()

    print("Loading ACTS material and geometry map")
    with open("material-map.json") as f:
        atlas_geometry = json.load(f)
    geometry_r = CylinderR("geometry-map.json") 
    atlas_geometry = atlas_geometry["Surfaces"]["entries"]
    atlas_detector = Detector("ATLAS", atlas_geometry, geometry_r)

    print("Resue world volume from ITk gdml")
    atlas_detector.use_world(reg.worldVolume)
    
    # print("Creating world volume")
    # atlas_detector.create_world()

    print("Generating detector passive layers")
    atlas_detector.generate_layers()

    print("Adding sensitive volumes from ITk gdml")
    for each_key in reg.physicalVolumeDict:
        physical = reg.physicalVolumeDict[each_key]
        logical = physical.logicalVolume
        material_name = logical.material.name
        if "SiliconMat" in material_name:
            atlas_detector.add_volume(physical)
    atlas_detector.save(f"ATLAS.gdml")

if __name__ == "__main__":
    main()
