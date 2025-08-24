import pyg4ometry.geant4 as _g4
import periodictable
import pyg4ometry

class MaterialCreator:
    def __init__(self, reg, tolerance=0.04):
        self.materials = {}
        self.tolerance = tolerance
        self.reg = reg
        self.useapproximate = True

    def create_material(self, ZAFstring, massdensity):
        elements_string = ZAFstring.split(";")
        elements = []
        for each_element_string in elements_string:
            if len(each_element_string) == 0:
                continue
            each_element_string = each_element_string.replace("Z", "")
            each_z = float(each_element_string.split("A")[0])
            each_a = float(each_element_string.split("A")[1].split("F")[0])
            each_f = float(each_element_string.split("F")[1]) # molar fraction
            elements.append((each_z, each_a, each_f))
        elements.sort(key=lambda x: x[0])
        molarmass = sum([each[1] * each[2] for each in elements])
        massdensity = massdensity # g/cm3

        materialkey = ""
        for each in elements:
            materialkey += str(int(each[0])) + "_"

        if materialkey in self.materials and self.useapproximate:
            for each_material in self.materials[materialkey]:
                # compare mass density
                if abs(each_material[2] - massdensity)/massdensity > 0.05:
                    continue
                # compare fractions
                passed = True
                for i in range(len(each_material[1])):
                    if abs(each_material[1][i] - elements[i][2]) > self.tolerance:
                        passed = False
                        break
                if passed:
                    # print("Using existing material: ", each_material[0].name)
                    return each_material[0]

        if materialkey not in self.materials:
            self.materials[materialkey] = []
        this_id = len(self.materials[materialkey])
        # print("Creating new material: ", "material_" + materialkey + "id_" + str(this_id))
        this_material = _g4.MaterialCompound("material_" + materialkey + "id_" + str(this_id), massdensity, len(elements), self.reg)
        this_fractions = []
        for i in range(len(elements)):
            each_z, each_a, each_f = elements[i]
            each_symble = "NONE"
            each_name = "NONE"
            for el in periodictable.elements:
                if el.number == int(each_z):
                    each_symble = el.symbol
                    each_name = el.name
                    break
            element = _g4.ElementSimple(each_name, each_symble, int(each_z), each_a)
            this_material.add_element_massfraction(element, each_f*each_a/molarmass)
            this_fractions.append(each_f)
        self.materials[materialkey].append((this_material, this_fractions, massdensity))
        return this_material


if __name__ == "__main__":
    reg = pyg4ometry.geant4.Registry()
    mat = MaterialCreator(reg)
    mat.create_material("Z1A1.00794F0.043;Z6A12.0107F0.793;Z7A14.0068F0.007;Z8A15.9994F0.121;Z14A28.0854F0.036;", 8.96)
    mat.create_material("Z1A1.00794F0.043;Z6A12.0107F0.793;Z7A14.0068F0.007;Z8A15.9994F0.121;", 8.96)
    mat.create_material("Z1A1.00794F0.043;Z6A12.0107F0.793;Z7A14.0068F0.007;Z8A15.9994F0.121;Z14A28.0854F0.033;", 8.8)
