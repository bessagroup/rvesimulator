
import re

import numpy as np


class rvesimulator2crateIO:

    def read_template(self, file_path) -> None:
        """Reads the file and returns its contents as a list of lines.

        file_path: str
            The path to the file.
        """
        with open(file_path, 'r') as f:
            self.lines = f.readlines()

    def replace_loadings(self,
                         strain: np.ndarray,
                         stress: np.ndarray = None,
                         prescription_index: np.ndarray = None):
        """Replace the Macroscale_Strain with the given strain. 

        strain: np.ndarray
            The strain to replace the Macroscale_Strain with. it has the shape
            (4, num_steps) where the first row is E11, the second row is E21,
            the third row is E12, and the fourth row is E22. The following rows
            should be replace by each row of the strain.
        stress: np.ndarray
            The stress to replace the Macroscale_Stress with. It has the shape
            (4, num_steps) where the first row is S11, the second row is S21,
            the third row is S12, and the fourth row is S22. The following rows
            should be replace by each row of the stress.
        prescription_index: np.ndarray
            The prescription index to replace the Macroscale_Prescription_Index.

        """
        # assert the shape of the strain, stress, and prescription_index, they
        # should have the same shape
        if stress is not None and prescription_index is not None:
            assert strain.shape == stress.shape == prescription_index.shape, \
                "Shape  strain, stress, and index should be the same."
        # assert the shape of the strain, it should have 4 rows
        assert strain.shape[0] == 4, "The strain should have 4 rows."

        # Get the line number of the Macroscale_Strain
        steps, line_num = self._extract_property("Macroscale_Strain")
        if not line_num:
            raise ValueError(
                "The Macroscale_Strain is not present in the input file.")
        # Replace the line with the new strain
        num_applied_steps = strain.shape[1]
        # replace this line with new number of steps
        self.lines[line_num-1] = f"Macroscale_Strain {num_applied_steps}\n"
        # replace the strain
        self.lines[line_num+0] = f"E11 " + \
            " ".join(["{:.8e}".format(x) for x in strain[0]]) + "\n"
        self.lines[line_num+1] = f"E21 " + \
            " ".join(["{:.8e}".format(x) for x in strain[1]]) + "\n"
        self.lines[line_num+2] = f"E12 " + \
            " ".join(["{:.8e}".format(x) for x in strain[2]]) + "\n"
        self.lines[line_num+3] = f"E22 " + \
            " ".join(["{:.8e}".format(x) for x in strain[3]]) + "\n"

        if stress is not None:
            steps, line_num = self._extract_property("Macroscale_Stress")
            if not line_num:
                raise ValueError(
                    "The Macroscale_Stress is not present in the input file.")
            # replace the line with the new stress
            num_applied_steps = stress.shape[1]
            # replace this line with new number of steps
            self.lines[line_num-1] = f"Macroscale_Stress {num_applied_steps}\n"
            # replace the stress
            self.lines[line_num+0] = f"S11 " + \
                " ".join(["{:.8e}".format(x) for x in stress[0]]) + "\n"
            self.lines[line_num+1] = f"S21 " + \
                " ".join(["{:.8e}".format(x) for x in stress[1]]) + "\n"
            self.lines[line_num+2] = f"S12 " + \
                " ".join(["{:.8e}".format(x) for x in stress[2]]) + "\n"
            self.lines[line_num+3] = f"S22 " + \
                " ".join(["{:.8e}".format(x) for x in stress[3]]) + "\n"
        else:
            stress = np.zeros_like(strain)
            steps, line_num = self._extract_property("Macroscale_Stress")
            if not line_num:
                raise ValueError(
                    "The Macroscale_Stress is not present in the input file.")
            # replace the line with the new stress
            num_applied_steps = stress.shape[1]
            # replace this line with new number of steps
            self.lines[line_num-1] = f"Macroscale_Stress {num_applied_steps}\n"
            # replace the stress
            self.lines[line_num+0] = f"S11 " + \
                " ".join(["{:.2e}".format(x) for x in stress[0]]) + "\n"
            self.lines[line_num+1] = f"S21 " + \
                " ".join(["{:.2e}".format(x) for x in stress[1]]) + "\n"
            self.lines[line_num+2] = f"S12 " + \
                " ".join(["{:.2e}".format(x) for x in stress[2]]) + "\n"
            self.lines[line_num+3] = f"S22 " + \
                " ".join(["{:.2e}".format(x) for x in stress[3]]) + "\n"

        # replace the prescription index
        # Get the line number of the Macroscale_Prescription_Index
        if prescription_index is not None:
            steps, line_num = self._extract_property(
                "Mixed_Prescription_Index")
            if not line_num:
                raise ValueError(
                    "The Mixed_Prescription_Index is not present in the input file.")
            # replace the line with the new prescription index
            num_applied_steps = prescription_index.shape[1]
            # replace this line with new number of steps
            self.lines[line_num -
                       1] = f"Mixed_Prescription_Index {num_applied_steps}\n"
            # replace the prescription index
            for i in range(prescription_index.shape[0]):
                self.lines[line_num+0] = " ".join(["{:.8e}".format(x)
                                                  for x in prescription_index[i]]) + "\n"
        else:
            prescription_index = np.zeros_like(strain)
            # set the values to integer
            prescription_index = prescription_index.astype(int)
            steps, line_num = self._extract_property(
                "Mixed_Prescription_Index")
            if not line_num:
                raise ValueError(
                    "The Mixed_Prescription_Index is not present in the input file.")
            # replace the line with the new prescription index
            num_applied_steps = prescription_index.shape[1]
            # replace this line with new number of steps
            self.lines[line_num -
                       1] = f"Mixed_Prescription_Index {num_applied_steps}\n"
            # replace the prescription index
            for i in range(prescription_index.shape[0]):
                # keep the integer format
                self.lines[line_num+i] = " ".join(["{:.0f}".format(x)
                                                  for x in prescription_index[i]]) + "\n"

    def replace_discretization_file(self, file_path):
        """Replace the Discretization_File with the given file_path.

        file_path: str
            The path to the new Discretization_File.
        """
        # Get the line number of the Discretization_File
        _, line_num = self._extract_property("Discretization_File")
        if not line_num:
            raise ValueError(
                "The Discretization_File is not present in the input file.")
        # Replace the line with the new file path
        self.lines[line_num] = f"{file_path}\n"

    def change_clustering_scheme(self,
                                 n_particle_clusters: int = 6,
                                 n_matrix_clusters: int = 12):
        """Change the number of particle and matrix clusters.

        n_particle_clusters: int
            The number of particle clusters.
        n_matrix_clusters: int
            The number of matrix clusters.
        """
        # Get the line number of the Clustering_Scheme
        _, line_num = self._extract_property("Number_of_Clusters")
        # change the number of clusters for material phase 1 (matrix)
        self.lines[line_num] = f"1 {n_matrix_clusters}\n"
        # change the number of clusters for material phase 2 (particle)
        self.lines[line_num+1] = f"2 {n_particle_clusters}\n"

    def change_matrix_material_properties(self,
                                          young_modulus: float,
                                          poisson_ratio: float,
                                          hardening_law_table: np.ndarray):
        """Change the matrix material properties.

        young_modulus: float
            The Young's modulus of the matrix.
        poisson_ratio: float
            The Poisson's ratio of the matrix.
        hardening_law_table: np.ndarray
            The hardening law table of the matrix.
        """
        # Get the line number of the Material_Phases
        num_mat, line_num = self._extract_property("Material_Phases")
        # change the young modulus and poisson ratio
        self.lines[line_num+2] = f"    E {young_modulus}\n"
        self.lines[line_num+3] = f"    v {poisson_ratio}\n"
        # change the hardening law table according to the new table
        self.lines[line_num +
                   4] = f"  isotropic_hardening piecewise_linear {hardening_law_table.shape[0]}\n"
        for i in range(hardening_law_table.shape[0]):
            self.lines[line_num+5+i] = f"    hp " + " ".join(
                ["{:.8e}".format(x) for x in hardening_law_table[i]]) + "\n"

    def change_particle_material_properties(self,
                                            young_modulus: float,
                                            poisson_ratio: float,
                                            ):
        """Change the particle material properties.

        young_modulus: float
            The Young's modulus of the particle.
        poisson_ratio: float
            The Poisson's ratio of the particle.
        """
        # Get the line number of the "2 elastic 1"
        _, line_num = self._extract_property("2 elastic 1")
        # change the young modulus and poisson ratio
        self.lines[line_num+1] = f"    E {young_modulus}\n"
        self.lines[line_num+2] = f"    v {poisson_ratio}\n"

    def change_size_of_rve(self, size: float):
        """Change the size of the RVE.

        size: np.ndarray
            The size of the RVE.
        """
        # Get the line number of the RVE_Size
        _, line_num = self._extract_property("RVE_Dimensions")
        # change the size of the RVE
        self.lines[line_num] = f"{size} {size}\n"

    def write_file(self, file_path):
        """Write the modified lines to a new file.

        file_path: str
            The path to the new file.
        """
        with open(file_path, 'w') as f:
            f.writelines(self.lines)

    def _extract_property(self, property_name):
        """Try to get the expected property from the input file.

        lines: list of strings
            The lines of the input file.
        property_name: str
            The name of the property to extract.
        """

        pattern = re.compile(f'^{property_name}\s*(.*)', re.IGNORECASE)

        for line_num, line in enumerate(self.lines, start=1):
            match = pattern.match(line)
            if match:
                return match.group(1).strip(), line_num
        return None, None
