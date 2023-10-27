"""Back of the napkin energetics calculator.

This calculator provides first-pass approximations of animal energy flow (kJ/day)
based on body-mass and temperature. 

"""
import tkinter as tk
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from math import exp, log10


class EnergyCalculator:
    """GUI-based energy calculator for forest ecosystems."""

    def __init__(self, master):
        """
        Initialize the calculator.

        Args:
            master (Tk object): Tkinter root object.
        """
        self.master = master
        self.master.title("Energy Calculator")
        self.master.geometry("800x800")

        # Create a Canvas widget
        self.canvas = tk.Canvas(self.master, width=800, height=800)
        self.scrollbar = tk.Scrollbar(
            self.master, orient="vertical", command=self.canvas.yview
        )
        self.scrollable_frame = tk.Frame(self.canvas)

        self.canvas.configure(yscrollcommand=self.scrollbar.set)

        self.scrollbar.pack(side="right", fill="y")
        self.canvas.pack(side="left", fill="both", expand=True)

        self.canvas.create_window((400, 0), window=self.scrollable_frame, anchor="n")

        self.scrollable_frame.bind(
            "<Configure>",
            lambda e: self.canvas.configure(scrollregion=self.canvas.bbox("all")),
        )

        # Initialize empty list to store animal group data
        self.animal_data = []

        # Create and place widgets
        self.label = tk.Label(
            self.scrollable_frame, text="Forest Energetics Calculator"
        )
        self.label.pack(side="top", fill="both", expand="yes")

        self.lbl_result = tk.Label(self.scrollable_frame, text="")
        self.lbl_result.pack()

        # Entry box for temperature
        self.temperature_entry = tk.Entry(self.scrollable_frame)
        self.temperature_entry.pack()
        self.temperature_entry.insert(0, "Temperature (°C)")

        self.add_button = tk.Button(
            self.scrollable_frame, text="Add Animal Group", command=self.add_group
        )
        self.add_button.pack()

        self.calculate_button = tk.Button(
            self.scrollable_frame, text="Calculate", command=self.calculate_energy
        )
        self.calculate_button.pack()

        self.stats_button = tk.Button(
            self.scrollable_frame, text="Display Statistics", command=self.display_stats
        )
        self.stats_button.pack()

        self.canvas_frame = tk.Frame(self.scrollable_frame)
        self.canvas_frame.pack()

    def add_group(self):
        """Add a new animal group with its respective parameters in the GUI."""

        # Frame to hold widgets for this animal group
        group_frame = tk.Frame(self.scrollable_frame)
        group_frame.pack(pady=10)

        # Dropdown menu for animal group name
        OPTIONS = ["Mammals", "Birds", "Insects", "Frogs"]
        group_name_var = tk.StringVar(self.scrollable_frame)
        group_name_var.set(OPTIONS[0])
        group_menu = tk.OptionMenu(group_frame, group_name_var, *OPTIONS)
        group_menu.grid(row=0, column=0)

        # Create a StringVar to hold the type of the measurement
        measurement_type = tk.StringVar()
        measurement_type.set("Mass (kg)")

        # Entry box for body mass
        mass_entry = tk.Entry(group_frame)
        mass_entry.grid(row=0, column=1)
        mass_entry.insert(0, "Body Mass (kg)")

        # Entry box for number of individuals
        number_entry = tk.Entry(group_frame)
        number_entry.grid(row=0, column=2)
        number_entry.insert(0, "Number")

        METABOLIC_OPTIONS = ["Endothermic", "Ectothermic"]
        metabolic_var = tk.StringVar(self.scrollable_frame)
        metabolic_var.set(METABOLIC_OPTIONS[0])
        metabolic_menu = tk.OptionMenu(group_frame, metabolic_var, *METABOLIC_OPTIONS)
        metabolic_menu.grid(row=0, column=4)

        # Button to remove this animal group
        remove_button = tk.Button(
            group_frame, text="Remove", command=lambda: group_frame.destroy()
        )
        remove_button.grid(row=0, column=3)

        # Create the Entry box for mass/length
        measure_entry = tk.Entry(group_frame)
        measure_entry.grid(row=0, column=1)
        measure_entry.insert(0, measurement_type.get())

        # Update the label and clear the entry when animal type changes
        def update_measurement(*args):
            animal_type = group_name_var.get()
            if animal_type == "Insects":
                measurement_type.set("Length (cm)")
            else:
                measurement_type.set("Mass (kg)")
            measure_entry.delete(0, "end")
            measure_entry.insert(0, measurement_type.get())

        # Listen for changes to the animal type dropdown
        group_name_var.trace("w", update_measurement)

        # ... (rest of the code for add_group)
        self.animal_data.append(
            {
                "name": group_name_var,
                "measure": measure_entry,
                "number": number_entry,
                "metabolic_type": metabolic_var,
                "measurement_type": measurement_type,
            }
        )

    def calculate_energy(self):
        """Calculate the total energy consumed by all animal groups and update the result label."""
        total_energy = 0.0
        group_energies = []  # To store specific group energies

        # Fetch temperature from the widget
        temperature = float(self.temperature_entry.get())

        for group in self.animal_data:
            group_energy = self.calculate_group_energy(
                group["name"].get(),
                float(group["measure"].get()),
                int(group["number"].get()),
                group["metabolic_type"].get().lower(),
                group["measurement_type"].get(),
                temperature,
            )
            total_energy += group_energy
            group_energies.append(f"{group['name'].get()}: {group_energy} kJ/day")

        specific_energies = "\n".join(group_energies)
        self.lbl_result.config(
            text=f"Total Energy: {total_energy} kJ/day\nSpecific Group Energies:\n{specific_energies}"
        )

    def metabolic_rate(self, mass, metabolic_type, temperature=25.0):
        """
        Calculate the metabolic rate based on mass, metabolic_type type, and temperature.

        Args:
            mass (float): Body mass in kg.
            metabolic_type (str): 'Endothermic' or 'Ectothermic'.
            temperature (float, optional): Environmental temperature in degrees Celsius. Defaults to 25.0.

        Returns:
            float: Metabolic rate.
        """

        # Es = 3.7 * 10 ** (-2)  # energy to mass conversion constant (g/kJ)
        sig = 0.5  # proportion of time-step with temp in active range (toy)
        Ea = 0.69  # aggregate activation energy of metabolic reactions
        kB = 8.617333262145e-5  # eV/K
        mass_g = mass * 1000  # convert mass to grams

        if metabolic_type == "endothermic":
            Ib, bf = (4.19 * 10**10, 0.69)  # field metabolic constant and exponent
            If, bb = (9.08 * 10**11, 0.7)  # basal metabolic constant and exponent
            Tk = 310.0  # body temperature of the individual (K)
            return (sig * If * exp(-(Ea / (kB * Tk)))) * mass_g**bf + (
                (1 - sig) * Ib * exp(-(Ea / (kB * Tk)))
            ) * mass_g**bb
        elif metabolic_type == "ectothermic":
            Ib, bf = (4.19 * 10**10, 0.69)  # field metabolic constant and exponent
            If, bb = (1.49 * 10**11, 0.88)  # basal metabolic constant and exponent
            Tk = temperature + 274.15  # body temperature of the individual (K)
            return (sig * If * exp(-(Ea / (kB * Tk)))) * mass_g**bf + (
                (1 - sig) * Ib * exp(-(Ea / (kB * Tk)))
            ) * mass_g**bb

    def calculate_group_energy(
        self, name, measure, number, metabolic_type, measurement_type, temperature
    ):
        """
        Calculate the total energy consumption of a particular animal group.

        Args:
            name (str): Name of the animal group.
            measure (float): Measurement (could be mass or length).
            number (int): Number of individuals in the group.
            metabolic_type (str): 'Endothermic' or 'Ectothermic'.
            measurement_type (str): 'Mass (kg)' or 'Length (cm)'.
            temperature (float): Environmental temperature in degrees Celsius.

        Returns:
            float: Total energy consumption for the group.
        """
        if measurement_type == "Length (cm)":
            # Use your allometric equation here to get mass from length
            mass = self.length_to_mass(measure)
        else:
            mass = measure

        metabolic_rate = self.metabolic_rate(mass, metabolic_type, temperature)
        return metabolic_rate * number  # We may adjust this formula as necessary

    def length_to_mass(self, length):
        """
        Convert an insect's length to mass using allometric scaling.

        The formula used is:
        mass = 10^(intercept + slope * log10(length))

        Parameters:
        length (float): The length of the insect in cm.

        Returns:
        mass: The estimated mass of the insect in kg.

        Equation assumes length in mm and mass in mg, conversions are applied.

        Notes:
        - The intercept and slope from Sohlström et al 2018
        """
        intercept = -0.792
        slope = 2.181
        log_length = log10(length * 10)  # Convert cm to mm for the equation
        log_mass = intercept + slope * log_length
        mass = 10**log_mass * 1e-6  # Convert mg to kg for the output
        return mass

    def display_stats(self):
        """
        Display statistics on energy consumption in a Matplotlib plot embedded in the Tkinter window.

        The length to mass transformation is handled in calculate_group_energy

        """

        # Create a new figure and axis for Matplotlib
        fig, ax = plt.subplots()

        # Create lists to store names and total energies for plotting
        group_names = []
        group_energies = []

        temperature = float(
            self.temperature_entry.get()
        )  # Fetch temperature from the widget

        for group in self.animal_data:
            name = group["name"].get()
            measure = float(group["measure"].get())
            number = int(group["number"].get())
            metabolic_type = group["metabolic_type"].get().lower()
            measurement_type = group["measurement_type"].get()

            # Reusing calculate_group_energy function
            total_energy = self.calculate_group_energy(
                name, measure, number, metabolic_type, measurement_type, temperature
            )

            group_names.append(name)
            group_energies.append(total_energy)

            # Individual bar for this group
            ax.bar(name, total_energy)

        # Basic statistics
        mean_energy = np.mean(group_energies)
        median_energy = np.median(group_energies)
        std_dev_energy = np.std(group_energies)

        stats_text = f"Mean Energy: {mean_energy}\nMedian Energy: {median_energy}\nStd Dev: {std_dev_energy}"
        self.lbl_result.config(text=stats_text)

        # Add labels and title
        ax.set_ylabel("Total Energy (kJ/day)")
        ax.set_title("Total energy by animal group")

        # Display the Matplotlib plot on the Tkinter window
        canvas = FigureCanvasTkAgg(fig, master=self.canvas_frame)
        canvas.draw()
        canvas.get_tk_widget().pack()


if __name__ == "__main__":
    """Create and run the EnergyCalculator GUI."""
    root = tk.Tk()
    calculator = EnergyCalculator(root)
    root.mainloop()
