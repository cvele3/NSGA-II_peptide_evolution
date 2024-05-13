import RandomGenerator
import FitnessFunctionScraper


class Peptide:
    def __init__(self, peptide_list, peptide_string, ff_amp_probability):
        self.peptide_list = peptide_list
        self.peptide_string = peptide_string
        self.ff_amp_probability = ff_amp_probability

    
    def calculate(lower_length, upper_length, number_of_peptides):
        peptides = RandomGenerator.generate_random_peptides(lower_length, upper_length, number_of_peptides)

        list_of_peptide_objects = []

        with open('in.txt', 'w') as file:
            for peptide in peptides:
                peptide_string = ''.join(peptide)
                file.write(f'>{peptide_string}\n{peptide_string}\n')


        peptide_and_ff_amp_probability = FitnessFunctionScraper.scrape_fitness_function()

        for peptide, ff_amp_probability in peptide_and_ff_amp_probability:
            peptide_list = list(peptide)
            peptide_string = peptide
            ff_amp_probability = ff_amp_probability
            print(peptide_list, peptide_string, ff_amp_probability)
            list_of_peptide_objects.append(Peptide(peptide_list, peptide_string, ff_amp_probability))
            

